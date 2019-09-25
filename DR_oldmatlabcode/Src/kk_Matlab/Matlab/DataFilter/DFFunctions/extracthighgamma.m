function extracthighgamma(directoryname, fileprefix, day, tetrode, nstd)

%function extracthighgamma(directoryname, fileprefix, day, tetrode, min_suprathresh_duration, nstd)
% Reads in the high gamma files from the specified day and tetrodes and
% extracts all of the large gamma events from those tetrodes.
% 
% High gamma events are defined as 400ms windows around low gamma amplitude
% maxima. Gamma events must occur during valid theta times and must occur
% at least 100ms apart from each other. This extraction is based on the
% methods described in Colgin LL et al., Nature 2009.
%
%directoryname - example '/data99/user/animaldatafolder/', a folder
%                containing processed matlab data for the animal
%
%fileprefix	- folder name where the day's data is stored
%
%day		- the day to process
%
%tetrode_list	- the tetrode(s) to process.
%			-1 indicates all tetrodes should be processed
%           0 indicates that the tetrode with the most cells will
%               be processed
%			#(s) indicate the set of tetrodes to proces
%           a string describing which tetrodes in tetinfo to use
%               Ex: 'isequal($area,''CA1'')'
%
%nstd		- the number of standard dev that gamma must be from mean to
%			be detected. Start with 2.
%
% Outputs:
%gamma 	- structue with various fields, including the following which
%			describe each event.
%	starttime - time of beginning of gamma
%	endtime	  - time of end of gamma
%	midtime   - time of midpoint of energy of event
%	peak	  - peak height of waveform)
%	maxthresh - the largest threshold in stdev units at which this gamma
%			would still be detected.
%	startind  - index of start time in gamma structure
%	endind    - index of end time in gamma structure
%	midind    - index of middle time in gamma structure
%	posind    - index into pos structure for the midpoint of this gamma
% 	posinterp - interpolation factor value between adjacent position
% 		    elements


% define the standard deviation for the Gaussian smoother which we
% apply before thresholding (this reduces sensitivity to spurious
% flucutations in the gamma envelope)

smoothing_width = 0.01; % 10 ms

if tetrode == 0
    %select eeg with most cells
    numcells = findnumcell(directoryname, fileprefix, day, 2);
    [i tet] = max(numcells) ;
elseif isstr(tetrode)
    load(sprintf('%s/%stetinfo.mat',directoryname,fileprefix));
    tmptet = evaluatefilter(tetinfo,tetrode);
    tet = unique(tmptet(tmptet(:,1)==day,3));
    clear tetinfo
elseif tetrode == -1
    tet = 1:1000;
else
    tet = tetrode;
end

clear tetrode

dsz = '';
if (day < 10)
    dsz = '0';
end

% load the positionf file
eval(['load ', directoryname, '/', fileprefix, 'pos', dsz, num2str(day), '.mat']); %load pos file
load '/home/mcarr/Src/Matlab/Filters/thetafilter.mat'

% move to the EEG directory
pushd([directoryname,'/EEG']);

% go through each tetrode
for j = 1:length(tet)
    t = tet(j);
    tmpflist = dir(sprintf('*highgamma%02d-*-%02d.mat', day, t)); 
    tmptheta = dir(sprintf('*eeg%02d-*-%02d.mat', day, t));
    
    % Compute the baseline and day for the whole day
    allrenv = [];
    for i = 1:length(tmpflist)
    	load(tmpflist(i).name);
        % get the epoch number
        dash = find(tmpflist(i).name == '-');
        e = str2num(tmpflist(i).name((dash(1)+1):(dash(2)-1)));
            
        % convert the gamma envelope field to double
        temprenv = double(highgamma{day}{e}{t}.data(:,3));
        % smooth the envelope:
        samprate = highgamma{day}{e}{t}.samprate;
        kernel = gaussian(smoothing_width*samprate, ceil(8*smoothing_width*samprate));
        temprenv = smoothvect(temprenv, kernel);
        allrenv = [allrenv; temprenv];
    end
    totalbaseline = mean(allrenv);
    totalstdev = std(allrenv);
    clear allrenv clear temprenv
    
    for i = 1:length(tmpflist)
       % load the lowgamma and theta files
        load(tmpflist(i).name);
        
        % get the epoch number
        dash = find(tmpflist(i).name == '-');
        e = str2num(tmpflist(i).name((dash(1)+1):(dash(2)-1)));
        for k = 1:length(tmptheta)
            dash = find(tmptheta(k).name == '-');
            te = str2num(tmptheta(k).name((dash(1)+1):(dash(2)-1)));
            if te == e
                load(tmptheta(k).name)
            end
        end
        
        % convert the gamma envelope field to double
        gamma = double(highgamma{day}{e}{t}.data(:,3));
        time = 1:length(gamma);
        
        % smooth the envelope:
        samprate = highgamma{day}{e}{t}.samprate;
        kernel = gaussian(smoothing_width*samprate, ceil(8*smoothing_width*samprate));
        renv = smoothvect(gamma, kernel);
        clear gamma
        
        % calculate the threshold in uV units
        baseline = mean(renv);
        stdev = std(renv);
        thresh = baseline + nstd * stdev;       
       
        % find times that exceed threshold
        tmpind = time(renv > thresh);
        
        % create 160ms windows around each time that exceeds threshold
        win = round(0.160 * samprate);
        while ~isempty(tmpind) && tmpind(1) <win+5
            tmpind(1) = [];
        end
        while ~isempty(tmpind) && tmpind(end)>length(renv)-win-5
            tmpind(end) = [];
        end
        tmpmax = zeros(length(tmpind),1);
         for k = 1:length(tmpind)
            tmpmax(k) = max(renv(tmpind(k)-win:tmpind(k)+win));
            tmpind(k) = tmpind(k) - win -1 + find(tmpmax(k) == renv(tmpind(k)-win:tmpind(k)+win));
        end
        
        % get rid of duplicate max values
        valid = diff([tmpmax; 0]) ~= 0;
        tmpmax = tmpmax(valid);
        tmpind = tmpind(valid);
        
        % require that maxima be separated by 100ms
        win = round(0.100 * samprate);
        validind = diff([tmpind 0]) > win;
        tmpind = tmpind(validind);        
        tmpmax = tmpmax(validind);
        
        
        %%% COMMENTED OUT BY KK 3.14.12
        
%         % require that maxima occur during valid theta times
%         % calculate valid theta times
%         phase = angle(hilbert(filtfilt(thetafilter,1,eeg{day}{e}{t}.data')));
%         ttime = geteegtimes(eeg{day}{e}{t});
%         validtheta = getvalideegtimes(phase,ttime);
%         valid = validtheta(tmpind);
%         tmpind = tmpind(logical(valid));
%         tmpmax = tmpmax(logical(valid));
%         clear eeg phase ttime validtheta valid
                
        if ~isempty(tmpmax)

            % make sure that window does not include the edges
            win = round(0.200 * samprate);
            while tmpind(1)-win < 0
                tmpind(1) = [];
                tmpmax(1) = [];
            end
            while tmpind(end)+win > length(time)
                tmpind(end) = [];
                tmpmax(end) = [];
            end
            % create 400ms windows around each valid index
            gamma.startind = tmpind -win;
            gamma.endind = tmpind + win;
            gamma.midind = tmpind;
            gamma.starttime = highgamma{day}{e}{t}.starttime + gamma.startind / samprate;
            gamma.endtime = highgamma{day}{e}{t}.starttime + gamma.endind / samprate;
            gamma.midtime = highgamma{day}{e}{t}.starttime + gamma.midind / samprate;
            gamma.peak = tmpmax;
            gamma.maxthresh = (tmpmax - baseline) / stdev;

            % position time element
            postime = pos{day}{e}.data(:,1);

            gamma.posind = floor(interp1q(postime,(1:length(postime))',gamma.midtime'));
        else
            gamma.startind = [];
            gamma.endind = [];
            gamma.midind = [];
            gamma.starttime = [];
            gamma.endtime = [];
            gamma.midtime = [];
            gamma.peak = [];
            gamma.energy = [];
            gamma.posind = [];
        end
        
        gamma.timerange = [0 length(renv)/samprate] + highgamma{day}{e}{t}.starttime;
        gamma.samprate = highgamma{day}{e}{t}.samprate;
        gamma.threshold = thresh;
        gamma.baseline = baseline;
        gamma.std = stdev;
        gamma.totalbaseline = totalbaseline;
        gamma.totalstd = totalstdev;
        
        %gamma
        clear feeg
        gammah{day}{e}{t} = gamma;
        clear highgamma;
        clear gamma
    end
end
save(sprintf('%s/%sgammah%02d.mat', directoryname, fileprefix, day), 'gammah');

popd;

