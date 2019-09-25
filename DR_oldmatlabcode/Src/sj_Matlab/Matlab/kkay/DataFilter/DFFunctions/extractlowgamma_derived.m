function extractlowgamma_derive(directoryname, fileprefix, day, tetrode, min_suprathresh_duration, nstd, varargin)

%% simple modification of common shared extractripple -- 6.19.12

%function extractripples(directoryname, fileprefix, day, tetrode_list, %%min_suprathresh_duration, nstd, maxpeakval, options)
%
%	Reads in the ripple files from the specified day and tetrodes and
%	extracts all of the ripples from that tetrodes.
%
%	assumes position data stores in pos file in animdirectory
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
%			 0 indicates that the tetrode with the most cells will
%			 	be processed
%			 #(s) indicate the set of tetrodes to proces
%min_suprathresh_duration
%		- the time (in seconds) which the signal
%       must remain above threshold to be counted as as ripple; this guards
%       against short spikes in signal (e.g. noise) as being counted as
%       ripples. Set min_suprathreshold_duration to some small value, like
%       0.015 s.
%
%nstd		- the number of standard dev that ripple must be from mean to
%			be detected. Start with 2.
%
%
%options	'stdev', stdev   sets the size of the standard deviation used to
%				allow direct comparison across epochs
%       	'baseline', b   sets the size of the baseline used to
%				allow direct comparison across epochs
%           'maxpeakval, m	- ripples with maximal peaks above this value
%				are exluded.  Use this avoid detecting noise
%				events. Default 1000
%           'samethreshperday' - 0 or 1, default 0, 
%               0 calculates baseline and threshold per epoch for each
%               tetrode
%               1 uses same threshold for all sessions on the same day on 
%               the same tetrode.  Calculates this threshold based on the 
%               baseline and stdev for all sessions (run and sleep) for
%               the entire day.    AS added 1-12-10
%
% Outputs:
%ripples 	- structue with various fields, including the following which
%			describe each ripple.
%	starttime - time of beginning of ripple
%	endtime	  - time of end of ripple
%	midtime   - time of midpoint of energy of event
%	peak	  - peak height of waveform)
%	maxthresh - the largest threshold in stdev units at which this ripple
%			would still be detected.
%	energy	  - total sum squared energy of waveform
%	startind  - index of start time in ripple structure
%	endind    - index of end time in ripple structure
%	midind    - index of middle time in ripple structure
%	posind    - index into pos structure for the midpoint of this ripple
% 	posinterp - interpolation factor value between adjacent position
% 		    elements
%
%

stdev = 0;
baseline = 0;
maxpeakval = 1000;
samethreshperday = 0;
for option = 1:2:length(varargin)-1
    if isstr(varargin{option})
        switch(varargin{option})
            case 'stdev'
                stdev = varargin{option+1};
            case 'baseline'
                baseline = varargin{option+1};
            case 'maxpeakval'
                maxpeakval = varargin{option+1};
            case 'samethreshperday'
                samethreshperday = varargin{option+1};
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end
    else
        error('Options must be strings, followed by the variable');
    end
end




% define the standard deviation for the Gaussian smoother which we
% apply before thresholding (this reduces sensitivity to spurious
% flucutations in the ripple envelope)
smoothing_width = 0.004; % 4 ms

d = day;

if tetrode == 0
    %select eeg with most cells
    numcells = findnumcell(directoryname, fileprefix, d, 2);
    [i tet] = max(numcells) ;
else
    tet = tetrode;
end

dsz = '';
if (d < 10)
    dsz = '0';
end

% load the positionf file
eval(['load ', directoryname, '/', fileprefix, 'pos', dsz, num2str(d), '.mat']); %load pos file



% move to the EEG directory
cd([directoryname,'/EEG']);
if (tet == -1)
    % set up to handle all tetrodes
    tet = 1:1000;
end
% go through each tetrode
for t = tet
    tmpflist = dir(sprintf('*lowgamma%02d-*-%02d.mat', day, t));
    
    %AS added following 1-12-10 to allow all epochs to sue same threshold
    if samethreshperday == 1
        allrenv = [];
        for i = 1:length(tmpflist)
            load(tmpflist(i).name);
            % get the epoch number
            dash = find(tmpflist(i).name == '-');
            e = str2num(tmpflist(i).name((dash(1)+1):(dash(2)-1)));
            
            % convert the ripple envelope field to double
            temprenv = double(lowgamma{d}{e}{t}.data(:,3));
            allrenv = [allrenv; temprenv];
        end
        baseline = mean(allrenv);
        stdev = std(allrenv);
        thresh = baseline + nstd * stdev;
    end
    
    
    for i = 1:length(tmpflist)
        % load the ripple file
        load(tmpflist(i).name);
        % get the epoch number
        dash = find(tmpflist(i).name == '-');
        e = str2num(tmpflist(i).name((dash(1)+1):(dash(2)-1)));
        
        % convert the ripple envelope field to double
        renv = double(lowgamma{d}{e}{t}.data(:,3));
        
        % smooth the envelope:
        samprate = lowgamma{d}{e}{t}.samprate;
        kernel = gaussian(smoothing_width*samprate, ceil(8*smoothing_width*samprate));
        renv = smoothvect(renv, kernel);
        % find the ripples
        % calculate the duration in terms of samples
        mindur = round(min_suprathresh_duration * samprate);
        
        % calculate the threshold in uV units
        if samethreshperday == 0 % is calculating threshold per epoch
            baseline = mean(renv);
            stdev = std(renv);
            thresh = baseline + nstd * stdev;
        end
        
        % AS commented out below and replace with 3 lines above.  1-12-10
        % % calculate the threshold in uV units
        %         if (stdev == 0)
        %             stdev = std(renv);
        %         end
        %         if (baseline == 0)
        %             baseline = mean(renv);
        %         end
        %         thresh = baseline + nstd * stdev;
        
        
        % extract the events if this is a valid trace
        if (thresh > 0) & any(find(renv<baseline))
            tmprip = extractevents(renv, thresh, baseline, 0, mindur, 0)';
            % Assign the fields
            % start and end indeces
            rip.startind = tmprip(:,1);
            rip.endind = tmprip(:,2);
            % middle of energy index
            rip.midind = tmprip(:,8);
            
            %convert the samples to times for the first three fields
            rip.starttime = lowgamma{d}{e}{t}.starttime + rip.startind / samprate;
            rip.endtime = lowgamma{d}{e}{t}.starttime + rip.endind / samprate;
            rip.midtime = lowgamma{d}{e}{t}.starttime + rip.midind / samprate;
            rip.peak = tmprip(:,3);
            rip.energy = tmprip(:,7);
            rip.maxthresh = (tmprip(:,9) - baseline) / stdev;
            
            % position time element
            postime = pos{d}{e}.data(:,1);
            
            idxs = (1:length(postime))';
            rip.posind = floor(interp1q(postime,idxs,rip.midtime));
        else
            rip.startind = [];
            rip.endind = [];
            rip.midind = [];
            rip.starttime = [];
            rip.endtime = [];
            rip.midtime = [];
            rip.peak = [];
            rip.energy = [];
            rip.posind = [];
        end
        if any(find(renv<baseline))==0
            warning(['No below baseline values in data.  Fields left blank, ',tmpflist(i).name])
        end
        
        rip.timerange = [0 length(renv)/samprate] + lowgamma{d}{e}{t}.starttime;
        rip.samprate = lowgamma{d}{e}{t}.samprate;
        rip.threshold = thresh;
        rip.baseline = baseline;
        rip.std = stdev;
        rip.minimum_duration = min_suprathresh_duration;
        
        %ripple
        clear feeg
        gammal{d}{e}{t} = rip;
        clear rip;
    end
end
save(sprintf('%s/%slowgamma%02d.mat', directoryname, fileprefix, d), 'gammal');
end
