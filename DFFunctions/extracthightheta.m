function extracthightheta(directoryname, fileprefix, day, tetrode, min_suprathresh_duration, nstd, maxpeak)

%function extracthightheta(directoryname, fileprefix, day, tetrode_list, 
%                min_suprathresh_duration, nstd, maxpeakval)
%
%	Reads in the theta files from the specified day and tetrodes and
%	extracts all of the high theta periods from that tetrodes.
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
%       must remain above threshold to be counted as as theta; this guards
%       against short spikes in signal (e.g. noise) as being counted as
%       hightheta. Set min_suprathreshold_duration to some small value, like
%       0.015 s.
%
%nstd		- the number of standard dev that theta must be from mean to 
%			be detected. Start with 2.
%
%maxpeakval	- hightheta events with maximal peaks above this value are exluded.  Use
%			this avoid detecting noise events. Start with 1000
% Outputs:
%hightheta 	- structue with various fields, including the following which
%			describe each theta. 
%	starttime - time of beginning of theta
%	endtime	  - time of end of theta
%	midtime   - time of midpoint of energy of event
%	peak	  - peak height of waveform)
%	maxthresh - the largest threshold in stdev units at which this theta 
%			would still be detected.
%	energy	  - total sum squared energy of waveform
%	startind  - index of start time in theta structure
%	endind    - index of end time in theta structure
%	midind    - index of middle time in theta structure
%	posind    - index into pos structure for the midpoint of this theta
% 	posinterp - interpolation factor value between adjacent position
% 		    elements
%    
%	


% define the standard deviation for the Gaussian smoother which we
% apply before thresholding (this reduces sensitivity to spurious
% flucutations in the theta envelope)
smoothing_width = 1; % 1 second

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
currdir = pwd;
cd(fullfile(directoryname,'/EEG'));
if (tet == -1)
    % set up to handle all tetrodes 
    tet = 1:1000;
end
% go through each tetrode 
for t = tet
    tmpflist = dir(sprintf('*theta%02d-*-%02d.mat', day, t));
    for i = 1:length(tmpflist)
        % load the theta file
        load(tmpflist(i).name);
        % get the epoch number
	dash = find(tmpflist(i).name == '-');
	e = str2num(tmpflist(i).name((dash(1)+1):(dash(2)-1)));

	% convert the theta envelope field to double
   if size(theta{d}{e}{t}.data,2) == 2
      renv = double(theta{d}{e}{t}.data(:,2));
   else
      renv = double(theta{d}{e}{t}.data(:,3));
   end

	% smooth the envelope:
	samprate = theta{d}{e}{t}.samprate;
	kernel = gaussian(smoothing_width*samprate, ceil(8*smoothing_width*samprate));
	renv = smoothvect(renv, kernel);
	% find the hightheta
	% calculate the duration in terms of samples
	mindur = round(min_suprathresh_duration * samprate);

	% calculate the threshold in uV units 
	baseline = mean(renv);
	stdev = std(renv);
	thresh = baseline + nstd * stdev;

	% extract the events if this is a valid trace
	if (thresh > 0) 
	    tmphightheta = extractevents(renv, thresh, baseline, 0, mindur, 0)';
	    % Assign the fields 
	    % start and end indeces
	    htheta.startind = tmphightheta(:,1);
	    htheta.endind = tmphightheta(:,2);
	    % middle of energy index
	    htheta.midind = tmphightheta(:,8);

	    %convert the samples to times for the first three fields
	    htheta.starttime = theta{d}{e}{t}.starttime + htheta.startind / samprate;
	    htheta.endtime = theta{d}{e}{t}.starttime + htheta.endind / samprate;
	    htheta.midtime = theta{d}{e}{t}.starttime + htheta.midind / samprate;
	    htheta.peak = tmphightheta(:,3);
	    htheta.energy = tmphightheta(:,7);
	    htheta.maxthresh = tmphightheta(:,9) / stdev;

	    % position time element
       if ~isempty(pos{d}{e}.data)
          postime = pos{d}{e}.data(:,1);
          idxs = (1:length(postime))';
          htheta.posind = floor(interp1q(postime,idxs,htheta.midtime));
       else
          htheta.posind = [];
       end
	else
	    htheta.startind = [];
	    htheta.endind = [];
	    htheta.midind = [];
	    htheta.starttime = [];
	    htheta.endtime = [];
	    htheta.midtime = [];
	    htheta.peak = [];
	    htheta.energy = [];
	    htheta.posind = [];
	    htheta.maxthresh = [];
	end

	htheta.timerange = [0 length(renv)/samprate] + theta{d}{e}{t}.starttime;
	htheta.samprate = theta{d}{e}{t}.samprate;
	htheta.threshold = thresh;
	htheta.baseline = baseline;
	htheta.std = std(renv);
	htheta.minimum_duration = min_suprathresh_duration;

	%theta
	clear feeg
	hightheta{d}{e}{t} = htheta;
	clear htheta;
    end
end
save(sprintf('%s/%shightheta%02d.mat', directoryname, fileprefix, d), 'hightheta');
cd(currdir);
