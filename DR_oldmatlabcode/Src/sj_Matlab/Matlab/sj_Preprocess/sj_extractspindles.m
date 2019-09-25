function sj_extractspindles(directoryname, fileprefix, day, tetrode, min_suprathresh_duration, nstd, varargin)
% Shantanu - Nov 2012. From extractripples. Following Matt Jones, Neuron 2012 and Siapas

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
%       0.1 s.
%
%nstd		- the number of standard dev that ripple must be from mean to
%			be detected. Start with 2.
%
%
%options	'stdev', stdev   sets the size of the standard deviation used to
%				allow direct comparison across epochs
%       	'baseline', b   sets the size of the baseline used to
%				allow direct comparison across epochs
%               'maxpeakval, m	- ripples with maximal peaks above this value
%				are exluded.  Use this avoid detecting noise
%				events. Default 1000
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
for option = 1:2:length(varargin)-1
    if isstr(varargin{option})
        switch(varargin{option})
            case 'stdev'
                stdev = varargin{option+1};
            case 'baseline'
                baseline = varargin{option+1};
            case 'maxpeakval'
                maxpeakval = varargin{option+1};
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end
    else
        error('Options must be strings, followed by the variable');
    end
end

if (directoryname(end) == '/')
    directoryname = directoryname(1:end-1);
end


% define the standard deviation for the Gaussian smoother which we
% apply before thresholding (this reduces sensitivity to spurious
% flucutations in the ripple envelope)

% WHAT SHOULD THIS BE - Talk with loren. Say 20 ms for now
%smoothing_width = 0.004; % 4 ms
smoothing_width = 0.02; % 20 ms

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
pushd([directoryname,'/EEG']);
if (tet == -1)
    % set up to handle all tetrodes
    tet = 1:1000;
end
% go through each tetrode
for t = tet
    tmpflist = dir(sprintf('*spindlegnd%02d-*-%02d.mat', day, t));
    for i = 1:length(tmpflist)
        % load the spindle file
        load(tmpflist(i).name);
        % get the epoch number
        dash = find(tmpflist(i).name == '-');
        e = str2num(tmpflist(i).name((dash(1)+1):(dash(2)-1)));
        
        % Check if spinpple band exists - may not be true for reference(Shantanu-Jun2013)
        if isfield(spindle{d}{e}{t},'data')
            % convert the spindle envelope field to double
            renv = double(spindle{d}{e}{t}.data(:,3));
            
            % smooth the envelope:
            samprate = spindle{d}{e}{t}.samprate;
            kernel = gaussian(smoothing_width*samprate, ceil(8*smoothing_width*samprate));
            renv = smoothvect(renv, kernel);
            % find the spindles
            % calculate the duration in terms of samples
            mindur = round(min_suprathresh_duration * samprate);
            
            % calculate the threshold in uV units
            if (stdev == 0)
                stdev = std(renv);
            end
            if (baseline == 0)
                baseline = mean(renv);
            end
            thresh = baseline + nstd * stdev;
            
            % extract the events if this is a valid trace
            if (thresh > 0) & any(find(renv<baseline))
                tmprip = extractevents(renv, thresh, baseline, 0, mindur, 0)';
                % Assign the fields
                % start and end indeces
                spin.startind = tmprip(:,1);
                spin.endind = tmprip(:,2);
                % middle of energy index
                spin.midind = tmprip(:,8);
                
                %convert the samples to times for the first three fields
                spin.starttime = spindle{d}{e}{t}.starttime + spin.startind / samprate;
                spin.endtime = spindle{d}{e}{t}.starttime + spin.endind / samprate;
                spin.midtime = spindle{d}{e}{t}.starttime + spin.midind / samprate;
                spin.peak = tmprip(:,3);
                spin.energy = tmprip(:,7);
                spin.maxthresh = (tmprip(:,9) - baseline) / stdev;
                
                % position time element
                postime = pos{d}{e}.data(:,1);
                
                idxs = (1:length(postime))';
                spin.posind = floor(interp1q(postime,idxs,spin.midtime));
            else
                spin.startind = [];
                spin.endind = [];
                spin.midind = [];
                spin.starttime = [];
                spin.endtime = [];
                spin.midtime = [];
                spin.peak = [];
                spin.energy = [];
                spin.posind = [];
            end
            if any(find(renv<baseline))==0
                warning(['No below baseline values in data.  Fields left blank, ',tmpflist(i).name])
            end
            
            spin.timerange = [0 length(renv)/samprate] + spindle{d}{e}{t}.starttime;
            spin.samprate = spindle{d}{e}{t}.samprate;
            spin.threshold = thresh;
            spin.baseline = baseline;
            spin.std = stdev;
            spin.minimum_duration = min_suprathresh_duration;
            
            %spindle
            clear feeg
            spindles{d}{e}{t} = spin;
            clear spin;
            
        end % end isfield data
        
        % If commented, baseline and std from ep1 used for all epochs
        stdev=0; %%% IMP - I want to get standard deviation from each epoch separately
        baseline=0; %%% And baseline too
        
    end
end

save(sprintf('%s/%sspindles%02d.mat', directoryname, fileprefix, d), 'spindles');
%save(sprintf('%s/%sspindles%02d.mat', directoryname, fileprefix, d), 'spindles');
popd;
