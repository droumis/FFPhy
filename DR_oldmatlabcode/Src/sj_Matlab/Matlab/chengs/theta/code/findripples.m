function ripplesnippets = findripples(ripplestruct,minimum_duration,nstd,varargin)
%function ripplesnippets = findripples(ripplestruct,threshold,baseline,minimum_duration,nstd,varargin)
%
%   ripplesnippets = findripples( ... 
%       ripplestruct, ...
%       threshold,baseline,min_suprathreshold_duration, ...
%       [time])
%       nstd:   threshold this many std above baseline 
%
% Scan the ripple eeg struct to snip out start and stop time indices 
% for ripple events. Events whose amplitude envelope exceed threshold
% will be snipped out over points which are above baseline (using a 
% watershed algorithm which scans from the threshold-crossings down to 
% baseline). Overlapping snippets will be combined intelligently, and 
% snippets which do not exceed threshold for more than minimum_duration 
% will be culled.
%
% threshold and baseline are to be expressed in eeg amplitude units
% min_suprathreshold_duration is the time (in seconds) which the signal 
% must remain above threshold to be counted as as ripple; this guards
% against short spikes in signal (e.g. noise) as being counted as 
% ripples. Set min_suprathreshold_duration to some small value, like
% 0.015 s.
%
% posstruct is an optional posstruct. if the argument is provided, 
% the function looks up the index to the nearest position sample which
% precedes the ripple
%
% Returns a ripplesnippets struct with some rather self-explanatory text and 
% scalar fields, and a .data field with a lot of stuff.
%
% Each cell of ripplesnippets.data is an Nx9 matrix, in which each row 
% corresponds to a ripple event. The first column is the start time of 
% each ripple, the second column is the end time, the third column is
% the time of the "middle" (useful for arbitrary alignment of multiple 
% ripple events). The fourth column is the peak height of the envelope 
% for the ripple event. The fifth, sixth, and seventh columns are indices 
% into the ripplestruct.data which correspond to the start, end, and 
% middle times in columns 1-3. The eight column is an index into the 
% nearest preceding position sample in the corresponding posparams struct
% for the epoch, and the ninth column is a position interpolation factor 
% for interpolating between that position sample and the next one.
% These position indices are aligned to the middle of the ripple and allow 
% for quick look-up of the animal's location etc. at time of the ripple.
%
% Written by smk.

if ~numel(varargin)
    warning('You did not provide a posstruct argument. Some data fields will be left blank.');
end

% define the standard deviation for the Gaussian smoother which we
% apply before thresholding (this reduces sensitivity to spurious 
% flucutations in the ripple envelope)
SMOOTHING_WIDTH = 0.002; % 2 ms

% take the Hilbert-transform amplitude envelope of the ripple trace and
% smooth with a Gaussian with standard deviation of SMOOTHING_WIDTH
% which spans 4 SD (note that this uses the gausswin function
% of the Signal Processing toolbox)
gaussian_kernel = gausswin(ceil(8*SMOOTHING_WIDTH*ripplestruct.samprate),4);
% normalize
gaussian_kernel = gaussian_kernel/sum(gaussian_kernel);
% filter the ripple-filtered trace
smoothed_envelope = filtfilt(gaussian_kernel,[1],ripplestruct.data(:,3));

baseline= mean(smoothed_envelope);
dev= std(smoothed_envelope);
threshold= baseline+nstd*dev;

if threshold<baseline; warning('threshold<baseline'); end

fprintf(1, 'baseline= %.2f, threshold= %.2f, std= %.2f\n', baseline, threshold, dev);
%keyboard

% now determine whether the smoothed Hilbert envelope of the ripple trace
% falls above, below, or right on this threshold within each LIA interval
% we keep this information in a cell array called flags, whose organization 
% parallels that of LIA_intervals
flags = sign(smoothed_envelope - threshold);
% within each cell of flags (i.e. for each LIA interval), we determine which
% samples fall above or below threshold according to the flags in flags
% if there exist zero entries (i.e. where the envelope exactly meets the 
% threshold), we apply a greedy neighbor iteration to flag them as +1 (above) 
% or -1 (below)
if ~all(flags)
    while true % keep looping until the break condition is met below
	    borderline_samples = find(flags==0); % find the indices to all 
                                             % elements of ripplestruct.data(:,1) 
        if isempty(borderline_samples); break; end
                                             % that are exactly at threshold
        % discard endpoint entries if they correspond to very first or very 
        % last eeg samples
    	if (borderline_samples(1)==1)
	        borderline_samples(1) = [];
		elseif (borderline_samples(end)==length(ripplestruct.data))
		    borderline_samples(end) = [];	
    	end
        % for each zero flag in flags, check the immediately preceding and 
        % following neighbor flags. if either of these is +1 (i.e. above 
        % threshold), then we "greedily" bump up this element to +1 as well 
        % (it makes sense if you think about it carefully).
        % the + sign here at the front of the parentheses serves to convert 
        % the logical values of the boolean equality tests to numerics note
        % that I am using the single pipe | OR operator, not the double OR
        % this is because I want to compute element-by-element OR to return 
        % a vector of results, not a scalar over the whole vector
        flags_to_flip_positive = +((flags(borderline_samples(:)-1)==1) | ... 
        (flags(borderline_samples(:)+1)==1));
        if any(flags_to_flip_positive)
            flags(borderline_samples(:))=flags_to_flip_positive;
        else
        % we keep doing this until there is no more possibility of change 
        % in the elements of flags
            break
        end
    end
    % and once we exit the while block, all remaining zero flags are 
    % converted to -1
    flags(flags==0) = -1;
end

% now proceed, having gotten rid of all those pesky zero values
% compute first-order differences. ascending threshold-crossings give diff of 
% +2, whereas descending threshold-crossings give diff of -2. note that 
% ascending and descending crossings must alternate; this makes it easy to
% identify each interval. Also need to be careful about the boundary cases 
% (what happens at the beginning or end of the ripple trace)
threshold_crossings = diff(flags);
% the +1 on the next line is to account for the way diff works
ascending_crossings = find(threshold_crossings == +2) + 1; 
% if the trace begins above threshold, then remove first element of 
% descending_crossings
descending_crossings = find(threshold_crossings == -2);
if (descending_crossings(1) < ascending_crossings(1)) 
    % we don't want to grab a truncated ripple at the beginning of the eeg trce
    descending_crossings(1) = [];
% likewise at the end of the trace
elseif (ascending_crossings(end) > descending_crossings(end)) 
    ascending_crossings(end) = [];
end

% we should end up with equal numbers of ascending and descending crossings
if length(ascending_crossings)~=length(descending_crossings)
    error('There is a bug in the calculation of threshold-crossings!')
end

    
% okay, now that we have the threshold-crossings of minimum duration 
% collected, we now have a basic handle on all ripple events. However, there 
% remain two problems: (1) we haven't snipped each ripple event in its 
% entirety, so we may be missing the beginning and end which fall below our 
% threshold; and (2) we may be snipping out adjacent portions of the same 
% ripple event. We need to iron out these issues and then populate the output 
% array.

% construct a look-up vector of indices where the Hilbert envelope falls 
% below the watershed level (which should be set so that it is somewhere in 
% the "noise")
edge_samples = find(smoothed_envelope < baseline); 
% ripples that are too close to the beginning or end of the record 
% may be truncated before they are able to descend into the baseline noise
% we check whether any of the ripple snippets that have been identified 
% are cut off at the endpoints
while true
    if isempty(find(edge_samples > descending_crossings(end)))
        descending_crossings(end) = [];
        ascending_crossings(end) = [];
    elseif isempty(find(edge_samples < ascending_crossings(1)))
        ascending_crossings(1) = [];
        descending_crossings(1) = [];
    else
        break
    end
end

% again, check that we end up with equal numbers of ascending and descending 
% crossings
if length(ascending_crossings)~=length(descending_crossings)
    error('There is a bug in the culling of truncated ripples!')
end
% now we allocate an array to hold the output 
output = NaN*ones(length(ascending_crossings),9);

% now scan before and after each above-treshold ripple top (think of an 
% iceberg) to find the nearest sample that is below noise level
for j = 1:size(output,1)
    output(j,5) = ... 
    max(edge_samples(edge_samples < ascending_crossings(j))) - 1;
    output(j,6) = ... 
    min(edge_samples(edge_samples > descending_crossings(j))) - 1;
end

% test for overlap and combine/delete as necessary
while 1
    can_quit_iterating = true;
    for j = 1:(size(output,1)-1)
        if (output(j,6) >= output(j+1,5))
            % if I've found something to change, I can't quit iterating
            can_quit_iterating = false; 
            output(j,6) = output(j+1,6);
            output(j+1,5) = NaN; 
            % note that comparisons against NaN always give false; this
            % bypasses the test condition on the next go-through
        end
    end
    % now delete the rows correspond to intervals that have already been 
    % accounted for in overlap
    [rows,cols] = find(isnan(output(:,5)));
    output(rows,:) = [];
    if can_quit_iterating
    % this is true only if I've scanned the entire dataset without 
    % changing anything
        break
    end
end

% now remove those ripples whose threshold-crossings are too brief - the 
% ripple envelope should remain over threshold for a certain duration of 
% time, otherwise, we don't treat it as a real ripple.
while 1
    can_quit_iterating = 1;
    for j = 1:size(output,1)
        above_threshold = find( ...
        smoothed_envelope(output(j,5):output(j,6)) > threshold );
        if (length(above_threshold) < minimum_duration*ripplestruct.samprate)
            % if I've found something to change, I can't quit iterating
            can_quit_iterating = 0; 
            output(j,5) = NaN; % put in sentinels to target deletion of rows
        end
    end
    % now delete ripples that do not satisfy minimum_duration
    rows = find(isnan(output(:,5)));
    output(rows,:) = [];
    if can_quit_iterating 
    % this is true only if all ripples scanned are compliant
        break
    end
end

% compute the "median" (in time) of the energy of each ripple event
% also, compute the peak height of the envelope within the ripple event
for j = 1:size(output,1)
    % grab a slice of the ripple trace, square the values, compute cumulative 
    % sums
    ripple_energy_cumulatives = cumsum( ... 
    ripplestruct.data(output(j,5):output(j,6),1).^2 );
    % find the index of the element of ripple_energy_cumulatives which is 
    % closest to 1/2 of the total, and assign this to be the "center" of the 
    % ripple event in the second column of the output array
    [value, index] = min( ...
    abs(ripple_energy_cumulatives - ripple_energy_cumulatives(end)/2) );
    output(j,7) = output(j,5) + index - 1;
    output(j,4) = max(ripplestruct.data(output(j,5):output(j,6),3));
end

% figure out the actual times (in seconds) which correspond to these first and
% last eeg samples for each ripple snippet
timestamps = ripplestruct.starttime + ...
    ( (1:size(ripplestruct.data,1))' - 1 )/ripplestruct.samprate;
output(:,1) = timestamps(output(:,5));
output(:,2) = timestamps(output(:,6));
output(:,3) = timestamps(output(:,7));

% now figure out the index to the position sample which most closely precedes
% the ripple time, and an interpolation factor for how far to the next 
% position sample

% if the posstruct argument is specified, find the indices to the nearest
% position samples which immediately precede the ripples and compute 
% an interpolation factor for how far in between this position sample 
% and the next we need to interpolate to get an exact position estimate
% which corresponds to the time of the ripple
if numel(varargin)
    time=varargin{1};
    postimestamps = time(:,1);
    idxs = (1:length(postimestamps))';
    % I'm doing a clever trick here with interp1q to extract indices 
    % ... think about it if it's not obvious
    output(:,8) = floor(interp1q(postimestamps,idxs,output(:,3)));
    % figure out a factor between 0 and 1 to determine relative 
    % weighting of the posindex and posindex+1 position samples 
    % pos_interp=1 means use the posindex+1 position sample, 
    % whereas pos_interp=0 means use posindex position sample
    valid_idxs = find( ...
        ~isnan(output(:,8)) & ... 
        (output(:,8) > 0) & ...
        (output(:,8) < size(time,1)) );
    output(valid_idxs,9) = ... 
        ( output(valid_idxs,3) - ... 
          time(output(valid_idxs,8),1) ) ...
        ./ ...
        ( time(output(valid_idxs,8)+1,1) - ...
         time(output(valid_idxs,8),1) );
end

% now package output into a struct
ripplesnippets = struct('descript',strvcat( ...
  'array indices for start center, end of ripple events occurring in dataset:', ...
  ripplestruct.descript) );
% define the time window over which we are counting the ripples
ripplesnippets.timerange = [ timestamps(1) timestamps(end) ];
ripplesnippets.samprate = ripplestruct.samprate;
ripplesnippets.threshold = threshold;
ripplesnippets.std = dev; 
ripplesnippets.baseline = baseline; 
ripplesnippets.minimum_duration = minimum_duration;
ripplesnippets.fields = 'starttime1 endtime2 midtime3 peakheight4 startidx5 endidx6 mididx7 posidx8 posinterp9';
ripplesnippets.data = output;

