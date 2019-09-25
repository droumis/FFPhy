function spikes = ss_dejitter(spikes, method, maxshift)

% SS_DEJITTER  Library of methods aligning waveform peaks.
%     SPIKES = SS_DEJITTER(SPIKES) takes and returns a spike-sorting object
%     SPIKES and attempts to dejitter the spike waveforms using the 'com'
%     method with a maximum shift of 3 (see below).
%
%     'Dejittering' refers to the registration of digitally sampled/thresholded 
%     waveforms.  The sampling/thresholding procedure introduces noise due to
%     (1) variation in relative timing between an analog waveform and the
%     sample clock and (2) variation in the time of threshold crossing due to
%     noise.  The resulting misalignment even in otherwise identical analog 
%     waveforms is known as jitter.  Dejittering involves estimating the location
%     of a fiducial (e.g., the central peak) for each waveform and aligning these
%     markers across waveforms.
%
%     SPIKES = SS_DEJITTER(SPIKES, METHOD, MAXSHIFT) also takes a METHOD argument
%     (default: 'com') that selects the dejittering method and MAXSHIFT (default: 3),
%     specifying the maximum shift to be applied during dejittering.
% 
%     The number of samples per waveform will be altered by this function when 
%     samples are made invalid by realignment (i.e., the amount of shift requires
%     extrapolating past the boundaries of the original waveform) and the returned
%     SPIKES object will have its waveforms clipped to exclude invalid regions.
%     The MAXSHIFT argument limits this by specifying the maximum allowed shift
%     (in samples).  Waveforms requiring more than a MAXSHIFT samples realignment
%     are taken to be outliers and will not be shifted.
% 
%     All methods require the definition of threshold crossing time (threshT) and
%     threshold level vector (threshV) with both low and high thresholds (use +/-
%     Inf for either the low or the high threshold if it was not used in extracting
%     the waveforms).  Valid choices for 'method' are:
%
%         'peak' upsamples (via cubic spline interpolation) the region of each waveform
%         after threshT that remains below/above threshold and chooses the fiducial
%         as min/max value (depending on polarity) of the upsampled waveform in this
%         region.  The resolution is limited by the upsample factor (currently 10).
% 
%         'com' also considers the region following threshT that remains below/above
%         threshold.  It finds the center of mass of the resulting points using
%         (sum_t[t * (v_t - thresh)] / sum_t[(v_t - thresh)]).
%         
%     The interpolation to align the waveforms is always performed using a cubic
%     spline, regardless of the method use to find the fiducial.
%
% References:
%     Fee MS et al (1996).  J. Neurosci Methods (69): 175-88
%     Sahani M (1999).  PhD thesis, Pasadena, CA: Caltech.
%
% Last Modified: sbm, 8/18/03

starttime = clock;

%%%%%%%%%% DEFAULTS
if (nargin < 2)
    method = 'com';  % default
end
if (nargin < 3)
    maxshift = 3;    % default max shift
end


%%%%%%%%%% SPIKES FIELDS WHICH ARE BEING CHANGED: KEEP AN ORIGINAL COPY

spikes.oriwaveforms = spikes.waveforms;
spikes.orithreshT = spikes.threshT;
spikes.orithreshV = spikes.threshV;


%%%%%%%%%% ARGUMENT CHECKING

if (~isfield(spikes, 'waveforms') | (size(spikes.waveforms, 1) < 1))
    error('SS:waveforms_undefined', 'The SS object does not contain any waveforms!');
elseif (~isfield(spikes, 'threshT'))
	error('SS:threshT_undefined', 'The SS object must define the sample index of threshold crossing.');
elseif (~isfield(spikes, 'threshV') | (length(spikes.threshV) ~= 2))
	error('SS:threshV_error', 'The SS object must define both low and high thresholds.  Use +/- Inf if only one of the two was used.'); 
elseif ((spikes.threshV(2) - spikes.threshV(1)) < 0)
    error('SS:threshV_illegal', 'The SS object high threshold must be greater than its low threshold.');
end

%%%%%%%%%% CONSTANTS
spline_chunk = 5000;   % # of splines to do at a time; keeps memory use down
upsample = 10;    % amount of upsampling to find extrema using 'peak' method
numspikes = size(spikes.waveforms, 1);
numsamples = size(spikes.waveforms, 2);

maxchunk = floor(((numspikes - 1) / spline_chunk));

%%%%%%%%%% FIND FIDUCIALS
switch (lower(method)),
    case 'peak',   % Fiducials are location of max dist from threshold during central peak
        fiducials = zeros(numspikes, 1);
        % We do this by upsampling waveforms and using 'thresholded_peaks', but to
        % avoid the memory hit of upsampling all waveforms, we do this in pieces
        % ... the chunks are big enough that the for loop overhead should be minimal.
        for group = 0:maxchunk
            startind = (group * spline_chunk) + 1;
            if (group < maxchunk)
                inds = [startind:((group + 1) * spline_chunk)];  % next chunk
            else
                inds = [startind:numspikes];                     % leftover spikes
            end
            upspikes.waveforms = spline(1:numsamples, spikes.waveforms(inds,:), ...
                                        1:1/upsample:numsamples); % upsample
            upspikes.threshT = (spikes.threshT * upsample) - upsample + 1; % rescale time
            upspikes.threshV = spikes.threshV;
            [h, w, peakLocs] = thresholded_peaks(upspikes);  % get peak locs ...
            fiducials(inds) = (peakLocs(:,2) + (upsample - 1)) ./ upsample; % ... & unscale time
        end
        clear upspikes;  % clean up
        
    case 'com',
        % first, get a mask over the peaks
        isThreshLo = (spikes.waveforms(:, spikes.threshT) < spikes.threshV(1));
        [h, w, pl, mask] = thresholded_peaks(spikes);

        % next, subtract appropriate thresholds, making all values in the peak
        % have positive sign, since we're about to treat them as 'mass'
        waves = spikes.waveforms;
        waves(~isThreshLo, :) = waves(~isThreshLo, :) - spikes.threshV(2);
        waves( isThreshLo, :) = spikes.threshV(1) - waves( isThreshLo, :);
        waves = waves .* mask;
  
        % now COM is straightforward
        fiducials = (waves * [1:numsamples]') ./ sum(waves, 2);

        clear mask waves;  % manual garbage collection
    otherwise,
        error('SS:method_invalid', 'Unknown dejittering method.');
end

%%%%%%%%%% REALIGN FIDUCIALS
% We line up the fiducials around the sample that requires the least shifting.
target = round(mean(fiducials));

% Determine shifted indices for each waveform
shifts = fiducials - target;
shifts(abs(shifts) > maxshift) = 0;   % big shifts are outliers; don't alter these
resample_inds = repmat([1:numsamples], [numspikes, 1]) + repmat(shifts, [1, numsamples]);

% Which regions are invalid?
left_valid = max(find(any(resample_inds < 1))) + 1;
right_valid = min(find(any(resample_inds > numsamples))) - 1;
if (left_valid >= right_valid)
    error('SS:waveform_invalidated', 'Realignment invalidates all samples.  Try reducing MAXSHIFT.');
end

% Although 'spline' can find independent interpolants for each row of a matrix
% (warning: the help for 'spline' is misleading on rows vs. cols), it does not
% allow each row to be resampled on its own grid.  To get around this (since the
% for-loop approach is ~4x slower), we obtain the piecewise polynomial coefficients
% along each waveform, pad to handle endpoints, and then string them all together
% into one long series.  This will let us compute the shifted waveforms in one 
% fell swoop . . . the only caveat is that requests to extrapolate beyond
% [1:numsamples] become meaningless; but we're not interested in these anyway.
pp = spline([1:numsamples], spikes.waveforms);

% The coefficients from 'spline' come out ordered by column rather than by
% row/waveform, so we reorder.
pp.coefs = reshape(pp.coefs, numspikes, numsamples-1, []);
pp.coefs = permute(pp.coefs, [2 1 3]);

% 'ppeval' uses the right spline at endpoints; i.e., if we just strung together
% the splines as is, the last sample of each waveform would be interpolated
% using the first polynomial from the next waveform.  So we add a DC polynomial
% to the end of each spike whose value equals to the last sample of the spike.
% Now after stringing together, interpolation at ((r-1)*numsamples + [1:numsamples])
% gives the rth spike (within roundoff error).
padzeros = zeros(1, numspikes, 4);
pp.coefs = cat(1, pp.coefs, padzeros);
pp.coefs(numsamples,:,4) = spikes.waveforms(:,end)';

% Make the 'strung-together' version.
pp.coefs = reshape(pp.coefs, [], 4);

% Rework the remaining information in the piecewise polynomial to be consistent.
pp.pieces = numsamples * numspikes;
pp.dim = 1;
pp.breaks = [1:(pp.pieces+1)];

% Change resample indices to correspond to the strung-together piecewise polynomial.
offset = ([1:numspikes] - 1)' * numsamples;
resample_inds = resample_inds + repmat(offset, [1, numsamples]);

% OK, finally.  Do the resample & clip the invalid regions
resampled = ppval(pp, resample_inds);
spikes.waveforms = resampled(:, left_valid:right_valid);

% Also, the threshold crossing time has changed due to data clipping
spikes.threshT = spikes.threshT - left_valid + 1;

% Worse, the threshold time/value may have changed while shifting.  This
% will make it hard to identify peaks later, so we need new values.
meanlevel = mean(spikes.waveforms(:));   % estimate of baseline
notallowed = [1:(spikes.threshT-maxshift-1),(spikes.threshT+maxshift+1):(size(spikes.waveforms, 2))];

deflection = mean(abs(spikes.waveforms - meanlevel), 1);  % mean deviation from baseline
deflection(notallowed) = 0;
[junk, spikes.threshT] = max(deflection);  % max deflect near old threshT marks new threshold marker
threshValues = spikes.waveforms(:, spikes.threshT);  % values at new marker
spikes.threshV(1) = max([threshValues(threshValues < meanlevel); -Inf]); % value closest to baseline but below it
spikes.threshV(2) = min([threshValues(threshValues > meanlevel); Inf]); % value closest to baseline but above it 

spikes.tictoc.dejitter = etime(clock, starttime);
