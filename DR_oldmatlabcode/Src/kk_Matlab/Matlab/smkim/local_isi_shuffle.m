function shuffled = local_isi_shuffle(unit,window_range,num_shuffles)
%LOCAL_ISI_SHUFFLE Shuffle spike trains, preserving local rate and ISI distribution
%
%   SHUFFLED = LOCAL_ISI_SHUFFLE(UNIT,WINDOW_SIZE,NUM_SHUFFLES) shuffles the
%   spike train in UNIT. UNIT must be a scalar struct containing single-unit
%   spike data, of the type that validates with IS_UNIT (written by SMK).
%   WINDOW_SIZE must be a 2-element row vector of real finite positive
%   floating-point scalars, such that 
%
%       2*WINDOW_SIZE(1) >= WINDOW_SIZE(2) > WINDOW_SIZE(1)
%
%   WINDOW_SIZE specifies, in seconds, the timescale on which the spike train is
%   to be shuffled. The timerange of UNIT is divided into contiguous windows
%   whose durations are drawn from a uniform distribution on WINDOW_SIZE, and
%   spikes are shuffled locally over each window in a manner which preserves
%   interspike intervals. WINDOW_SIZE should be chosen to preserve overall
%   slow modulations of firing rate. Shuffling will be poor if the spiking rate
%   is not greater than ~2 spikes/window.
%
%   NUM_SHUFFLES must be a positive real integer.
%
%   The output SHUFFLED is a NUM_SHUFFLES-by-1 column vector struct containing
%   metadata fields inherited from UNIT, but with shuffled spike times in the
%   'timestamp' field of each element.
%
%   Reference:
%   [1] Rivlin-Etzion M., Ritov Y., Heimer G., Bergman H., Bar-Gad I. (2006)
%       Local shuffling of spike trains boosts the accuracy of spike train
%       spectral analysis. _Journal of Neurophyiology_ 95:3245-3256.
%
%Depends on:
%   RANDSAMPLE (MATLAB Statistics Toolbox)
%   IS_UNIT (written by SMK)
%
%Written by SMK, 2010 January 29.
%

TS_PER_SEC = 1e4;

if (exist('is_unit') ~= 2)
  error(['LOCAL_ISI_SHUFFLE depends on m-file IS_UNIT ' ...
      '(written by SMK)']);
end
if (exist('randsample') ~= 2)
  error(['LOCAL_ISI_SHUFFLE depends on m-file RANDSAMPLE ' ...
      '(MATLAB Statistics Toolbox)']);
end

if ~is_unit(unit) || ~isscalar(unit)
  error('UNIT must be a scalar struct containing single-unit spike data');
end

if ~isreal(window_range) || ~isfloat(window_range) || ...
    any(~isfinite(window_range(:))) || ~isequal(size(window_range),[1 2]) || ...
    ~all(window_range > 0) || (window_range(2) <= window_range(1)) || ...
    (window_range(2) > 2*window_range(1))
  error('WINDOW_SIZE is not valid. See the help for LOCAL_ISI_SHUFFLE.');
end

if ~isscalar(num_shuffles) || ~isreal(num_shuffles) || ...
    ~isfinite(num_shuffles) || (num_shuffles <= 0) || ...
    (round(num_shuffles) ~= num_shuffles)
  error('NUM_SHUFFLES must be a positive real integer');
end

% allocate output
REQUIRED_FIELDS = { ...
    'uid'                     , ...
    'subject'                 , ...
    'day'                     , ...
    'epoch'                   , ...
    'environment'             , ...
    'timerange'               , ...
    'tetrode'                 , ...
    'depth'                   , ...
    'hemisphere'              , ...
    'region'                  , ...
    'reference'               , ...
    'passbands'               , ...
    'thresholds'              , ...
    'Fs'                      , ...
    'samples_before_trigger'  , ...
    'amplitude_cutoff'        , ...
    'timestamp'               , ...
    'samples'                 , ...
    'clustnum'                , ...
    'shape'                   , ...
    'cluster_quality'         , ...
    'sources'                 };
fn = fieldnames(unit);
shuffled = repmat(rmfield(unit,fn(~ismember(fn,REQUIRED_FIELDS))), ...
    [num_shuffles 1]);

% Get original spike times and ISIs
t = double(unit.timestamp)/TS_PER_SEC;
isi = diff(t);
for i = 1:num_shuffles
  % Bins of random duration drawn from window_range; to ensure adequate coverage,
  % we generate more than we need
  bin_edges = double(unit.timerange(1))/TS_PER_SEC + ...
      cumsum(window_range(1) + diff(window_range) .* ...
      rand([ceil(double(diff(unit.timerange))/TS_PER_SEC/min(window_range)) 1]));
  % Trim excess so that the randomized bins exactly coincide with the total
  % duration of observations
  bin_edges = bin_edges(1 : ...
      find(bin_edges <= double(unit.timerange(end))/TS_PER_SEC,1,'last'));
  bin_edges(end) = double(unit.timerange(end))/TS_PER_SEC;
  n_bins = numel(bin_edges)-1;
  % Now for each bin, find the spikes in the local 3-bin neighborhood that are
  % nearest in time to the boundaries of the bin, and use these as stationary
  % "pivots" for shuffling intervening ISIs. If <2 spikes occur in the local
  % 3-bin neighborhood, then we can't shuffle without changing the ISIs so we
  % just do nothing.
  t_shuffled = t;
  for j = 1:n_bins
    % Find spikes that are closest to bin_edges([i, i+1])
    [junk, idx_start] = min(abs(t - bin_edges(j)));
    [junk, idx_end] = min(abs(t - bin_edges(j+1)));
    % Don't shuffle if an ordered pair of nearest spikes are not available in
    % the target bin or an adjacent bin, or if there are no intervening spikes
    % between this pair
    if isempty(idx_start) || isempty(idx_end) || (idx_end <= idx_start+1) || ...
        (t(idx_start) < bin_edges(j - double(j > 2))) || ...
        (t(idx_start) >= bin_edges(j+1)) || ...
        (t(idx_end) >= bin_edges(j + 1 + double(j < n_bins))) || ...
        (t(idx_end) < bin_edges(j))
      continue;
    end
    % Shuffle intervening spike times, preserving ISIs
    t_shuffled((idx_start+1):(idx_end-1)) = t(idx_start) + ...
        cumsum(randsample(isi(idx_start:(idx_end-1)),idx_end-idx_start-1));
  end
  shuffled(i).timestamp = uint32(t_shuffled * TS_PER_SEC);
end

% Add information about shuffling
[shuffled(:).shuffle_window_range] = deal(window_range);


