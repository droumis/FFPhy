function spectra = spike_spectrogram(unit,params)
%SPIKE_SPECTROGRAM Compute short time-window spectra of single-unit spike train.
%
%   SPECTRA = SPIKE_SPECTROGRAM(UNIT,PARAMS) takes a struct UNIT containing
%   single-unit spike data (of the type which validates with IS_UNIT) and
%   computes the spectrum of the spike train in short time windows using the
%   parameters specified in PARAMS. 
%
%   PARAMS must be a scalar struct with the following fields:
%     'discretization_timestep': a nonzero uint32 scalar, specifying the
%         resolution with which to discretize the spike train (in timestamp
%         units of 1e-4 seconds)
%     'window_size': a real finite floating-point scalar, specifying the
%         duration of the short time windows in seconds. window_size/TW is the
%         halfwidth of frequency resolution of the spectral estimate. An error
%         will be raised if discretization_timestep is not a divisor of 
%         window_size*1e4.
%     'TW': positive real scalar >= 1, specifying the time-bandwidth product.
%         Larger values give lower variance at the expense of spectral
%         resolution. floor(2*TW - 1) Slepian tapers are applied to the data
%         and the resulting spectral estimates are averaged.
%     'nfft': a finite positive real integer, specifying the number of points
%         in the FFT. The computation is fastest when nfft is a power of two,
%         and slowest when nfft has large prime factors.
%     'frequency_band': a 2-element row vector containing positive real
%         floating-point values, specifying the range of frequencies at which
%         the spectral estimate is to be computed. frequency_band(2) must be
%         greater than frequency_band(1) and less than or equal to the Nyquist
%         frequency (which is set by discretization_timestep). An error will be
%         raised if the band is narrower than the spacing of the frequency grid
%         that is specified by nfft.
%
%     Additionally, PARAMS must have one or the other of the following fields
%     (not both):
%       'window_overlap': a real scalar on the interval [0,1), specifying the
%           fractional overlap between successive short time windows.
%       'window_centers': a cell array of the same size as UNIT, in which each
%           cell contains a vector of uint32 timestamps that specify the
%           center-times of the short time windows. An error will be raised if
%           these times are not compatible with window_size and the available
%           timerange of the spike data. If UNIT is a scalar struct, then the
%           packaging cell can be omitted and window_centers can simply be a
%           uint32 vector.
%     If 'window_overlap' is specified, then the window centers are spaced
%     evenly to tile the timerange of the data as completely as possible (with
%     loss of fractional-window data at the start and end).
%
%   SPECTRA is a struct of the same size as UNIT with the following fields:
%       'uid': inherited from UNIT
%       'subject': inherited from UNIT
%       'day': inherited from UNIT
%       'epoch': inherited from UNIT
%       'environment': inherited from UNIT
%       'tetrode': inherited from UNIT
%       'depth': inherited from UNIT
%       'hemisphere': inherited from UNIT
%       'region': inherited from UNIT
%       'reference': inherited from UNIT
%       'passbands': inherited from UNIT
%       'thresholds': inherited from UNIT
%       'clustnum': inherited from UNIT
%       'cluster_quality': inherited from UNIT
%       'discretization_timestep': inherited from PARAMS
%       'window_size': inherited from PARAMS
%       'window_overlap': inherited from PARAMS if defined, otherwise NaN
%       'nfft': inherited from PARAMS
%       'TW': inherited from PARAMS
%       'num_tapers': equal to (2*TW - 1)
%       'timestamp': vector of uint32 timestamps, specifying the window
%           centers. equal to PARAMS.window_centers if it was defined.
%       'frequency': vector of real nonnegative floating-point values,
%           specifying the frequency grid of the spectral estimates
%       'power_spectral_density': a matrix of size 
%           [numel(timestamp), numel(frequency)], containing estimated power
%           spectral density for each time window, expressed in spikes/second.
%       'jackknife_power_spectral_density': a 3-dimensional array of size 
%           [numel(timestamp), numel(frequency), number_of_jackknife_samples]
%           (empty if only one taper was used)
%       'spike_count': vector whose size is [numel(timestamp) 1], counting the
%           number of spikes that contribute to the spectral esimate for each
%           short time window
%
%   References:
%   [1] Percival D.B., Walden A.T. (1993) _Spectral analysis for physical
%       applications: multitaper and conventional univariate techniques_.
%       Cambridge University Press.
%   [2] Mitra P.P., Pesaran B. (1999) Analysis of dynamic brain imaging data.
%       _Biophysical Journal_ 76:691-708.
% 
%Depends on:
%   IS_UNIT (written by SMK)
%   IS_INTERVALS (written by SMK)
%   ISMEMBER_INTERVALS (written by SMK)
%   DPSS (MATLAB Signal Processing toolbox)
%
%Written by SMK, 2010 January 30.
%

TS_PER_SEC = 1e4;

if (exist('is_unit') ~= 2)
  error(['SPIKE_SPECTROGRAM depends on m-file IS_UNIT ' ...
      '(written by smk)']);
end
if (exist('is_intervals') ~= 2)
  error(['SPIKE_SPECTROGRAM depends on m-file IS_INTERVALS ' ...
      '(written by smk)']);
end
if (exist('ismember_intervals') ~= 2)
  error(['SPIKE_SPECTROGRAM depends on m-file ISMEMBER_INTERVALS ' ...
      '(written by smk)']);
end
if (exist('dpss') ~= 2)
  error(['SPIKE_SPECTROGRAM depends on m-file DPSS ' ...
      '(MATLAB Signal Processing toolbox)']);
end

if ~is_unit(unit)
  error('UNIT is not a valid struct containing single-unit spike data');
end
if ~isstruct(params) || ~isscalar(params) || ~all(isfield(params, ...
    {'discretization_timestep', 'window_size', 'TW', 'nfft', ...
    'frequency_band'}))
  error('PARAMS is not a scalar struct with the required fields');
end
if ~xor(isfield(params,'window_overlap'),isfield(params,'window_centers'))
  error(['PARAMS must have either a ''window_overlap'' field or a ' ...
      '''window_centers'' field, but not both']);
end
if ~isscalar(params.discretization_timestep) || ...
    ~isreal(params.discretization_timestep) || ...
    ~isa(params.discretization_timestep,'uint32') || ...
    ~(params.discretization_timestep > 0)
  error(['PARAMS.discretization_timestep must be a real positive uint32 ' ...
      'scalar']);
end
% compute effective sampling rate
Fs = TS_PER_SEC/double(params.discretization_timestep);

if ~isscalar(params.window_size) || ~isreal(params.window_size) || ...
    ~isfloat(params.window_size) || ~isfinite(params.window_size) || ...
    ~(params.window_size > 0)
  error(['PARAMS.window_size must be a real finite positive floating-point ' ...
      'scalar']);
end
% Raise an error if window_size is longer than the data timerange in any
% element of unit
if any(arrayfun(@(u) diff(double(u.timerange)/TS_PER_SEC) < ...
    params.window_size,unit))
  error(['PARAMS.window_size must be shorter than the timerange of the ' ...
      'spike data']);
end
if ~isscalar(params.TW) || ~isnumeric(params.TW) || ~isreal(params.TW) || ...
    ~isfinite(params.TW) || ~(params.TW >= 1)
  error('PARAMS.TW must be a real finite positive scalar');
end

if ~isscalar(params.nfft) || ~isreal(params.nfft) || ...
    ~isfinite(params.nfft) || (params.nfft < 2) || ...
    ~(round(params.nfft) == params.nfft)
  error('PARAMS.nfft must be a real finite integer-valued scalar >1');
end
if ~isequal(size(params.frequency_band),[1 2]) || ...
    ~is_intervals(params.frequency_band) || ...
    ~all(params.frequency_band >= 0) || ...
    ~all(params.frequency_band < ...
    TS_PER_SEC/double(params.discretization_timestep)/2);
  error(['PARAMS.frequency_band must be a 2-element row vector which ' ...
      'specifies a frequency interval between 0 and Nyquist']);
end

if isfield(params,'window_overlap')
  % Check params.window_overlap
  if ~isscalar(params.window_overlap) || ~isreal(params.window_overlap) || ...
      ~isfloat(params.window_overlap) || ~isfinite(params.window_overlap) || ...
      ~ismember_intervals(params.window_overlap,[0 1])
    error(['PARAMS.window_overlap must be a real finite positive ' ...
        'floating-point scalar greater than or equal to zero and less ' ...
        'than one']);
  end
  % Spacing between window centers
  step_size = params.window_size * (1 - params.window_overlap);
  % Generate vector of window centers for each element of UNIT
  params.window_centers = cell(size(unit));
  for i = 1:numel(unit)
    num_windows = 1 + floor( ...
        (diff(double(unit(i).timerange)/TS_PER_SEC) - params.window_size) / ...
        step_size);
    % Divide the data timerange into evenly-spaced overlapping windows,
    % discarding fractional-window leftovers at the start and end
    params.window_centers{i} = unit(i).timerange(1) + uint32(TS_PER_SEC * ...
        ( rem(diff(double(unit(i).timerange)/TS_PER_SEC),step_size)/2 + ...
        params.window_size/2 + step_size*(0:(num_windows-1))' ));
  end
else
  assert(isfield(params,'window_centers'));
  % Repackage params.window_centers as a cell scalar if necessary
  if isscalar(unit) && ~iscell(params.window_centers)
    params.window_centers = {params.window_centers};
  end
  if ~iscell(params.window_centers) || ...
      ~isequal(size(unit),size(params.window_centers)) || ...
      ~all(cellfun(@(c) isvector(c) && isa(c,'uint32') && ...
      (size(c,1) >= 1) && (size(c,2) == 1),params.window_centers))
    error(['PARAMS.window_centers must be a cell array whose size matches ' ...
        'that of UNIT, in which each cell is a column vector of uint32 ' ...
        'timestamps']);
  end
  % Check each cell of params.window_centers against each element of UNIT
  for i = 1:numel(unit)
    if (unit(i).timerange(1) > min(params.window_centers{i}) - ...
        uint32(TS_PER_SEC*params.window_size/2)) || ...
        (unit(i).timerange(end) < max(params.window_centers{i}) + ...
        uint32(TS_PER_SEC*params.window_size/2))
      error(['PARAMS.window_centers and PARAMS.window_size are together ' ...
          'not compatible with the available data timerange']);
    end
  end
  % Assign a sentinel value to params.window_overlap
  params.window_overlap = NaN;
end

% Bin edges for discretizing spike times relative to window center, expressed
% in uint32 timestamp units but cast as double for compatibility with HISTC
Fs = TS_PER_SEC/double(params.discretization_timestep);
bins_per_window = round(Fs*params.window_size);
bin_edges = linspace(-TS_PER_SEC*params.window_size/2, ...
    +TS_PER_SEC*params.window_size/2,bins_per_window+1);

% For a given TW, there exist 2*TW-1 Slepian tapers whose spectral
% concentration is close to optimal
num_tapers = round(2*params.TW - 1);
% Construct Slepian tapers at resolution that matches the discretization of
% spike times
[tapers, eigenweights] = dpss(bins_per_window,params.TW,num_tapers);
% normalization and reorient for later computational convenience
tapers = sqrt(Fs) .* permute(tapers,[3 1 2]);
if (num_tapers > 1)
  assert(isequal(size(tapers),[1, bins_per_window, num_tapers]));
else
  assert(isequal(size(tapers),[1, bins_per_window]));
end

% FFT frequency grid
frequency = linspace(0,Fs,1+params.nfft);
frequency(end) = [];
% subset of the FFT frequencies that span the desired frequency band
f_idx = find(ismember_intervals(frequency,params.frequency_band));
if (f_idx(1) > 1)
  f_idx = [f_idx(1)-1, f_idx];
end
if (f_idx(end) < numel(frequency))
  f_idx = [f_idx, f_idx(end)+1];
end

% allocate output struct
FIELDS_TO_INHERIT = { ...
    'uid'             , ...
    'subject'         , ...
    'day'             , ...
    'epoch'           , ...
    'environment'     , ...
    'tetrode'         , ...
    'depth'           , ...
    'hemisphere'      , ...
    'region'          , ...
    'reference'       , ...
    'passbands'       , ...
    'thresholds'      , ...
    'clustnum'        , ...
    'cluster_quality' };
fn = fieldnames(unit);
spectra = rmfield(unit,fn(~ismember(fn,FIELDS_TO_INHERIT)));
[spectra(:).discretization_timestep] = deal(params.discretization_timestep);
[spectra(:).window_size] = deal(params.window_size);
[spectra(:).window_overlap] = deal(params.window_overlap);
[spectra(:).nfft] = deal(params.nfft);
[spectra(:).TW] = deal(params.TW);
[spectra(:).timestamp] = deal(params.window_centers{:});
[spectra(:).frequency] = deal(frequency(f_idx));

for i = 1:numel(spectra)
  % Count the spikes that occur within each bin of each short time window.
  % spike_counts is transposed so that each row corresponds to a short time
  % window and each column corresponds to a bin within the respective window.
  spike_count = histc(bsxfun(@minus, ...
      double(unit(i).timestamp),double(params.window_centers{i}')), ...
      bin_edges)';
  % HISTC counts out-of-bounds elements in the last position, which we don't
  % need.
  spike_count(:,end) = [];
  assert(isequal(size(spike_count), ...
      [numel(params.window_centers{i}), bins_per_window]));
  % Compute FFT of the projections of the data onto the Slepian tapers and
  % subtract DC component
  J = fft(bsxfun(@times,spike_count,tapers),params.nfft,2) - ...
      bsxfun(@times,mean(spike_count,2),fft(tapers,params.nfft,2));
  % Select frequencies of interest
  J = J(:,f_idx,:);
  % Power is the squared modulus of the fourier transform
  S = J .* conj(J);
  if (num_tapers > 1)
    assert(isequal(size(S), ...
        [numel(params.window_centers{i}), numel(f_idx), num_tapers]));
  else
    assert(isequal(size(S),[numel(params.window_centers{i}), numel(f_idx)]));
  end
  % Average the per-taper estimates to get the multitaper estimate
  spectra(i).power_spectral_density = mean(S,3);
  % Leave one out to compute jackknife estimates
  if (num_tapers > 1)
    spectra(i).jackknife_power_spectral_density = nan(size(S));
    for j = 1:num_tapers
      spectra(i).jackknife_power_spectral_density(:,:,j) = ...
          mean(S(:,:,(1:num_tapers ~= j)),3); 
    end
  else
    % empty array
    spectra(i).jackknife_power_spectral_density = ...
        nan([size(spectra(i).power_spectral_density), 0]);
  end
  spectra(i).spike_count = sum(spike_count,2);
end

