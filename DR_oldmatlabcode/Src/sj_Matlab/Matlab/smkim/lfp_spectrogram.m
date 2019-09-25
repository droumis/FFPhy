function spectra = lfp_spectrogram(lfp,params,varargin)
%LFP_SPECTROGRAM Compute short time-window spectra of continuous LFP data
%
%   SPECTRA = LFP_SPECTROGRAM(LFP,PARAMS) takes a struct LFP containing
%   continuous LFP data (of the type which validates with IS_CONTINUOUS and
%   computes the spectrum of the LFP in short time windows using the parameters
%   specified in PARAMS. 
%
%   PARAMS must be a scalar struct with the following fields:
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
%       'window_centers': a cell array of the same size as LFP, in which each
%           cell contains a vector of uint32 timestamps that specify the
%           center-times of the short time windows. An error will be raised if
%           these times are not compatible with window_size and the available
%           timerange of the LFP data. If LFP is a scalar struct, then the
%           packaging cell can be omitted and window_centers can simply be a
%           uint32 vector.
%     If 'window_overlap' is specified, then the window centers are spaced
%     evenly to tile the timerange of the data as completely as possible (with
%     loss of fractional-window data at the start and end).
%
%   SPECTRA is a struct of the same size as LFP with the following fields:
%       'subject': inherited from LFP
%       'day': inherited from LFP
%       'epoch': inherited from LFP
%       'environment': inherited from LFP
%       'electrode': inherited from LFP
%       'channel': inherited from LFP
%       'depth': inherited from LFP
%       'hemisphere': inherited from LFP
%       'region': inherited from LFP
%       'reference': inherited from LFP
%       'passband': inherited from LFP
%       'Fs': inherited from LFP
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
%       (optionally: 'filter_object' if one was passed as an input)
%
%   References:
%   [1] Percival D.B., Walden A.T. (1993) _Spectral analysis for physical
%       applications: multitaper and conventional univariate techniques_.
%       Cambridge University Press.
%   [2] Mitra P.P., Pesaran B. (1999) Analysis of dynamic brain imaging data.
%       _Biophysical Journal_ 76:691-708.
%
%   SPECTRA = LFP_SPECTROGRAM(LFP,PARAMS,HD) computes the spectra after
%   filtering the LFP signal using a MATLAB dfilt object Hd. This is useful for
%   prewhitening the signal to reduce its dynamic range (and thereby reduce
%   spectral leakage in the spectral estimate). The filter HD is applied
%   *twice*, once forwards and once backwards. It is the user's responsibility
%   to correctly design the filter. When LFP_SPECTROGRAM is called in this
%   manner, SPECTRA will have a field named 'filter_object' which contains a
%   copy of Hd.
%
%Depends on:
%   IS_CONTINUOUS (written by SMK)
%   FILTFILTHD (written by Malcolm Lidierth, MATLAB Central File ID #17061; see 
%       http://www.mathworks.com/matlabcentral/fileexchange/17061-filtfilthd)
%   IS_INTERVALS (written by SMK)
%   ISMEMBER_INTERVALS (written by SMK)
%   DPSS (MATLAB Signal Processing toolbox)
%
%Written by SMK, 2010 January 30.
%

TS_PER_SEC = 1e4;

if (exist('is_continuous') ~= 2)
  error('LFP_SPECTROGRAM depends on m-file IS_CONTINUOUS (written by smk)');
end
if (exist('filtfilthd') ~= 2)
  error(['LFP_SPECTROGRAM depends on m-file FILTFILTHD ' ...
      '(written by Malcolm Lidierth, MATLAB Central File ID #17061)']);
end
if (exist('is_intervals') ~= 2)
  error(['LFP_SPECTROGRAM depends on m-file IS_INTERVALS ' ...
      '(written by smk)']);
end
if (exist('ismember_intervals') ~= 2)
  error(['LFP_SPECTROGRAM depends on m-file ISMEMBER_INTERVALS ' ...
      '(written by smk)']);
end
if (exist('dpss') ~= 2)
  error(['LFP_SPECTROGRAM depends on m-file DPSS ' ...
      '(MATLAB Signal Processing toolbox)']);
end

if ~is_continuous(lfp)
  error(['LFP does not appear to be a valid LFP data struct']);
end
Fs = lfp(1).Fs;
if ~all(arrayfun(@(l) isequal(Fs,l.Fs),lfp))
  error('All element of LFP must share the same ''Fs'' field');
end

if isempty(varargin)
  prefilter = false;
else
  if (length(varargin) > 1)
    error('Too many extra arguments');
  end
  Hd = varargin{1};
  % TODO: Use meta-classes to interrogate class identity
  if ~strncmp('dfilt.',class(Hd),6)
    error('HD must be a MATLAB Signal Processing Toolbox dfilt object');
  end
  prefilter = true;
end

if ~isstruct(params) || ~isscalar(params) || ~all(isfield(params, ...
    {'window_size', 'TW', 'nfft', 'frequency_band'}))
  error('PARAMS is not a scalar struct with the required fields');
end
if ~xor(isfield(params,'window_overlap'),isfield(params,'window_centers'))
  error(['PARAMS must have either a ''window_overlap'' field or a ' ...
      '''window_centers'' field, but not both']);
end

if ~isscalar(params.window_size) || ~isreal(params.window_size) || ...
    ~isfloat(params.window_size) || ~isfinite(params.window_size) || ...
    ~(params.window_size > 0)
  error(['PARAMS.window_size must be a real finite positive floating-point ' ...
      'scalar']);
end
% Raise an error if window_size is longer than the data timerange in any
% element of LFP
if any(arrayfun(@(l) double(l.timestamp(end) - l.timestamp(1))/TS_PER_SEC < ...
    params.window_size,lfp))
  error(['PARAMS.window_size must be shorter than the timerange of the ' ...
      'LFP data']);
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
    ~all(params.frequency_band < Fs/2)
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
  % Generate vector of window centers for each element of LFP
  params.window_centers = cell(size(lfp));
  for i = 1:numel(lfp)
    duration = ...
        double(lfp(i).timestamp(end) - lfp(i).timestamp(1))/TS_PER_SEC;
    num_windows = 1 + floor((duration - params.window_size) / step_size);
    % Divide the data timerange into evenly-spaced overlapping windows,
    % discarding fractional-window leftovers at the start and end
    params.window_centers{i} = lfp(i).timestamp(1) + uint32(TS_PER_SEC * ...
        ( rem(duration,step_size)/2 + ...
        params.window_size/2 + step_size*(0:(num_windows-1))' ));
  end
else
  assert(isfield(params,'window_centers'));
  % Repackage params.window_centers as a cell scalar if necessary
  if isscalar(lfp) && ~iscell(params.window_centers)
    params.window_centers = {params.window_centers};
  end
  if ~iscell(params.window_centers) || ...
      ~isequal(size(lfp),size(params.window_centers)) || ...
      ~all(cellfun(@(c) isvector(c) && isa(c,'uint32') && ...
      (size(c,1) >= 1) && (size(c,2) == 1),params.window_centers))
    error(['PARAMS.window_centers must be a cell array whose size matches ' ...
        'that of LFP, in which each cell is a column vector of uint32 ' ...
        'timestamps']);
  end
  % Assign a sentinel value to params.window_overlap
  params.window_overlap = NaN;
end

% Look up LFP samples that fall within each window and confirm that the windows
% do not exceed the range of LFP timestamps
samples_per_window = round(Fs*params.window_size);
% For each window, we want to find the LFP samples that occur at symmetric time
% offsets before and after the window center
offsets = TS_PER_SEC*((1:samples_per_window) - mean(1:samples_per_window))/Fs;
samples_idx = cell(size(lfp));
for i = 1:numel(lfp)
  samples_idx{i} = interp1(double(lfp(i).timestamp), ...
      1:numel(lfp(i).timestamp), ...
      bsxfun(@plus,double(params.window_centers{i}),offsets),'nearest');
  % We expect the LFP samples that fall around each window_center to be exactly
  % spaced at unit increments. If they are not, then this indicates that the
  % window is clipped at the start or end of the LFP record. (We already
  % checked that the timestamps are almost evenly spaced in
  % monotonically-increaseing order.)
  if ~any(all(diff(samples_idx{i},1,2) == 1,2),1)
    error(['PARAMS.window_centers and PARAMS.window_size are together ' ...
        'not compatible with the available data timerange']);
  end
end

% For a given TW, there exist 2*TW-1 Slepian tapers whose spectral
% concentration is close to optimal
num_tapers = round(2*params.TW - 1);
% Construct Slepian tapers at resolution that matches the discretization of
% spike times
[tapers, eigenweights] = dpss(samples_per_window,params.TW,num_tapers);
% normalization and reorient for later computational convenience
tapers = sqrt(Fs) .* permute(tapers,[3 1 2]);
if (num_tapers > 1)
  assert(isequal(size(tapers),[1, samples_per_window, num_tapers]));
else
  assert(isequal(size(tapers),[1, samples_per_window]));
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
    'subject'         , ...
    'day'             , ...
    'epoch'           , ...
    'environment'     , ...
    'electrode'       , ...
    'channel'         , ...
    'depth'           , ...
    'hemisphere'      , ...
    'region'          , ...
    'reference'       , ...
    'passband'        };
fn = fieldnames(lfp);
spectra = rmfield(lfp,fn(~ismember(fn,FIELDS_TO_INHERIT)));
[spectra(:).Fs] = deal(Fs);
[spectra(:).window_size] = deal(params.window_size);
[spectra(:).window_overlap] = deal(params.window_overlap);
[spectra(:).nfft] = deal(params.nfft);
[spectra(:).TW] = deal(params.TW);
[spectra(:).timestamp] = deal(params.window_centers{:});
[spectra(:).frequency] = deal(frequency(f_idx));

for i = 1:numel(spectra)
  if (prefilter)
    % Filter in double precision for numerical stability
    try
      signal = filtfilthd(Hd,double(lfp(i).samples),'reflect');
      assert(all(isfinite(signal(:))));
    catch
      error(['Could not filter, or filter application resulted in ' ...
          'numerical instability']);
    end
  else
    signal = double(lfp(i).samples);
  end
  % Compute FFT of the projections of the data onto the Slepian tapers
  J = fft(bsxfun(@times,signal(samples_idx{i}),tapers),params.nfft,2)/Fs;
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
end

if (prefilter)
  [spectra(:).filter_object] = deal(Hd);
end


