function spectra = windowed_spike_phase_spectra(unit,lfp,windows,params)
%WINDOWED_SPIKE_PHASE_SPECTRA Compute spike phase spectra of single-unit spike
%train with respect to oscillatory LFP in defined windows
%
%   SPECTRA = WINDOWED_SPIKE_PHASE_SPECTRA(UNIT,LFP,WINDOWS,PARAMS)
%   takes a struct UNIT containing single-unit spike data (of the type which
%   validates with IS_UNIT) and oscillatory LFP data containing filtered LFP
%   data (of the type which validates with IS_CONTINUOUS, and with a 'phase'
%   field that contains phase values in radians) and computes the "spike phase
%   spectrum" (Mizuseki et al., 2009) on time windows defined in the struct
%   WINDOWS using the parameters specified in PARAMS. This procedure is repeated
%   for jittered spike trains.
%
%   WINDOWS must a struct of the same size as UNIT and LFP with the following
%   field:
%     'timerange': a column cell vector in which each cell contains a 1-by-2
%         uint32 timestamp vector which specifies the start and end of the
%         window.
%   All other fields of WINDOWS must be column vectors whose number of elements
%   match the number of cells in the 'timerange' field. These will be inherited
%   as a meta-data field of the same name in the output struct SPECTRA. (Field
%   names that conflict with reserved fields of SPECTRA will raise an error.)
%   For example, if each window corresponds to a trial, then you can include
%   information about the trial type/outcome/etc.
%
%   PARAMS must be a scalar struct with the following fields:
%     'taper_family': a string specifying the family of taper; recognized values
%         are 'slepian' or 'sinusoidal'
%     'W': positive real scalar >= 1, specifying the bandwidth. Tapers are
%         selected according to the length of each window to achieve spectral
%         concentration within this bandwidth; windows of different duration may
%         have different numbers of tapers.
%     'bins_per_cycle': a real positive integer scalar, specifying the number 
%         of phase bins to discretize the LFP cycle. This should be large
%         relative to the frequency band of interest to prevent aliasing.
%     'nfft': a finite positive real integer, specifying the number of points
%         in the FFT. The computation is fastest when nfft is a power of two,
%         and slowest when nfft has large prime factors.
%     'frequency_band': a 2-element row vector containing positive real
%         floating-point values, specifying the range of frequencies (in units
%         of 1/cycle) at which the spectral estimate is to be computed.
%         frequency_band(2) must be greater than frequency_band(1) and less than
%         or equal to the Nyquist frequency (= 0.5*bins_per_cycle). An error
%         will be raised if the band is narrower than the spacing of the
%         frequency grid that is specified by nfft.
%     'num_jitters': a real positive integer, the number of jitters to perform
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
%       'lfp': struct containing various metadata fields inherited from LFP
%       'taper_family': inherited from PARAMS
%       'W': inherited from PARAMS
%       'bins_per_cycle': inherited from PARAMS
%       'nfft': inherited from PARAMS
%       'frequency': vector of real nonnegative floating-point values,
%           specifying the frequency grid of the spectral estimates, expressed
%           in units of inverse cycles 
%       'timerange': inherited from WINDOWS
%       'tapers': a cell vector whose number of elements equals the number
%           of windows. Each cell contains the Fourier transforms of the tapers
%           used to estimate the spectrum on that window.
%       'tapered_FT': a cell vector whose number of elements equals the number
%           of windows. Each cell contains a 2-dimensional array of size
%           [numel(frequency), num_tapers], containing the tapered Fourier
%           transforms of spike phase train on that window, expressed in units
%           ofspikes/sqrt(cycle). The average squared modulus of the Fourier
%           transform over all tapers is the multitaper estimate of the power
%           spectral density.
%       'jittered_tapered_FT': a cell vector of the same size as the
%           'tapered_FT' field. Each cell contains a 2-dimensional array of
%           size [numel(frequency), num_tapers, num_jitters], containing the
%           tapered Fourier transforms for the phase-jittered spike trains.
%       'spike_count': vector of the same size as the 'tapered_FT' field,
%           counting the number of spikes that contribute to the spectral
%           esimate for each window
%       'num_cycles': vector of the same size as the 'tapered_FT' field,
%           counting the number of LFP oscillation cycles that occured in each
%           window.
%       'overall_mean_rate': vector of the same size as the 'tapered_FT' field,
%           with number of spikes divided by duration (in seconds) of the window
%
%   References:
%   [1] Mizuseki K., Sirota A., Pastalkova E., Buzsaki G. (2009) Theta
%       oscillations provide temporal windows for local circuit computation in
%       the entorhinal-hippocampal loop. _Neuron_ 64:267-280.
%   [2] Geisler C., Diba K., Pastalkova E., Mizuseki K., Royer S., Buzsaki G.
%       (2010) Temporal delays among place cells determine the frequency of
%       population theta oscillations in the hippocampus. PNAS. 
%       doi: 10.1073/pnas.0912478107
%   [3] Jarvis M.R., Mitra P.P. (2001) Sampling properties of the spectrum and
%       coherency of sequences of action potentials. _Neural Computation_
%       13:717-749.
%   [4] Riedel K.S., Sidorenko A. (1995) Minimum bias multiple taper spectral
%       estimation. _IEEE Transactions on Signal Processing_ 43:188-195.
%   [5] Walden A.T., McCoy E.J., Percival D.B. (1995) The effective bandwidth of
%       a multitaper spectral estimator. _Biometrika_ 82:201-214.
%
%Depends on:
%   IS_UNIT (written by SMK)
%   IS_CONTINUOUS (written by SMK)
%   STRUCT_CMP (written by SMK)
%   IS_INTERVALS (written by SMK)
%   IS_TIMERANGE (written by SMK)
%   DIFF_INTERVALS (written by SMK)
%   UNION_INTERVALS (written by SMK)
%   INTERSECT_INTERVALS (written by SMK)
%   ISMEMBER_INTERVALS (written by SMK)
%   DPSS (MATLAB Signal Processing toolbox)
%
%Written by SMK, 2010 January 30.
%

BYTES_PER_COMPLEX_SINGLE = 8;
TS_PER_SEC = 1e4;

if (exist('is_unit') ~= 2)
  error(['WINDOWED_SPIKE_PHASE_SPECTRA depends on m-file IS_UNIT ' ...
      '(written by smk)']);
end
if (exist('is_continuous') ~= 2)
  error(['WINDOWED_SPIKE_PHASE_SPECTRA depends on m-file IS_CONTINUOUS ' ...
      '(written by smk)']);
end
if (exist('struct_cmp') ~= 2)
  error(['WINDOWED_SPIKE_PHASE_SPECTRA depends on m-file STRUCT_CMP ' ...
      '(written by smk)']);
end
if (exist('is_intervals') ~= 2)
  error(['WINDOWED_SPIKE_PHASE_SPECTRA depends on m-file IS_INTERVALS ' ...
      '(written by smk)']);
end
if (exist('is_timerange') ~= 2)
  error(['WINDOWED_SPIKE_PHASE_SPECTRA depends on m-file IS_TIMERANGE ' ...
      '(written by smk)']);
end
if (exist('diff_intervals') ~= 2)
  error(['WINDOWED_SPIKE_PHASE_SPECTRA depends on m-file DIFF_INTERVALS ' ...
      '(written by smk)']);
end
if (exist('union_intervals') ~= 2)
  error(['WINDOWED_SPIKE_PHASE_SPECTRA depends on m-file UNION_INTERVALS ' ...
      '(written by smk)']);
end
if (exist('intersect_intervals') ~= 2)
  error(['WINDOWED_SPIKE_PHASE_SPECTRA depends on m-file INTERSECT_INTERVALS ' ...
      '(written by smk)']);
end
if (exist('ismember_intervals') ~= 2)
  error(['WINDOWED_SPIKE_PHASE_SPECTRA depends on m-file ISMEMBER_INTERVALS ' ...
      '(written by smk)']);
end
if (exist('dpss') ~= 2)
  error(['WINDOWED_SPIKE_PHASE_SPECTRA depends on m-file DPSS ' ...
      '(MATLAB Signal Processing toolbox)']);
end

if ~is_unit(unit)
  error('UNIT is not a valid struct containing single-unit spike data');
end
if ~is_continuous(lfp) || ~isfield(lfp,{'phase'}) || ...
    ~all(arrayfun(@(l) isvector(l.phase) && isfloat(l.phase) && ...
    isreal(l.phase) && all((l.phase >= -pi) & (l.phase <= +pi)),lfp))
  error(['LFP is not a valid struct containing continuous LFP data with ' ...
      'an instantaneous phase field']);
end
if ~isequal(size(unit),size(lfp),size(windows))
  error('UNIT, LFP, and WINDOWS must be the same size');
end
if ~all(struct_cmp(lfp(1),lfp, ...
    {'subject','electrode','reference','depth','region','hemisphere', ...
    'Fs','environment'}))
  error(['Elements of LFP must correspond to the same single electrode ' ...
      'and environment and must have the same sampling rate and filter ' ...
      'specifications']);
end
if ~isequal(size(unit),size(lfp)) || ~all(arrayfun(@(u,l) ...
    struct_cmp(u,l,{'subject','day','epoch','environment','hemisphere'}), ...
    unit,lfp))
  error(['UNIT and LFP must match in size and must have ' ...
      'matching ''subject'',''day'',''epoch'',''environment'', ' ...
      '''hemisphere'' metadata']);
end
if ~all(arrayfun(@(u,l) isempty(diff_intervals( ...
    u.timerange,[l.timestamp(1), l.timestamp(end)])),unit,lfp))
  error(['Timestamps of LFP samples in LFP must span the timerange of ' ...
      'spike data in the corresponding element of UNIT']);
end

if ~isstruct(windows) || ~isfield(windows,'timerange') || ...
    ~all(arrayfun(@(x) iscell(x.timerange) && ...
    all(cellfun(@(c) is_timerange(c) && (size(c,1)==1),x.timerange)),windows))
  error('''timerange'' field of WINDOWS is missing or invalid');
end
if ~all(arrayfun(@(u,w) isempty(w.timerange) || ...
    (isequal(union_intervals(w.timerange{:}), ...
    intersect_intervals(u.timerange,union_intervals(w.timerange{:}))) && ...
    isequal(u.timerange, ...
    union_intervals(u.timerange,union_intervals(w.timerange{:})))), ...
    unit,windows))
  error(['WINDOWS lie outside the timerange of the corresponding ' ...
      'element of UNIT']);
end
WINDOW_METADATA_FIELDS = fieldnames(windows);
WINDOW_METADATA_FIELDS(strcmp('timerange',WINDOW_METADATA_FIELDS)) = [];
if ~all(cellfun(@(c) all(arrayfun(@(x) isvector(x.(c)) && ...
    isequal(size(x.(c)),size(x.timerange)),windows(:))), ...
    WINDOW_METADATA_FIELDS))
  error(['Optional metadata fields of WINDOWS must be column vectors ' ...
      'whose lengths match the number of rows in ''timerange''']);
end

if ~isstruct(params) || ~isscalar(params) || ~all(isfield(params, ...
    {'taper_family','W','bins_per_cycle','nfft','frequency_band', ...
    'num_jitters'}))
  error('PARAMS is not a scalar struct with the required fields');
end

if ~ischar(params.taper_family) || ...
    ~any(strcmp(params.taper_family,{'slepian','sinusoidal'}))
  error(['PARAMS.taper_family is not a recognized string value ' ...
      '(either ''slepian'' or ''sinusoidal''']);
end
if ~isscalar(params.W) || ~isnumeric(params.W) || ~isreal(params.W) || ...
    ~isfinite(params.W) || ~(params.W > 0)
  error('PARAMS.W must be a real finite positive scalar');
end
if ~isscalar(params.bins_per_cycle) || ~isreal(params.bins_per_cycle) || ...
    ~isfinite(params.bins_per_cycle) || ~(params.bins_per_cycle > 0) || ...
    (params.bins_per_cycle ~= round(params.bins_per_cycle))
  error(['PARAMS.bins_per_cycle must be a real positive integer scalar']);
end
if ~isscalar(params.nfft) || ~isreal(params.nfft) || ...
    ~isfinite(params.nfft) || (params.nfft < 2) || ...
    ~(round(params.nfft) == params.nfft)
  error('PARAMS.nfft must be a real finite integer-valued scalar >1');
end
if ~is_intervals(params.frequency_band) || ...
    ~isequal(size(params.frequency_band,1),1) || ...
    ~all(params.frequency_band >= 0) || ...
    ~all(params.frequency_band < params.bins_per_cycle/2)
  error(['PARAMS.frequency_band must be a 2-element row vector which ' ...
      'specifies a frequency interval between 0 and Nyquist']);
end
if ~isscalar(params.num_jitters) || ~isreal(params.num_jitters) || ...
    ~isfinite(params.num_jitters) || ~(params.num_jitters > 0) || ...
    (round(params.num_jitters) ~= params.num_jitters)
  error(['PARAMS.num_jitters must be a positive real integer']);
end
    
% FFT frequency grid
frequency = linspace(0,params.bins_per_cycle,1+params.nfft)';
frequency(end) = [];
% subset of the FFT frequencies that span the desired frequency band
f_idx = find(ismember_intervals(frequency,params.frequency_band));
if (f_idx(1) > 1)
  f_idx = [f_idx(1)-1; f_idx];
end
if (f_idx(end) < numel(frequency))
  f_idx = [f_idx; f_idx(end)+1];
end

% Allocate output struct
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
LFP_FIELDS_TO_INHERIT = { ...
    'electrode'       , ...
    'channel'         , ...
    'depth'           , ...
    'hemisphere'      , ...
    'region'          , ...
    'reference'       , ...
    'passband'        , ...
    'Fs'              , ...
    'filter_object'   };
lfp_fn = fieldnames(lfp);
for i = 1:numel(spectra)
  spectra(i).lfp = rmfield(lfp(i), ...
      lfp_fn(~ismember(lfp_fn,LFP_FIELDS_TO_INHERIT)));
end
[spectra(:).taper_family] = deal(params.taper_family);
[spectra(:).W] = deal(params.W);
[spectra(:).bins_per_cycle] = deal(params.bins_per_cycle);
[spectra(:).nfft] = deal(params.nfft);
[spectra(:).frequency] = deal(frequency(f_idx));
[spectra(:).timerange] = deal(windows(:).timerange);
% Reserve required fields so that we can check for namespace collision with
% fields inherited from WINDOWS
[spectra(:).spike_count] = deal([]);
[spectra(:).num_cycles] = deal([]);
[spectra(:).overall_mean_rate] = deal([]);
[spectra(:).tapers] = deal([]);
[spectra(:).tapered_FT] = deal([]);
[spectra(:).jittered_tapered_FT] = deal([]);
% Inherit optional fields from WINDOWS
if any(ismember(WINDOW_METADATA_FIELDS,fieldnames(spectra)))
  error(['Optional metadata field(s) of WINDOWS has naming conflict with ' ...
      'required fields of SPECTRA']);
end
for f = 1:numel(WINDOW_METADATA_FIELDS)
  [spectra(:).(WINDOW_METADATA_FIELDS{f})] = ...
      deal(windows(:).(WINDOW_METADATA_FIELDS{f}));
end

for i = 1:numel(spectra)
  % Compute unwrapped oscillatory phase for each LFP sample
  lfp_cycle = unwrap(lfp(i).phase)/(2*pi);
  % For an ideal monocomponent oscillation, the unwrapped phase is strictly
  % monotonically increasing. But for real data, the unwrapped phase often has
  % local defects (due to real phase slips in the oscillation as well as noise
  % in the phase estimation). We report the proportion of these phase defects to
  % the user.
  monotonicity_violation_fraction = ...
      nnz(diff(lfp_cycle) <= 0)/(numel(lfp_cycle)-1);
  if (monotonicity_violation_fraction > 0)
    warning(['%.3f%% of consecutive LFP samples do not exhibit strictly ' ... 
        'monotonically increasing phase'],100*monotonicity_violation_fraction);
  end
  % Get unwrapped LFP cycle at the beginning and end of each window
  num_windows = numel(windows(i).timerange);
  window_cycle = interp1(double(lfp(i).timestamp),lfp_cycle, ...
      double(cell2mat(windows(i).timerange)));
  % For each window, construct bin edges for discretizing spike phases and
  % tapers with matching support
  bin_edges = cell([num_windows, 1]);
  tapers = cell([num_windows, 1]);
  ntapers = zeros([num_windows, 1]);
  for j = 1:num_windows
    T = window_cycle(j,end) - window_cycle(j,1);
    nbins = round(params.bins_per_cycle*T);
    if (params.nfft < nbins)
      error(['PARAMS.nfft is smaller than the number of discrete bins ' ...
          'in the binned spike phase train.']);
    end
    bin_edges{j} = linspace(window_cycle(j,1),window_cycle(j,end),1+nbins)';
    % Construct tapers. Bandwidth formulas are from Walden et al. (1995).
    switch params.taper_family
    case 'slepian'
      % Construct Slepian tapers at resolution that matches the discretization of
      % spike times. For a given TW, there exist 2*TW-1 Slepian tapers whose
      % spectral concentration is close to optimal
      ntapers(j) = floor(2*T*params.W - 1);
      if (ntapers(j) < 2)
        error('params.W is too small for this window size.');
      end
      % Normalize tapers to unit impulse and convert to single precision.
      tapers{j} = single(sqrt(params.bins_per_cycle) .* ...
          dpss(nbins,T*params.W,ntapers(j)));
    case 'sinusoidal'
      % Sinusoidal tapers are close to optimal for minimizing local bias within
      % bandwidth (Riedel & Sidorenko, 1995). Convert to single precision.
      ntapers(j) = floor(2*T*params.W*(nbins+1)/nbins - 1);
      if (ntapers(j) < 2)
        error('params.W is too small for this window size.');
      end
      tapers{j} = zeros([nbins, ntapers(j)],'single');
      for k = 1:ntapers(j)
        tapers{j}(:,k) = sqrt(params.bins_per_cycle) * ...
            sqrt(2/(nbins+1)) * sin(pi*k*(1:nbins)/(nbins+1));
      end
    otherwise
      error('PARAMS.taper_family is not a recognized value');
    end
    clear('nbins');
  end
  assert(isequal(ntapers,cellfun(@(c) size(c,2),tapers)));
  % spectra(i).tapers contains the *Fourier transforms* of tapers
  spectra(i).tapers = cellfun(@(c) fft(c,params.nfft,1),tapers, ...
      'UniformOutput',false);
  assert(all(cellfun(@(c) isa(c,'single'),tapers)) && ...
      all(cellfun(@(c) isa(c,'single'),spectra(i).tapers)));
  % Interpolate instantaneous phase at each spike time
  spike_cycle = interp1(double(lfp(i).timestamp),lfp_cycle, ...
      double(unit(i).timestamp));
  % Jitter the phase of each spike independently by a random offset that is
  % drawn from uniform distribution on [-pi,+pi]. This effectively jitters the
  % phase *intervals* between consecutive spikes with a triangle distribution
  % centered around zero and supported on [-2*pi,+2*pi].
  jittered_spike_cycle = bsxfun(@plus,spike_cycle, ...
      0.5-rand([numel(spike_cycle), params.num_jitters]));
  % Bin the spikes that fall within each window. spike_count{j} is the histogram
  % of spikes that fall within the jth window. Likewise,
  % jittered_spike_count{j}(:,1,k) is the histogram of jittered spikes from the
  % kth jittered spike train that fall within the jth window. We assume that
  % these counts can be contained in uint8.
  spike_count = cell([num_windows, 1]);
  jittered_spike_count = cell([num_windows, 1]);
  for j = 1:num_windows
    spike_count{j} = zeros([numel(bin_edges{j})-1, 1],'uint8');
    tmp = histc(spike_cycle(ismember_intervals( ...
        spike_cycle,window_cycle(j,:))),bin_edges{j},1);
    spike_count{j} = tmp(1:end-1);
    jittered_spike_count{j} = ...
        zeros([numel(bin_edges{j})-1, 1, params.num_jitters],'uint8');
    for k = 1:params.num_jitters
      tmp = histc(jittered_spike_cycle(ismember_intervals( ...
          jittered_spike_cycle(:,k),window_cycle(j,:)),k),bin_edges{j},1);
      % The second dimension is singleton, and the third dimension corresponds
      % to jitters. This array shape is for computational convenience with
      % BSXFUN.
      jittered_spike_count{j}(:,1,k) = tmp(1:end-1);
    end
  end
  % Check that we didn't overflow uint8
  assert(all(cellfun(@(c) ~any(c(:) == intmax('uint8')),spike_count)) && ...
      all(all(cellfun(@(c) ~any(c(:) == intmax('uint8')), ...
      jittered_spike_count),1),2));
  % Count total number of spikes on each window 
  spectra(i).spike_count = cellfun(@sum,spike_count);
  % Number of LFP cycles spanned by each window 
  spectra(i).num_cycles = cellfun(@range,bin_edges);
  % Overall mean firing rate
  spectra(i).overall_mean_rate = cellfun(...
      @(c1,c2) sum(c1)/(double(c2(end) - c2(1))/TS_PER_SEC), ...
      spike_count,windows.timerange);
  % Clear some variables that we no longer need to free memory
  clear('lfp_cycle','spike_cycle','jittered_spike_cycle','window_cycle', ...
      'bin_edges');
  % To save memory, we compute in single precision and save the FFT only over
  % the frequencies of interest
  spectra(i).tapered_FT = cell([num_windows, 1]);
  spectra(i).jittered_tapered_FT = cell([num_windows, 1]);
  for j = 1:num_windows
    % Progress meter
    status_msg = sprintf('processing window %d out of %d',j,num_windows);
    fprintf(1,'%s',status_msg);
    pause(0.001);
    fprintf(repmat('\b',size(status_msg)));
    % No need to compute if there are no spikes
    if (spectra(i).spike_count == 0)
      spectra(i).tapered_FT{j} = ...
          zeros([numel(f_idx), ntapers(j)],'single');
      spectra(i).jittered_tapered_FT{j} = ...
          zeros([numel(f_idx), ntapers(j), params.num_jitters],'single');
      continue;
    end
    tmp = single(spike_count{j});
    % spectra(i).tapers{j} is the FT of tapers{j}
    tmp = fft(bsxfun(@times,tmp,tapers{j}),params.nfft,1) - ...
        bsxfun(@times,mean(tmp,1),spectra(i).tapers{j});
    spectra(i).tapered_FT{j} = tmp(f_idx,:);
    clear('tmp');

    % Vectorized operation is faster, but we loop over jitters to avoid
    % exceeding available memory
    spectra(i).jittered_tapered_FT{j} = ...
        zeros([numel(f_idx), ntapers(j), params.num_jitters],'single');
    for k = 1:params.num_jitters
      tmp = single(jittered_spike_count{j}(:,:,k));
      % spectra(i).tapers{j} is the FT of tapers{j}
      tmp = fft(bsxfun(@times,tmp,tapers{j}),params.nfft,1) - ...
          bsxfun(@times,mean(tmp,1),spectra(i).tapers{j});
      spectra(i).jittered_tapered_FT{j}(:,:,k) = tmp(f_idx,:,:);
      clear('tmp');
    end
  end
  pause(0.001);
  fprintf(1,'\n');
  assert(all(cellfun(@(c) size(c,1)==numel(spectra(i).frequency), ...
      spectra(i).tapered_FT)) && ...
      all(cellfun(@(c) size(c,1)==numel(spectra(i).frequency), ...
      spectra(i).jittered_tapered_FT)))
  assert(all(cellfun(@(c) size(c,3)==params.num_jitters, ...
      spectra(i).jittered_tapered_FT)));
end

