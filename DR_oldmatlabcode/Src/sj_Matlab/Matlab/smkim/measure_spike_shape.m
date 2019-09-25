function shape = measure_spike_shape(spike,interp_factor,memlimit)
%MEASURE_SHAPE Interpolate and measure width, height of spike waveforms.
%   SHAPE = MEASURE_SHAPE(SPIKE,INTERP_FACTOR,MEMLIMIT) computes
%   various measures on spike SPIKE after interpolating them by
%   INTERP_FACTOR (with piecewise cubic spline) and returns the results in a
%   struct SHAPE. This function only handles negative-going spikes. The
%   computation is done in steps that use approximately MEMLIMIT bytes at once
%   (not exact); MEMLIMIT should be set to avoid out-of-memory errors.
%
%   SPIKE must be a struct array with the following fields:
%     subject: inherited from SPIKE_CONFIG
%     day: inherited from SPIKE_CONFIG
%     tetrode: inherited from SPIKE_CONFIG
%     depth: inherited from SPIKE_CONFIG
%     region: inherited from SPIKE_CONFIG
%     reference: inherited from SPIKE_CONFIG
%     passbands: inherited from SPIKE_CONFIG
%     thresholds: inherited from SPIKE_CONFIG
%     timestamp: a Nx1 vector of uint32 (in units of 1e-4 seconds)
%     samples: 40x4xN 3-dimensional array of int16; samples(:,I,J) is the
%       vector of 40 samples recorded on the Ith channel during the Jth
%       threshold-trigger event
%     Fs: inherited from SPIKE_CONFIG
%     samples_before_trigger: number of samples that were saved before the
%       threshold-crossing; equals 8 for an NSpike *.tt file
%     amplitude_cutoff: AMPLITUDE_CUTOFF
%     sources: {FILENAME, SPIKE_CONFIG.sources{:}}
%   This struct can be produced from a NSpike *.tt file using READ_SPIKE or from a
%   Matclust file using READ_MATCLUST.
%
%   WARNING: This function will fail if the .samples of SPIKE are sign-inverted!
%   Check that the waveforms are initially negative-going; you may have to
%   multiply by -1 to undo the effect of an inverting amplifier.
%
%   INTERP_FACTOR must be a positive integer. If set to 1, no interpolation is
%   done.
%
%   SHAPE has the following fields:
%     trough_amplitude: waveform amplitude at the trough, including negative sign
%     peak_amplitude: waveform amlitude at the peak
%     half_trough_time_before: time, measured relative to the time of the first
%       sample, at which the waveform reaches half-trough before the trough
%     trough_time: time, measured relative to the time of the first sample, at
%       which the waveform reaches the trough
%     half_trough_time_after: time, measured relative to the time of the first
%       sample, at which the waveform reaches half-trough after the trough
%     zero_crossing_time: time, measured relative to the time fo the first
%       sample, at which the waveform crosses zero after the trough and before
%       the peak
%     half_peak_time_before: time, measured relative to the time of the first
%       sample, at which the waveform reaches half-peak before the peak
%     peak_time: time, measured relative to the time of the first sample, at
%       which the waveform reaches the peak
%     half_peak_time_after: time, measured relative to the time of the first
%       sample, at which the waveform reaches half-peak after the peak
%
%   References:
%
%   Henze D.A., Wittner L., Buzsaki G. (2002) Single granule cells reliably
%   discharge targets in the hippocampal CA3 network in vivo. _Nature
%   Neuroscience_ 5:790-795.
%
%   Bartho P., Hirase H., Monconduit L., Zugaro M., Harris K.D., Buzsaki G.
%   (2004) Characterization of neocortical principal cells and interneurons by
%   network interactions and extracellular features. _Journal of
%   Neurophysiology_ 92:600-608.
%
%   Berke J.D., Okatan M., Skurski J., Eichenbaum H.B. (2004) Oscillatory
%   entrainment of striatal neurons in freely moving rats. _Neuron_ 43:883-896.
%   
%Depends on:
%   IS_SPIKE (written by smk)
%   MEASURE_SPIKE_SHAPE_MEX (written by smk)
%
%Written by smk, 2009 August 25
%

if (exist('is_spike') ~= 2)
  error('MEASURE_SPIKE_SHAPE depends on m-file IS_SPIKE (written by smk)');
end
if (exist('measure_spike_shape_mex') ~= 3)
  error(['MEASURE_SPIKE_SHAPE depends on mex-file MEASURE_SPIKE_SHAPE_MEX '...
      '(written by smk)']);
end
  
BYTES_PER_SINGLE = 4; % a MATLAB single is 4 bytes

% Check that SPIKE is a valid spike waveforms data
if ~is_spike(spike)
  error('SPIKE does not appear to be a valid spike data struct array');
end

if isempty(interp_factor) || ~isscalar(interp_factor) || ...
    ~isnumeric(interp_factor) || ~isreal(interp_factor) || ...
    ~(round(interp_factor) == interp_factor) || ...
    ~(interp_factor > 0) || ~isfinite(interp_factor)
  error('INTERP_FACTOR must be a positive integer');
end

if isempty(memlimit) || ~isscalar(memlimit) || ~isnumeric(memlimit) || ...
    ~isreal(memlimit) || ~isfinite(memlimit) || (memlimit < 1e7)
  error('MEMLIMIT must an integer at least 1e7');
end

% Initialize output with meta-data inherited from SPIKE
shape = rmfield(spike, ...
    {'samples','samples_before_trigger','Fs','amplitude_cutoff'});

for i = 1:numel(spike)
  points_per_window = size(spike(i).samples,1);
  num_channels = size(spike(i).samples,2);
  num_events = size(spike(i).samples,3);
  thresholds = single(spike(i).thresholds);
  % Convert thresholds to row vector
  if (size(thresholds,2) < size(thresholds,1))
    thresholds = thresholds';
  end
  % Allocate single-precision array for each field
  tmpval = nan([num_events num_channels],'single');
  shape(i).trough_amplitude = tmpval;
  shape(i).peak_amplitude = tmpval;
  shape(i).half_trough_time_before = tmpval;
  shape(i).trough_time = tmpval;
  shape(i).half_trough_time_after = tmpval;
  shape(i).zero_crossing_time = tmpval;
  shape(i).half_peak_time_before = tmpval;
  shape(i).peak_time = tmpval;
  shape(i).half_peak_time_after = tmpval;

  % Determine the chunk size that is compatible with memlimit
  chunk_sz = floor(memlimit / ...
      (BYTES_PER_SINGLE*num_channels*points_per_window*interp_factor));
  if (chunk_sz == 0)
    error(['memlimit of %d is too small to interpolate waveforms by a ' ...
        'factor of %d'],memlimit,interp_factor);
  end
  chunk_offsets = 1:chunk_sz:num_events;
  % Time zero is defined as the time of the first (original, uninterpolated)
  % sample in the window. The times of all spike features are calculated with
  % reference to this arbitrary zero.
  t = single((0:(points_per_window-1))')/spike(i).Fs;
  if (interp_factor > 1)
    t_interp = linspace(t(1),t(end),1 + interp_factor*(numel(t)-1))';
  end
  trigger_idx = 1 + interp_factor*spike(i).samples_before_trigger;
  % Process each chunk
  for j = 1:numel(chunk_offsets)
    k = chunk_offsets(j):min(num_events,chunk_offsets(j)+chunk_sz-1);
    if (interp_factor > 1)
      % If our signal were infinitely long, we would interpolate with a sinc
      % kernel (effectively what INTERPFT does). However, the spike waveforms
      % are defined in time-limited windows which have sharp cutoffs at the
      % edges. As a result, the output of INTERPFT has artifactual ringing.
      % Instead, we use SPLINE as an acceptable approximation.
      tmp = measure_spike_shape_mex( ...
          interp1(t,single(spike(i).samples(:,:,k)),t_interp,'spline'), ...
          t_interp,trigger_idx);
    else
      tmp = measure_spike_shape_mex(single(spike(i).samples(:,:,k)), ...
          t,trigger_idx);
    end
    % MEASURE_SPIKE_SHAPE_MEX returns a struct with data fields, which we copy
    % to the appropriate slices of the matching fields of shape(i)
    names = fieldnames(tmp);   
    for l = 1:numel(names)
      shape(i) = setfield(shape,{i},names{l},{k,1:num_channels},tmp.(names{l}));
    end
  end
end


