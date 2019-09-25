function [S,f] = spike_phase_spectrum(unit,lfp,timeranges,params)
%SPIKE_PHASE_SPECTRUM Estimate point-process spectrum with respect to oscillatory phase.
%
%   [S, F] = SPIKE_PHASE_SPECTRUM(UNIT,LFP,TIMERANGES,PARAMS) produces a
%   multitaper estimate of the "spike phase spectrum" (Mizuseki et al., 2009) of
%   a single-unit spiking relative to LFP oscillation. This measure is similar
%   to a conventional point-process spectrum, except that the spike times are
%   rescaled to match fluctuations in the LFP oscillation period and the
%   frequency axis is in units of 1/cycles instaed of 1/seconds.
%
%   UNIT must be a column vector struct of single-unit spike data which
%   validates with IS_UNIT. All elements of UNIT must correspond to the same
%   single unit and share the same uid field.
%
%   LFP must be a scalar struct or column vector struct of continuous data which
%   validates with IS_CONTINUOUS. Furthermore, the 'samples' field of LFP must
%   be complex floating-point, such that phase can be extracted with the ANGLE
%   function. This phase must be strictly monotonically increasing after mod
%   (2*pi) unwrapping! The number of elements in LFP must equal numel(UNIT).
%   Each element of UNIT and each element of LFP must share the same subject,
%   day, epoch and environment fields.
%
%   [OCC, RATE, POSITION_BIN_CENTERS, PHASE_BIN_CENTERS] = 
%   POSITION_PHASE_PLACE_FIELD(UNIT,LINEARIZED,LFP,TIMERANGES,TIMESTEP,PARAMS)
%   produces two-dimensional place field estimate(s) on the phase plane defined
%   by linearized position and LPF oscillation phase, given matching single-unit
%   spike data UNIT, linearized position data LINEARIZED, and filtered
%   Hilbert-transformed continuous data LFP. Separate estimates are made for the
%   different sets of time interval that specified by TIMERANGES.
%
%   UNIT must be a column vector struct of single-unit spike data which
%   validates with IS_UNIT. All elements of UNIT must correspond to the same
%   single unit and share the same uid field.
%
%   LFP must be a scalar struct or column vector struct of continuous data which
%   validates with IS_CONTINUOUS. Furthermore, the 'samples' field of LFP must
%   be complex floating-point, such that phase can be extracted with the ANGLE
%   function. The number of elements in LFP must equal numel(UNIT). Each element
%   of UNIT and each element of LFP must share the same subject, day, epoch and
%   environment fields.
%   
%   TIMERANGES must be a 2-dimensional cell array with numel(UNIT) rows. The
%   contents of TIMERANGES{i,j} must fall entirely within the time interval
%   given by UNIT(i,:).timerange. Each column of TIMERANGES is a set of time
%   filters that contributes to corresponding slices of OCC, such that
%   size(OCC,2) == size(RATE,2) == size(TIMERANGES,2).
%
%   TIMESTEP is the time bin on which the data are to be discretized, expressed
%   in uint32 timestamp units. Smaller values yield better estimates at greater
%   computational cost.
%
%   PARAMS is a struct with the following fields:
%     fpass: 2-element vector [fmin fmax] of the range of frequencies of
%       interest, expressed in 1/cycles
%     pad: padding factor for the FFT. -1 means no padding, 0 means padding to
%       the smallest power of 2 that accomodate the number of samples, +1
%       corresponds to the next highest power of 2, etc.
%     tapers: a three-element vector [W T p] where W is the bandwidth, T is the
%       window size (in the same units as sampling rate Fs), and p is an integer
%       such that (2*T*W - p) tapers are used.
%   Note: Do not be confused by the superficial similarity to a Chronux-style
%   parameter struct. An error will be raised if the PARAMS struct contains any
%   fields other than the ones defined above! 
%
%TODO: finite size correction. jackknife confidence bounds.
%
%   References:
%
%   Jarvisa M.R., Mitra P.P. (2001) Sampling properties of the spectrum and
%   coherency of sequences of action potential. _Neural Computation_ 13:717-749.
%
%   Mizuseki K., Sirota A., Pastalkova E., Buzsaki G. (2009) Theta oscillations
%   provide temporal windows for local circuit computation in the
%   entorhinal-hippocampal loop. _Neuron_ 64:267-280.
%
%Depends on:
%   IS_UNIT (written by smk)
%   IS_CONTINUOUS (written by smk)
%   DIFF_INTERVALS (written by SMK)
%   ISMEMBER_INTERVALS (written by smk)
%     (various Chronux functions)
%
%Written by SMK, 2009 November 21.
%

   
if (exist('is_unit') ~= 2)
  error(['POSITION_PHASE_PLACE_FIELD depends on m-file IS_UNIT ' ...
      '(written by smk)']);
end
if (exist('is_continuous') ~= 2)
  error(['POSITION_PHASE_PLACE_FIELD depends on m-file IS_CONTINUOUS ' ...
      '(written by smk)']);
end
if (exist('diff_intervals') ~= 2)
  error(['POSITION_PHASE_PLACE_FIELD depends on m-file DIFF_INTERVALS ' ...
      '(written by smk)']);
end
if (exist('ismember_intervals') ~= 2)
  error(['POSITION_PHASE_PLACE_FIELD depends on m-file ISMEMBER_INTERVALS ' ...
      '(written by smk)']);
end

if ~is_unit(unit) || (size(unit,2) > 1) || ...
    ~(isscalar(unit) || isvector(unit))
  error(['UNIT must be a scalar struct or column vector struct of ' ...
      'singe-unit spike data']);
end
if ~is_continuous(lfp) || (size(lfp,2) > 1) || ...
    ~(isscalar(lfp) || isvector(lfp)) || ...
    ~all(cellfun(@(c) isfloat(c) && ~isreal(c),{lfp.samples}))
  error(['LFP must be a scalar struct or column vector struct of Hilbert-' ...
      'transformed (i.e. complex-valued) bandpass-filtered LFP data']);
end
if ~iscell(timeranges) || (ndims(timeranges) > 2) || ...
    ~all(cellfun(@is_timerange,timeranges(:)))
  error('TIMERANGES must be a cell array of Nx2 uint32 timerange arrays');
end
if ~isequal(numel(unit),numel(lfp),size(timeranges,1))
  error(['number of elements in UNIT must match number of elements in ' ...
      'LFP and number of cells in TIMERANGES']);
end
for i = 1:numel(unit)
  if ~all(struct_cmp(lfp(i),unit(i), ...
      {'subject','day','epoch','environment'}))
    error(['each element of LFP must match the recording session ' ...
        'info for all elements in the matching row of UNIT']);
  end
  if ~(lfp(i).timestamp(1) <= unit(i).timerange(1)) || ...
      ~(lfp(i).timestamp(end) >= unit(i).timerange(end))
    error(['timestamp field of each element of LFP must span the ' ...
        'timerange field of each corresponding element of UNIT']);
  end
  for j = 1:size(timeranges,2)(timeranges,2)
    if ~isempty(diff_intervals(timeranges{i,j},unit(i).timerange))
      error('elements of TIMERANGES must agree with timerange fields of UNIT');
    end
  end
end

% Data preparation
for i = 1:numel(unit)
  unwrapped_phase{i} = unwrap(angle(lfp(i).samples));
  cycle{i} = interp(double(lfp(i).timestamp),unwrapped_phase{i}, ...
      double(unit(i).timestamp)) / (2*pi);
  % Map timeranges to unwrapped phase ranges
  phase_ranges{i} = 

interp(unwrapped_phase{i},double(lfp(i).timestamp), ...
      double(
    mtspectrumpt(unwrapped_p,params,1,timegrid)
  % what is a correct time grid? are these window starttimes or centers?

end
 
% Check for monotonicity of phase advance and warn if it exhibits slips
phase_slip_fraction = ...
    sum(cellfun(@(p) nnz(diff_angle(p) <= 0),unwrapped_phase)) ./ ...
    sum(cellfun(@(p) numel(p) - 1,unwrapped_phase));
if (phase_slip_fraction > 0)
  warning('%f of LFP samples exhibit backwards phase slip', ...
      phase_slip_fraction);
end

% Data preparation
inclusion_mask = cell([numel(unit) size(timeranges,2)]);
for i = 1:numel(unit)
  % First, discretize all data across the entire session
  timerange = unit(i,1).timerange;
  t{i} = double(timerange(1):timestep:timerange(end))';
  % Spike counts
  counts{i} = histc(double(unit(i).timestamp),t{i});
  % histc returns an extra element at the end which we don't need
  counts{i}(end) = [];
  % Resampled linearized position
  x{i} = interp1(double(linearized(i).timestamp), ...
      linearized(i).linearized_position,t{i}(1:end-1));
  % Resampled LFP oscillation phase
  p{i} = interp1_angle(double(lfp(i).timestamp),angle(lfp(i).samples), ...
      t{i}(1:end-1),pi);
  % Select only timebins that fall within timeranges{i,j}
  for j = 1:size(timeranges,2)
    inclusion_mask{i,j} = ...
        ismember_intervals(uint32(t{i}(1:end-1)),timeranges{i,j});
  end 
end



