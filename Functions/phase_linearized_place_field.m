function [occ, rate, position_bin_centers, phase_bin_centers] = phase_linearized_place_field(unit,linearized,lfp,timeranges,timestep,params)
%PHASE_LINEARIZED_PLACE_FIELD Estimate position-phase place field.
%
%   [OCC, RATE, POSITION_BIN_CENTERS, PHASE_BIN_CENTERS] = 
%   PHASE_LINEARIZED_PLACE_FIELD(UNIT,LINEARIZED,LFP,TIMERANGES,TIMESTEP,PARAMS)
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
%   LINEARIZED must be a scalar struct or a column vector struct of linearized
%   position data which validates with IS_LINEARIZED. The number of elements in
%   LINEARIZED must equal numel(UNIT). Each element of UNIT and each element of
%   LINEARIZED must share the same subject, day, epoch and environment fields.
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
%   PARAMS is a 2-dimensional struct array with the following fields:
%     PARAMS(1).TYPE must be the string 'linear', and PARAMS(2).TYPE must be the
%         string 'circular'
%     PARAMS(1).GRID is a vector of evenly-spaced grid points along the
%         linearized position dimension at which the firing rate is to be
%         estimated, expressed in the same units as linearized position, and
%         PARAMS(2).GRID is a vector of evenly-spaced grid points starting at
%         -pi and ending at +pi (phase in radians). The Cartesian product of
%         PARAMS(1).GRID and PARAMS(2).GRID defines bins in phase-position
%         space on which occupancy and firing rate are estimated.
%     PARAMS(1).DISPERSION is the standard deviation (sigma) of the Gaussian
%         smoothing kernal along the linearized position dimension, and
%         PARAMS(2).DISPERSION is the concentration parameter (kappa) of the von
%         Mises smoothing kernel along the phase dimension.
%
%   The output OCC is a 3-dimensional array of occupancy time (in seconds) per
%   bin, where each row of corresponds to a linearized position bin, each column
%   corresponds to a phase bin, and each slice OCC(:,:,i) contains occupancies
%   computed using data from the time intervals given by the cell array column
%   TIMERANGES{:,i}.
%
%   The output RATE is a 3-dimensional array of firing rate estimates (in
%   spikes/second) per bin. The dimensions of RATE correspond to the respective
%   dimensions of OCC.
%
%   POSITION_BIN_CENTERS is a vector of length size(OCC,1) of the centers of the
%   lniearized position bins. PHASE_BIN_CENTERS is a vector of length
%   size(OCC,2) of the centers of the phase bins. These are useful for plotting.
%
%Depends on:
%   IS_UNIT (written by smk)
%   IS_LINEARIZED (written by smk)
%   IS_CONTINUOUS (written by smk)
%   IS_TIMERANGE (written by smk)
%   DIFF_INTERVALS (written by smk)
%   ISMEMBER_INTERVALS (written by smk)
%   DENSITY_MAP (written by smk)
%
%Written by SMK, 2009 November 10.
%

TS_PER_SEC = 1e4;

if (exist('is_unit') ~= 2)
  error(['PHASE_LINEARIZED_PLACE_FIELD depends on m-file IS_UNIT ' ...
      '(written by smk)']);
end
if (exist('is_linearized') ~= 2)
  error(['PHASE_LINEARIZED_PLACE_FIELD depends on m-file IS_LINEARIZED ' ...
      '(written by smk)']);
end
if (exist('is_continuous') ~= 2)
  error(['PHASE_LINEARIZED_PLACE_FIELD depends on m-file IS_CONTINUOUS ' ...
      '(written by smk)']);
end
if (exist('is_timerange') ~= 2)
  error(['PHASE_LINEARIZED_PLACE_FIELD depends on m-file IS_TIMERANGE ' ...
      '(written by smk)']);
end
if (exist('diff_intervals') ~= 2)
  error(['PHASE_LINEARIZED_PLACE_FIELD depends on m-file DIFF_INTERVALS ' ...
      '(written by smk)']);
end
if (exist('ismember_intervals') ~= 2)
  error(['PHASE_LINEARIZED_PLACE_FIELD depends on m-file ISMEMBER_INTERVALS ' ...
      '(written by smk)']);
end
if (exist('density_map') ~= 2)
  error(['PHASE_LINEARIZED_PLACE_FIELD depends on m-file DENSITY_MAP ' ...
      '(written by smk)']);
end

if ~is_unit(unit) || (size(unit,2) > 1) || ...
    ~(isscalar(unit) || isvector(unit))
  error(['UNIT must be a scalar struct or column vector struct of ' ...
      'singe-unit spike data']);
end
if ~is_linearized(linearized) || (size(linearized,2) > 1) || ...
    ~(isscalar(linearized) || isvector(linearized))
  error(['LINEARIZED must be a scalar struct or column vector struct of ' ...
      'linearized position data']);
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
if ~isequal(numel(unit),numel(linearized),numel(lfp),size(timeranges,1))
  error(['number of elements in UNIT must match number of elements in ' ...
      'LINEARIZED and LFP and number of cells in TIMERANGES']);
end
for i = 1:numel(unit)
  if ~all(struct_cmp(linearized(i),unit(i), ...
      {'subject','day','epoch','environment'}))
    error(['each element of LINEARIZED must match the recording session ' ...
        'info for all elements in the matching row of UNIT']);
  end
  if ~all(struct_cmp(lfp(i),unit(i), ...
      {'subject','day','epoch','environment'}))
    error(['each element of LFP must match the recording session ' ...
        'info for all elements in the matching row of UNIT']);
  end
  if ~(linearized(i).timestamp(1) <= unit(i).timerange(1)) || ...
      ~(linearized(i).timestamp(end) >= unit(i).timerange(end))
    error(['timestamp field of each element of LINEARIZED must span the ' ...
        'timerange field of each corresponding element of UNIT']);
  end
  if ~(lfp(i).timestamp(1) <= unit(i).timerange(1)) || ...
      ~(lfp(i).timestamp(end) >= unit(i).timerange(end))
    error(['timestamp field of each element of LFP must span the ' ...
        'timerange field of each corresponding element of UNIT']);
  end
  for j = 1:size(timeranges,2)
    if ~isempty(diff_intervals(timeranges{i,j},unit(i).timerange))
      error('elements of TIMERANGES must agree with timerange fields of UNIT');
    end
  end
end

if ~isa(timestep,'uint32') || ~isscalar(timestep) || ~isreal(timestep) || ...
    ~(timestep > 0)
  error('TIMESTEP must be a real non-zero uint32 scalar');
end
if ~isstruct(params) || (numel(params) ~= 2)
  error('PARAMS must be a 2-element struct array');
end
if ~all(isfield(params,{'grid','type','dispersion'}))
  error('PARAMS is missing one or more required fields');
end
if ~strcmp('linear',params(1).type)
  error('PARAMS(1) type must be ''linear''');
end
if ~strcmp('circular',params(2).type)
  error('PARAMS(2) type must be ''circular''');
end
if (abs(params(2).grid(1) + pi) > eps(pi)) || ... 
    (abs(params(2).grid(end) - pi) > eps(pi))
  error('PARAMS(2).grid must start at -pi and end at +pi');
end
for i = 1:2
  if ~isnumeric(params(i).grid) || ~isvector(params(i).grid) || ...
     (numel(params(i).grid) < 3)
    error('PARAMS must contain numeric grid vector of at least 3 elements');
  end
  if any(diff(params(i).grid) <= 0) || ...
      any(abs(diff(diff(params(i).grid))) > 2*max(eps(params(i).grid)))
    error('grid spacings in PARAMS must be uniform monotonic increasing');
  end
  if ~(isscalar(params(i).dispersion) && isreal(params(i).dispersion) && ...
      isfinite(params(i).dispersion) && (params(i).dispersion > 0)) && ...
      ~isequal(param(i).dispersion,[])
    error(['PARAMS dispersion field must be either a real positive finite ' ...
      'scalar or the empty matrix []']);
  end
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

% Now for each set of timeranges, compute occupancy and place fields
Fs = TS_PER_SEC/double(timestep);
sz = [numel(params(1).grid)-1, numel(params(2).grid)-1, size(timeranges,2)];
occ = nan(sz);
rate = nan(sz);
% Function signature for density_map: 
%   [OCC,D] = DENSITY_MAP(X,Y,FS,PARAMS)
%
for j = 1:size(timeranges,2)
  x_concat = zeros([0,1]);
  p_concat = zeros([0,1]);
  counts_concat = zeros([0,size(unit,2)]);
  for i = 1:numel(unit)
    x_concat = [x_concat; x{i}(inclusion_mask{i,j})];
    p_concat = [p_concat; p{i}(inclusion_mask{i,j})];
    counts_concat = [counts_concat; counts{i}(inclusion_mask{i,j})];
  end
  [occ(:,:,j), rate(:,:,j)] = ...
      density_map([x_concat, p_concat],counts_concat,Fs,params);
end

% POSITION_BIN_CENTERS is a colums vector
position_bin_centers = 0.5*(params(1).grid(1:end-1) + params(1).grid(2:end));
if (size(position_bin_centers,2) > size(position_bin_centers,1))
  position_bin_centers = position_bin_centers';
end

% PHASE_BIN_CENTERS is a row vector
phase_bin_centers = 0.5*(params(2).grid(1:end-1) + params(2).grid(2:end));
if (size(phase_bin_centers,1) > size(phase_bin_centers,2))
  phase_bin_centers = phase_bin_centers';
end

