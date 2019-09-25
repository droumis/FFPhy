function ratemap = phase_rate_modulation(unit,lfp,timerange,params)
%PHASE_RATE_MODULATION Estimate oscillatory LFP modulation of single-unit firing rate
%
%   RATEMAP = PHASE_RATE_MODULATION(UNIT,LFP,TIMERANGES,PARAMS) produces
%   a phase-dependent firing-rate estimate(s), given a single-unit spike data
%   scalar struct UNIT and matching LFP data LFP. Separate estimates are made
%   for the different sets of time intervals that are specified in the cell
%   array TIMERANGES.
%
%   UNIT must be a scalar or vector struct of single-unit spike data which
%   validates with IS_UNIT. All elements of UNIT must share the same
%   'subject' and 'environment' fields.
%
%   LFP must be a scalar or vector struct of continuus LFP data which validates
%   with IS_CONTNUOUS, whose subject, day, epoch, and environment fields match
%   the corresponding elements of UNIT. All elements of LFP must share the same
%   'subject', 'environment', and 'units' fields. LFP must have a 'phase' field.
%   
%   TIMERANGES must be a cell vector of the same length as UNIT and LFP. Each
%   cell must contain a Nx2 matrix of uint32 timestamp intervals (of the type
%   which validates with IS_TIMERANGE), and these intervals must fall entirely
%   within the time ranges of UNIT and LINEARIZED. If UNIT and LFP are scalar
%   structs, then TIMERANGES may be an Nx2 uint32 matrix without the packaging
%   cell.
%
%   PARAMS is a scalar struct with the following fields:
%     'discretization_timestep': must be a uint32 scalar, which specifies the
%         time discretization resolution. Smaller values yield better estimates
%         at greater computational cost.
%     'phase_bin_edges': a vector of evenly-spaced grid points along the phase
%         dimension at which the firing rate is to be estimated, expressed in
%         radians.
%     'phase_kappa' is the (unitless) concentration parameter of the von Mises
%         smoothing kernel. Enter empty matrix [] for no smoothing.
%
%   The output RATEMAP is a scalar struct array. RATEMAP has the following
%   fields:
%     'subject': inherited from UNIT and LFP
%     'environment': inherited from UNIT and LFP
%     'days': inherited from UNIT and LFP
%     'epochs': inherited from UNIT and LFP
%     'timerange': inherited from TIMERANGES
%     'discretization_timestep': inherited from PARAMS
%     'phase_kappa': inherited from PARAMS
%     'phase_bin_centers': a vector of bin centers, determined from
%         PARAMS.phase_bin_edges
%     'rate': firing rate matrix; each row corresponds to a phase bin
%     'occupancy': occupancy matrix; each row corresponds to a phase bin
%     'unit': struct containing various metadata inherited from UNIT
%     'lfp': struct containing various metadata inherited from LFP
%
%Depends on:
%   IS_UNIT (written by smk)
%   IS_CONTINUOUS (written by smk)
%   STRUCT_CMP (written by SMK)
%   IS_TIMERANGE (written by smk)
%   DIFF_INTERVALS (written by smk)
%   ISMEMBER_INTERVALS (written by smk)
%   DENSITY_MAP (written by smk)
%
%Written by SMK, 2009 November 10.
%

TS_PER_SEC = 1e4;

if (exist('is_unit') ~= 2)
  error(['LINEARIZED_PLACE_FIELD depends on m-file IS_UNIT ' ...
      '(written by smk)']);
end
if (exist('is_continuous') ~= 2)
  error(['LINEARIZED_PLACE_FIELD depends on m-file IS_CONTINUOUS ' ...
      '(written by smk)']);
end
if (exist('struct_cmp') ~= 2)
  error(['LINEARIZED_PLACE_FIELD depends on m-file STRUCT_CMP ' ...
      '(written by smk)']);
end
if (exist('is_timerange') ~= 2)
  error(['LINEARIZED_PLACE_FIELD depends on m-file IS_TIMERANGE ' ...
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
if (exist('density_map') ~= 2)
  error(['LINEARIZED_PLACE_FIELD depends on m-file DENSITY_MAP ' ...
      '(written by smk)']);
end

if ~is_unit(unit) || ~(isscalar(unit) || isvector(unit))
  error(['UNIT must be a scalar or vector struct of singe-unit spike data']);
end
if ~all(struct_cmp(unit(1),unit, ...
    {'subject','uid','tetrode','reference','depth','region','hemisphere', ...
    'clustnum','environment'}))
  error(['Elements of UNIT must correspond to the same single unit and ' ...
      'environment']);
end
if ~is_continuous(lfp) || ~(isscalar(lfp) || isvector(lfp)) || ...
    ~isfield(lfp,'phase')
  error(['LFP must be a scalar or vector struct of continuous LFP data ' ...
        'with a ''phase'' field']);
end
if ~all(struct_cmp(lfp(1),lfp, ...
    {'subject','electrode','reference','depth','region','hemisphere', ...
    'Fs','environment'}))
  error(['Elements of LFP must correspond to the same single electrode ' ...
      'and environment and must have the same sampling rate and filter ' ...
      'specifications']);
end
if ~isequal(size(unit),size(lfp))
  error('UNIT and LFP must match in size');
end
if ~all(arrayfun(@(u,l) struct_cmp(u,l, ...
    {'subject','day','epoch','environment'}),unit,lfp))
  error(['Corresponding elements of UNIT and LFP must have ' ...
      'matching subject, day, epoch, environment metadata']);
end

if is_timerange(timerange) 
  timerange = {timerange};
end
if ~iscell(timerange) || ~all(cellfun(@is_timerange,timerange))
  error(['TIMERANGES must be a Nx2 matrix of uint32 timestamp intervals, ' ...
      'or a cell array of such timestamp interval matrices']);
end
if ~isvector(timerange) || (numel(timerange) ~= numel(unit))
  error(['Number of cells in TIMERANGES must match number of elements in ' ...
      'UNIT and LINEARIZED']);
end
if (size(timerange,1) < size(timerange,2))
  timerange = timerange';
end
for i = 1:numel(timerange)
  if ~isempty(diff_intervals(timerange{i},unit(i).timerange)) || ...
      ~isempty(diff_intervals(timerange{i}, ...
      [lfp(i).timestamp(1), lfp(i).timestamp(end)]))
    error(['All time intervals in TIMERANGES must fall within time range ' ...
        'of corresponding elements of UNIT, LINEARIZED and LFP']);
  end
end

if ~isstruct(params) || ~isscalar(params)
  error('PARAMS must be a scalar struct array');
end
if ~all(isfield(params, ...
    {'discretization_timestep','phase_bin_edges', 'phase_kappa'}));
  error('PARAMS is missing one or more required fields');
end
if ~isscalar(params.discretization_timestep) || ...
    ~isa(params.discretization_timestep,'uint32') || ...
    ~(params.discretization_timestep > 0)
  error(['''discretization_timestep'' field of PARAMS must be non-zero ' ...
      'uint32 scalar']);
end
if ~isnumeric(params.phase_bin_edges) || ...
    ~isvector(params.phase_bin_edges) || ...
    (numel(params.phase_bin_edges) < 3) || ...
    any((params.phase_bin_edges < -pi) | (params.phase_bin_edges > +pi))
  error(['''phase_bin_edges'' field of PARAMS must be a ' ...
      'numeric grid vector of least 3 elements, with values between ' ...
      '-pi and +pi']);
end
if any(diff(params.phase_bin_edges) <= 0) || ...
    any(abs(diff(diff(params.phase_bin_edges))) > ...
    sqrt(max(eps(params.phase_bin_edges))))
  error(['phase_bin_edges spacings in PARAMS must be ' ...
      'uniform monotonic increasing']);
end
if ~(isscalar(params.phase_kappa) && ...
    isreal(params.phase_kappa) && ...
    isfinite(params.phase_kappa) && ...
    (params.phase_kappa > 0)) && ...
    ~isequal(params.phase_kappa,[])
  error(['''phase_kappa'' field of PARAMS must be either ' ...
      'a real positive finite scalar or the empty matrix []']);
end

% Data preparation
include = cell([numel(timerange), 1]);
for i = 1:numel(timerange)
  % Discretize all data across the entire session
  t{i,1} = double(unit(i).timerange(1): ...
      params.discretization_timestep:unit(i).timerange(end))';
  % Select only timebins that fall within each element of timerange
  include{i,1} = find(ismember_intervals(uint32(t{i}(1:end-1)),timerange{i}));
  % Spike counts
  counts{i,1} = histc(double(unit(i).timestamp),t{i});
  % histc returns an extra element at the end which we don't need
  counts{i}(end) = [];
  % Resampled linearized position and lfp phase
  x{i,1} = interp1_angle(double(lfp(i).timestamp), ...
      lfp(i).phase,t{i}(1:end-1),pi);
end

% Concatenate across sessions and compute occupancy and place fields
counts = cell2mat(cellfun(@(arg1,arg2) arg1(arg2,:),counts,include, ...
    'UniformOutput',false));
x = cell2mat(cellfun(@(arg1,arg2) arg1(arg2,:),x,include, ...
    'UniformOutput',false));
Fs = TS_PER_SEC/double(params.discretization_timestep);
density_map_params = struct( ...
    'type',{'circular'}, ...
    'grid',{params.phase_bin_edges}, ...
    'dispersion',{params.phase_kappa});
[ratemap.occupancy, ratemap.rate] = ...
    density_map(x,counts,Fs,density_map_params);

% Add metadata
phase_bin_centers = 0.5 * ...
    (params.phase_bin_edges(1:end-1) + ...
    params.phase_bin_edges(2:end));
if (size(phase_bin_centers,1) > size(phase_bin_centers,2))
  phase_bin_centers = phase_bin_centers';
end
ratemap.phase_bin_centers = phase_bin_centers;
ratemap.phase_kappa = params.phase_kappa;
ratemap.discretization_timestep = params.discretization_timestep;
ratemap.timerange = timerange;
ratemap.subject = unit(1).subject;
ratemap.environment = unit(1).environment;
ratemap.days = {unit(:).day}';
ratemap.epochs = {unit(:).epoch}';
ratemap.unit = rmfield(unit, ...
    {'sources','subject','day','epoch','environment','passbands', ...
    'thresholds','Fs','samples_before_trigger','amplitude_cutoff', ...
    'timestamp','samples','timerange','maximum_amplitude_channel', ...
    'shape','isi_pre','isi_post','burst'});
ratemap.lfp = rmfield(lfp, ...
    {'sources','subject','day','epoch','environment','Fs', ...
    'timestamp','samples','phase'});


