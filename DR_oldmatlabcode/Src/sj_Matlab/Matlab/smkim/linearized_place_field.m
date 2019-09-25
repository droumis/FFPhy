function ratemap = linearized_place_field(unit,linearized,timerange,params)
%LINEARIZED_PLACE_FIELD Estimate linearized place field.
%
%   RATEMAP = LINEARIZED_PLACE_FIELD(UNIT,LINEARIZED,TIMERANGE,DIRECTION,PARAMS)
%   produces linearized place field estimate(s), given a single-unit spike data
%   scalar struct UNIT and matching linearized position data struct LINEARIZED,
%   TIMERANGE, running DIRECTION and PARAMS.
%
%   UNIT must be a scalar struct of single-unit spike data which validates with
%   IS_UNIT.
%
%   LINEARIZED must be a scalar struct of linearized position data which
%   validates with IS_LINEARIZED, whose subject, day, epoch, and environment
%   fields match those of UNIT.
%
%   TIMERANGE must be a cell vector (or 1-element cell array), where each cell
%   corresponds to a single trial. Each cell must contain a Nx2 matrix of uint32
%   timestamp intervals (of the type which validates with IS_TIMERANGE), and
%   these intervals must fall entirely within the time ranges of UNIT and
%   LINEARIZED, and intervals corresponding to different trials must not
%   overlap.
%
%   DIRECTION must assume the value either +1 or -1, indicating the direction of
%   running. This value is irrelevant to the computation, but it is included in
%   the output as meta-data. It is the caller's responsiblity to set this value
%   correctly.
%
%   PARAMS is a scalar struct with the following fields:
%     'discretization_timestep': must be a uint32 scalar, which specifies the
%         time discretization resolution. Smaller values yield better estimates
%         at greater computational cost.
%     'grid': a vector of evenly-spaced grid points along the linearized
%         position dimension at which the firing rate is to be estimated,
%         expressed in the same units as LINEARIZED.position.
%     'smoothing_width' is the standard deviation of the Gaussian smoothing
%         kernel, expressed in the same units as LINEARIZED.position. Enter
%         empty matrix [] for no smoothing.
%
%   The output RATEMAP is a scalar struct array. RATEMAP has the following
%   fields:
%     'uid': inherited from UNIT
%     'region': inherited from UNIT
%     'hemisphere': inherited from UNIT
%     'subject': inherited from UNIT and LINEARIZED
%     'tetrode': inherited from UNIT
%     'depth': inherited from UNIT
%     'reference': inherited from UNIT
%     'clustnum': inherited from UNIT
%     'cluster_quality': inherited from UNIT
%     'day': inherited from UNIT and LINEARIZED
%     'epoch': inherited from UNIT and LINEARIZED
%     'environment': inherited from UNIT and LINEARIZED
%     'timerange': inherited from TIMERANGE
%     'discretization_timestep': inherited from PARAMS
%     'grid': inherited from PARAMS
%     'smoothing_width': inherited from PARAMS
%     'kernel': the string 'gaussian'
%     'position_units': inherited from LINEARIZED
%     'direction': inherited from DIRECTION
%     'spikes': a struct with fields 'timestamp', 'x', 'trial', with the time,
%         linearized position and trial number of each spike
%     'firing_rate': estimated firing rate vector; each element corresponds to a
%         linearized position bin
%     'occupancy': occupancy vector; each element corresopnds to a linearized
%         position bin
%
%Depends on:
%   IS_UNIT (written by smk)
%   IS_LINEARIZED (written by smk)
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
if (exist('is_linearized') ~= 2)
  error(['LINEARIZED_PLACE_FIELD depends on m-file IS_LINEARIZED ' ...
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
if ~is_linearized(linearized) || ...
    ~(isscalar(linearized) || isvector(linearized))
  error(['LINEARIZED must be a scalar or vector struct of linearized ' ...
      'position data']);
end
if ~all(struct_cmp(linearized(1),linearized, ...
    {'subject','units','environment'}))
  error(['Elements of LINEARIZED must correspond to the same subject and ' ...
      'environment and must have the same units']);
end
if ~isequal(size(unit),size(linearized))
  error('UNIT and LINEARIZED must match in size');
end
if ~all(arrayfun(@(u,l) struct_cmp(u,l, ...
    {'subject','day','epoch','environment'}),unit,linearized))
  error(['Corresponding elements of UNIT and LINEARIZED must have ' ...
      'matching subject, day, epoch, environment metadata']);
end

if ~iscell(timerange) || ~isvector(timerange) || ...
    ~all(cellfun(@is_timerange,timerange))
  error(['TIMERANGE must be a cell vector whose cells contain ' ...
      'Nx2 matrix of uint32 timestamp intervals.']);
end
% Check that the contents of timerange are non-overlapping
t_concat = cell2mat(timerange);
if ~is_timerange(t_concat)
  error(['Time intervals corresponding to individual trials/passes ' ...
      'must not be overlapping']);
end
if ~isempty(diff_intervals(t_concat, ...
    [linearized.timestamp(1), linearized.timestamp(end)])) || ...
    ~isempty(diff_intervals(t_concat,unit.timerange))
  error(['All time intervals in TIMERANGE must fall within time ranges ' ...
      'of UNIT and LINEARIZED']);
end

if ~isscalar(direction) || ~isreal(direction) || ...
    ~any(direction == [+1, -1])
  error('DIRECTION must equal +1 (rightbound) or -1 (leftbound)');
end

if ~isstruct(params) || ~isscalar(params)
  error('PARAMS must be a scalar struct array');
end
if ~all(isfield(params, ...
    {'discretization_timestep','grid','smoothing_width'}))
  error('PARAMS is missing one or more required fields');
end
if ~isscalar(params.discretization_timestep) || ...
    ~isa(params.discretization_timestep,'uint32') || ...
    ~(params.discretization_timestep > 0)
  error(['''discretization_timestep'' field of PARAMS must be non-zero ' ...
      'uint32 scalar']);
end
if ~isnumeric(params.grid) || ~isvector(params.grid) || ...
    (numel(params.grid) < 3)
  error(['''grid'' field of PARAMS must be a numeric grid vector of ' ...
      'least 3 elements']);
end
if any(diff(params.grid) <= 0) || ...
    any(abs(diff(diff(params.grid))) > sqrt(max(eps(params.grid))))
  error(['grid spacings in PARAMS must be uniform monotonic increasing']);
end
if ~(isscalar(params.smoothing_sigma) && isreal(params.smoothing_sigma) && ...
    isfinite(params.smoothing_sigma) && (params.smoothing_sigma > 0)) && ...
    ~isequal(params.smoothing_sigma,[])
  error(['''smoothing_sigma'' field of PARAMS must be either ' ...
      'a real positive finite scalar or the empty matrix []']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DATA PREPARATION
num_trials = 0;
binned_data = struct( ...
    'x'      , [], ...
    'counts' , [], ...
    'trial'  , [] );
% Get time, position, phase, trial number for each spike (useful for generating
% scatterplots)
spikes = struct( ...
    'timestamp' , [], ...
    'x'         , [], ...
    'trial'     , [] );
bin_edges = (unit.timerange(1):params.discretization_timestep: ...
    unit.timerange(end))';
bin_centers = 0.5*(bin_edges(1:end-1) + bin_edges(2:end));
binned_data.counts = [ binned_data.counts; ...
    histc(double(unit.timestamp),bin_edges) ];
% histc returns an extra element at the end which we don't want
binned_data.counts(end) = [];
% Resampled linearized position
binned_data.x = [ binned_data.x; interp1(double(linearized.timestamp), ...
    double(linearized.position),double(bin_centers)) ];
% Label trials (allocate zeros here, and populate nonzero elements in the
% loop below)
binned_data.trial = [ binned_data.trial; zeros(size(bin_centers)) ];
for j = 1:numel(timerange)
  % Trials are numbered cumulatively over *all* elements of UNIT
  num_trials = num_trials + 1;
  spiketimes = unit.timestamp(ismember_intervals(unit.timestamp, ...
      timerange{j}));
  spikes.timestamp = [ spikes.timestamp; spiketimes ];
  spikes.x = [ spikes.x; ...
      interp1(double(linearized.timestamp), ...
      double(linearized.position), ...
      double(spiketimes),'spline') ];
  spikes.trial = [ spikes.trial; num_trials*ones(size(spiketimes)) ];
  binned_data.trial(numel(binned_data.trial) - numel(bin_centers) + ...
      find(ismember_intervals(bin_centers,timerange{j}))) = ...
      num_trials;
end

% Exclude bins whose trial number is zero
include_idx = find(binned_data.trial > 0);
binned_data.x = binned_data.x(include_idx);
binned_data.counts = binned_data.counts(include_idx);
binned_data.trial = binned_data.trial(include_idx);
% Sanity checks
assert(all(isfinite(binned_data.x)));
assert(all(isfinite(binned_data.counts) & (binned_data.counts >= 0)));
assert(all(isfinite(binned_data.trial) & (binned_data.trial > 0) & ...
    (round(binned_data.trial) == binned_data.trial)));

% Combine data across trials and compute occupancy and place fields
Fs = TS_PER_SEC/double(params.discretization_timestep);
density_map_params = struct( ...
    'type',{'linear'}, ...
    'grid',{params.grid}, ...
    'dispersion',{params.smoothing_sigma});
[occupancy, firing_rate] = ...
    density_map(binned_data.x,binned_data.counts,Fs,density_map_params);

% Add metadata
ratemap.uid = unit.uid;
ratemap.region = unit.region;
ratemap.hemisphere = unit.hemisphere;
ratemap.subject = unit.subject;
ratemap.tetrode = unit.tetrode;
ratemap.depth = unit.depth;
ratemap.reference = unit.reference;
ratemap.clustnum = unit.clustnum;
ratemap.cluster_quality = unit.cluster_quality;
ratemap.day = unit.day;
ratemap.epoch = unit.epoch;
ratemap.environment = unit.environment;
ratemap.direction = direction;
ratemap.timerange = timerange;
ratemap.discretization_timestep = params.discretization_timestep;
ratemap.grid = params.grid;
ratemap.kernel = 'gaussian';
ratemap.smoothing_width = params.smoothing_width;
ratemap.position_units = linearized.units;
ratemap.spikes = spikes;
ratemap.firing_rate = firing_rate;
ratemap.occupancy = occupancy;


