function unit = read_matclust(matclustfilename)
%READ_MATCLUST Extract sorted spikes during specified sessions from a matclust file.
%   UNIT = READ_MATCLUST(MATCLUSTFILENAME) returns a UNIT struct array that
%   contains spike times, waveforms and cluster quality measures for
%   each putative single unit.
%
%   MATCLUSTFILENAME must specify a standard Matclust file. This function
%   follows the file named in clustattrib.datafile to obtain spike waveforms;
%   an errors will be raised if this file is not available or if it does not
%   contain a valid spike data struct (of the format returned by READ_TT).
%   clustdata.customvar must have a 'session' field which is a struct array of
%   the format produced by READ_SESSION. For each cluster defined in the
%   Matclust structure, the time filters in which the cluster is defined are
%   compared against the session times in clustdata.customvar.session. For each
%   cluster, the UNIT struct array contains an element for each session
%   throughout which that cluster is defined.
%
%Depends on:
%   IS_SESSION (written by smk)
%   IS_SPIKE (written by smk)
%   MEASURE_SPIKE_SHAPE (written by smk)
%   MEASURE_CLUSTER_QUALITY (written by smk)
%   NANMEDIAN (MATLAB Statistics Toolbox)
%   
%Written by smk, 2009 February 14.
%

if (exist('is_session') ~= 2)
  error('READ_MATCLUST depends on m-file IS_SESSION (written by smk)');
end
if (exist('is_spike') ~= 2)
  error('READ_MATCLUST depends on m-file IS_SPIKE (written by smk)');
end
if (exist('measure_spike_shape') ~= 2)
  error(['READ_MATCLUST depends on m-file MEASURE_SPIKE_SHAPE ' ...
      '(written by smk)']);
end
if (exist('measure_cluster_quality') ~= 2)
  error(['READ_MATCLUST depends on m-file MEASURE_CLUSTER_QUALITY ' ...
      '(written by smk)']);
end
if (exist('nanmedian') ~= 2)
  error(['READ_MATCLUST depends on m-file NANMEDIAN ' ...
      '(in MATLAB Statistics Toolbox']);
end

MAX_NUM_TIMEFILTER = 32; % length of Matclust timefilternames cell array
TS_PER_SEC = 1e4;

% Verify that filename specifies a real file and includes path
if ~ischar(matclustfilename)
  error('filename must be a string');
elseif (exist(matclustfilename) ~= 2)
  error('filename %s does not refer to a valid file on search path', ...
      matclustfilename);
elseif ~isdir(regexp(matclustfilename,'.*/+(?=[^/]+$)','match','once'))
  error('filename %s does not include path',matclustfilename);
end
try
  clust = load(matclustfilename);
catch
  error('Could not load matclust file %s',matclustfilename);
end
try
  tmp = load(clust.clustattrib.datafile);
catch
  error('Could not load spike wavforms file %s',clust.clustattrib.datafile);
end

% Check that spike is valid spike waveforms struct array
if ~isfield(tmp,'spike') || ~is_spike(tmp.spike)
  error(['data loaded from %s does not appear to be a valid ' ...
      'spike waveforms data struct array'],clust.clustattrib.datafile);
end
spike = tmp.spike;
clear('tmp');
% Check whether all elements of sipke have identical subject/day
if (numel(unique({spike(:).subject})) > 1)
  error('spike data struct has elements whose subject fields differ');
end
if (numel(unique([spike(:).day])) > 1)
  error('spike data struct has elements whose day fields differ');
end

% Check that clust is valid
MATCLUST_REQUIRED_FIELDS = {'clustdata','clustattrib','graphattrib'};
if ~isstruct(clust) || ~all(isfield(clust,MATCLUST_REQUIRED_FIELDS))
  error('matclust data loaded from %s is missing required field', ...
      matclustfilename);
end
CLUSTDATA_REQUIRED_FIELDS = { ...
    'params'            , ...
    'names'             , ...
    'timefilterranges'  , ...
    'timefilternames'   , ...
    'timefiltermemmap'  , ...
    'UnitsPerSec'       };
if ~isstruct(clust.clustdata) || ~isscalar(clust.clustdata) || ...
    ~all(isfield(clust.clustdata,CLUSTDATA_REQUIRED_FIELDS))
  error('clustdata is missing required field');
end
if ~iscellstr(clust.clustdata.names) || ~(numel(clust.clustdata.names) > 1) ...
    (numel(unique(clust.clustdata.names)) ~= numel(clust.clustdata.names))
  error('clustdata.names is not a N+1 cell array of unique strings');
end
if ~isnumeric(clust.clustdata.params) || ~isreal(clust.clustdata.params) || ...
    ~isa(clust.clustdata.params,'double') || ...
    (ndims(clust.clustdata.params) ~= 2)
  error('clustdata.params is not a matrix of real doubles');
end    
if size(clust.clustdata.params,2) ~= numel(clust.clustdata.names)
  error('size mismatch between clustdata.params and clustdata.names');
end
if size(clust.clustdata.params,1) ~= size(spike.samples,3)
  error('size mismatch between clustdata.params and spike.samples');
end
if ~strcmp(clust.clustdata.names{1},'timestamp')
  error('first cell of clustdata.names must be "timestamp"');
end
if ~all(uint32(clust.clustdata.params(:,1)) == clust.clustdata.params(:,1))
  error('first column of clustdata.params must contain uint32 timestamps');
end
num_channels = size(spike.samples,2);
% clustdata.timefilternames must have at least two non-empty cells, and one of
% those must be the first cell (corresponding to the "all" time filter)
if ~iscell(clust.clustdata.timefilternames) || ...
    (numel(clust.clustdata.timefilternames) ~= MAX_NUM_TIMEFILTER) || ...
    (nnz(~cellfun(@isempty,clust.clustdata.timefilternames)) ~= ...
    numel(unique(clust.clustdata.timefilternames( ...
    ~cellfun(@isempty,clust.clustdata.timefilternames))))) || ...
    isempty(clust.clustdata.timefilternames{1}) || ...
    (nnz(~cellfun(@isempty,clust.clustdata.timefilternames)) <= 1)
  error('clustdata.timefilternames does not seem to be valid');
end
% clustdata.timefilterranges must have more than one row, and the first row
% must span the time intervals defined in the remaining rows. The remaining rows
% must form a valid timerange.
if ~isnumeric(clust.clustdata.timefilterranges) || ...  
    ~all(uint32(clust.clustdata.timefilterranges(:)) == ...
    clust.clustdata.timefilterranges(:))
    (ndims(clust.clustdata.timefilterranges) ~= 2) || ...
    ~(size(clust.clustdata.timefilterranges,1) > 1) || ...
    (size(clust.clustdata.timefilterranges,2) ~= 2) || ...
    ~all(clust.clustdata.timefilterranges(1,1) <= ...
    clust.clustdata.timefilterranges(:)) 
    ~all(clust.clustdata.timefilterranges(1,2) >= ...
    clust.clustdata.timefilterranges(:))
    ~is_timerange(uint32(clust.clustdata.timefilterranges(2:end,:))) || ...
  error('clustdata.timefilterranges does not seem to be valid');
end
if size(clust.clustdata.timefilterranges,1) ~= ...
    nnz(~cellfun(@isempty,clust.clustdata.timefilternames))
  error(['mismatch between clustdata.timefilterranges and number of ' ...
      'non-empty cells in clustdata.timefilternames']);
end
% ALERT: This code assuems that clustdata.timefiltermemmap indexes the time
% filters in the default trivial order [1; 2; 3; ...] with zero padding at the
% end; if any time filters were deleted or inserted, this condition will be
% violated, and this code fails
if ~isa(clust.clustdata.timefiltermemmap,'double') || ...
    ~isvector(clust.clustdata.timefiltermemmap) || ...
    (size(clust.clustdata.timefiltermemmap,2) ~= 1) || ...
    (numel(clust.clustdata.timefiltermemmap) ~= MAX_NUM_TIMEFILTER)
  error('clustdata.timefiltermemmap is not valid');
elseif (nnz(clust.clustdata.timefiltermemmap) ~= ...
    size(clust.clustdata.timefilterranges,1)) || ...
    ~all( clust.clustdata.timefiltermemmap( ...
    1:size(clust.clustdata.timefilterranges,1)) == ...
    (1:size(clust.clustdata.timefilterranges))' )
  error('clustdata.timefiltermemmap contains gaps or non-monotonic order');
end
% Check that clust.clustdata.customvar contains a valid session struct array
if ~isfield(clust.clustdata,'customvar') || ... 
    ~isstruct(clust.clustdata.customvar) ...
    ~isfield(clust.clustdata.customvar,'session') || ...
    ~is_session(clust.clustdata.customvar.session)
  error(['matclust struct does not contain a valid session data struct ' ...
      'array in its customvar field']);
else
  session = clust.clustdata.customvar.session;
end
% Check whether all elements of session have matching subject and day
if (numel(unique({session(:).subject})) > 1)
  error('session info struct array has elements whose subject fields differ');
end
if (numel(unique([session(:).day])) > 1)
  error('session info struct array has elements whose day fields differ');
end

% For each element of the session array, find matching time filters in the
% matclust clustdata struct.
session_timefilters = cell([numel(session), 1]);
for i = 1:numel(session)
  % Find all timefilterranges that fall within the session timerange
  j = find((clust.clustdata.timefilterranges(:,1) >= ...
    session(i).timerange(1)) & ...
    (clust.clustdata.timefilterranges(:,2) <= ...
    session(i).timerange(2)));
  % Remove the "all" filter from consideration (this can happen if the data
  % were acquired in a single recording session, in which case the "all"
  % interval is degenerate with a time filter)
  j(j == 1) = []; 
  % Sort the rows of timefilterranges that fall within the interval defined by
  % session(i).timerange
  t = sortrows(clust.clustdata.timefilterranges(j,:));
  % Check for agreement with the session time interval
  if ~all([t(1,1) t(end,end)] == session(i).timerange)
    error('matclust time filters do not agree with session time boundaries');
  end
  % If more than one timefilter falls within the interval, check that they
  % completely cover the interval and coincide at their boundaries
  if (size(t,1) > 1) && ~all(t(2:end,1) == t(1:end-1,2))
    error('matclust time filters do not agree with session time boundaries');
  end
  % Check the names of the matching time filters, and warn if they don't match
  % the name of the session
  if ~all(strcmp(session(i).epoch, ...
      regexp(clust.clustdata.timefilternames(j),'^\w+','match','once')))
    warning('clustdata.timefilternames do not match epoch name');
  end
  % Record these time filter indices for future lookup
  session_timefilters{i} = j;
end
% Sanity check: make sure that no two sessions correspond to the same time
% filters
for i1 = 1:numel(session)
  for i2 = i1:numel(session)
    if (i1 ~= i2) && ~isempty(intersect( ...
        session_timefilters{i1},session_timefilters{i2}))
      error('bug?! different sessions map onto the same time filters');
    end
  end
end
if clust.clustdata.UnitsPerSec ~= TS_PER_SEC
  error('clust.clustdata.UnitsPerSec is %f, does not match standard %f', ...
      clust.clustdata.UnitsPerSec,TS_PER_SEC);
end
CLUSTATTRIB_REQUIRED_FIELDS = {'clusters','filterindex'};
if ~isstruct(clust.clustattrib) || ...
    ~all(isfield(clust.clustattrib,CLUSTATTRIB_REQUIRED_FIELDS))
  error('clustattrib does not exist or is missing required clusters field');
end
% For each cluster, find matching time filters in the matclust clustdata struct.
cluster_timefilters = cell([numel(clust.clustattrib.clusters), 1]);
if isempty(clust.clustattrib.clusters)
  %pass 
elseif ~iscell(clust.clustattrib.clusters)
  error('clustattrib.clusters must be a cell array');
else
  for i = 1:numel(clust.clustattrib.clusters)
    if isempty(clust.clustattrib.clusters{i})
      %pass
    elseif ~isstruct(clust.clustattrib.clusters{i}) || ...
        ~isscalar(clust.clustattrib.clusters{i}) || ...
        ~all(isfield(clust.clustattrib.clusters{i}, ...
        {'defineaxes','index'})) || ...
        ~isvector(clust.clustattrib.clusters{i}.index) || ...
        ~isa(clust.clustattrib.clusters{i}.index,'uint32') || ...
        (numel(unique(clust.clustattrib.clusters{i}.index)) ~= ...
        numel(clust.clustattrib.clusters{i}.index)) || ...
        ~all(ismember(clust.clustattrib.clusters{i}.index, ...
        1:size(clust.clustdata.params,1))) || ...
        ~isnumeric(clust.clustattrib.clusters{i}.defineaxes) || ...
        ~isreal(clust.clustattrib.clusters{i}.defineaxes) || ...
        (ndims(clust.clustattrib.clusters{i}.defineaxes) ~= 2) || ...
        (size(clust.clustattrib.clusters{i}.defineaxes,2) ~= 3) || ...
        ~all(ismember(clust.clustattrib.clusters{i}.defineaxes(:,1), ...
        1:size(clust.clustdata.params,2))) || ...
        ~all(ismember(clust.clustattrib.clusters{i}.defineaxes(:,2), ...
        1:size(clust.clustdata.params,2))) || ...
        ~all(ismember(clust.clustattrib.clusters{i}.defineaxes(:,3), ...
        find(~cellfun(@isempty,clust.clustattrib.filterindex))))
      error('clustattrib.clusters{%d} is invalid',i);
    end
    % Iterate through defineaxes(:,3) to find time filters for this cluster;
    % this includes time filters in which an empty cluster polygon was defined
    % (i.e. unit did not fire any spikes)
    for j = 1:size(clust.clustattrib.clusters{i}.defineaxes,1)
      t = clust.clustattrib.filterindex{ ...
          clust.clustattrib.clusters{i}.defineaxes(j,3)}(1);
      if ~ismember(t,unique(clust.clustdata.timefiltermemmap( ...
          clust.clustdata.timefiltermemmap ~= 0)))
        error('clustattrib.clusters{%d}.defineaxes is invalid',i);
      else
        % Keep record of which time filters this cluster is defined over (i.e.,
        % in which time filters is a cluster box is drawn?)
        if isempty(cluster_timefilters{i})
          cluster_timefilters{i} = t;
        else
          cluster_timefilters{i} = unique([cluster_timefilters{i}; t]);
        end
      end
    end
  end
end

% Initialize empty UNIT struct
unit = struct( ...
    'uid'                       , {}, ...
    'sources'                   , {}, ...
    'subject'                   , {}, ...
    'tetrode'                   , {}, ...
    'depth'                     , {}, ...
    'region'                    , {}, ...
    'hemisphere'                , {}, ...
    'reference'                 , {}, ...
    'passbands'                 , {}, ...
    'thresholds'                , {}, ...
    'Fs'                        , {}, ...
    'samples_before_trigger'    , {}, ...
    'amplitude_cutoff'          , {}, ...
    'clustnum'                  , {}, ...
    'timestamp'                 , {}, ...
    'samples'                   , {}, ... 
    'day'                       , {}, ...
    'epoch'                     , {}, ...
    'environment'               , {}, ...
    'timerange'                 , {}, ...
    'overall_mean_rate'         , {}, ...
    'maximum_amplitude_channel' , {}, ...
    'shape'                     , {}, ...
    'cluster_quality'           , {});
% Steps to assembling a unit struct from a matclust file:
% 1. Inherit all data from spike waveforms
% 2. Inherit session meta-data and identify all threshold-trigger events in each
%    session
% 3. Replace 'samples' and 'timestamp' fields with data for threshold-trigger
%    events corresponding only to this unit, using cluster indices
%    (clust.clustattrib.clusters{clustnum}.index)
% 4. Compute waveform shape.
% 5. Compute cluster quality.

% Iterate over defined clusters
for i = 1:length(clust.clustattrib.clusters)
  if isempty(clust.clustattrib.clusters{i})
    continue;
  end
  % Iterate over sessions and determine whether cluster is defined in each
  % session
  for j = 1:numel(session)
    if ismember(1,cluster_timefilters{i})
      % cluster is defined for all times
    else
      overlap = intersect(session_timefilters{j},cluster_timefilters{i});
      if ~isempty(overlap)
        if isempty(setdiff(session_timefilters{j},overlap))
          % cluster is defined for all time filters that correspond to this
          % session
          disp(sprintf('%s cluster %d is defined for all of %s', ...
              matclustfilename,i,session(j).epoch));
        else
          % cluster is defined for some, but not all time filters that
          % correspond to this session (partial isolation); skip
          continue;
        end
      else
        continue;
      end
      % Construct unique identifier string for this putative single unit
      tmp_unit.uid = sprintf('%s_day%d_tetrode%02d_cluster%02d', ...
          spike.subject,spike.day,spike.tetrode,i);
      % Grab cluster number
      tmp_unit.clustnum = i;
      % Inherit data from spike waveforms struct
      spike_fieldnames = fieldnames(spike);
      for f = 1:numel(spike_fieldnames)
        tmp_unit.(spike_fieldnames{f}) = spike.(spike_fieldnames{f});
      end
      % Inherit data from session struct
      metadata_fields = {'subject','day','epoch','environment','timerange'};
      for f = 1:numel(metadata_fields)
        tmp_unit.(metadata_fields{f}) = session(j).(metadata_fields{f});
      end
      % sources field combines the name of the matclust file, the sources for
      % the spike file and the source for the session struct
      tmp_unit.sources = { ...
          matclustfilename, spike.sources{:}, session(j).source};
      % Indices of all threshold-crossing events that occur during this session
      event_idx = find( ...
          (uint32(clust.clustdata.params(:,1)) >= session(j).timerange(1)) & ...
          (uint32(clust.clustdata.params(:,1)) < session(j).timerange(end)));
      % Indices of the threshold-crossing events that occur during this session
      % *and* that belong to this single unit
      event_unit_idx = intersect(clust.clustattrib.clusters{i}.index,event_idx);
      % Modify timestamp and samples fields to contain only data for this unit
      tmp_unit.timestamp = tmp_unit.timestamp(event_unit_idx);
      tmp_unit.samples = tmp_unit.samples(:,:,event_unit_idx);
      % overall_mean_rate field is the number of spikes dividied by epoch duration
      tmp_unit.overall_mean_rate = numel(tmp_unit.timestamp) / ...
          duration_timerange(tmp_unit.timerange);
      % Measure spike shape
      tmp_unit.maximum_amplitude_channel = NaN; % we'll fill in this value later
      INTERPFACTOR = 10;
      MEMLIMIT = 1e8;
      tmp_unit.shape = measure_spike_shape(tmp_unit,10,MEMLIMIT);
      shape_fieldnames = fieldnames(tmp_unit.shape);
      for f = 1:numel(shape_fieldnames)
        % Remove redundant fields that are defined in the parent struct
        if isfield(tmp_unit,shape_fieldnames{f})
          tmp_unit.shape = rmfield(tmp_unit.shape,shape_fieldnames{f});
        end
      end
      % Measure cluster quality (be careful with indexing events within only
      % this session)
      FEATURE_IDX = 2:5;
      clust.clustdata.names{FEATURE_IDX};
      tmp_unit.cluster_quality = measure_cluster_quality( ...
          clust.clustdata.params(:,FEATURE_IDX), ...
          clust.clustdata.names(FEATURE_IDX),event_unit_idx, ...
          setdiff(event_idx,event_unit_idx));
    end
    % Append to the UNIT struct array
    unit(end+1,1) = orderfields(tmp_unit,unit);
  end
end

% now go back and process each unique unit across all epochs in which it was isolated
unique_uids = unique({unit(:).uid});
for i = 1:numel(unique_uids)
  j = find(strcmp(unique_uids{i},{unit(:).uid}));
  tmp_shape = [unit(j).shape];
  % Determine which channel has maximum trough amplitude
  [junk, maximum_amplitude_channel] = max( ...
          nanmedian(-vertcat(tmp_shape(:).trough_amplitude),1));  
  [unit(j).maximum_amplitude_channel] = deal(maximum_amplitude_channel);
  % Iterate through the (sub)fields of the shape field and keep only the column
  % that corresponds to maximum_amplitude_channel
  for k = 1:numel(j)
    shape_fields = fieldnames(unit(j(k)).shape);
    for l = 1:numel(shape_fields)
      unit(j(k)).shape.(shape_fields{l}) = ...  
          unit(j(k)).shape.(shape_fields{l})(:,maximum_amplitude_channel);
  end 
end

