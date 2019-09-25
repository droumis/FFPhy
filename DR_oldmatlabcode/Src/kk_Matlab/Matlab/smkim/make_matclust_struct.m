function clust = make_matclust_struct(filename,session)
%MAKE_MATCLUST_STRUCT Read waveforms file and create structs for Matclust.
%   CLUST = MAKE_MATCLUST_STRUCT(FILENAME, SESSION) reads a spike data struct
%   from FILENAME and returns a struct of structs that is used by Matclust.
%
%   load(FILENAME) must return a struct with the following fields:
%     timestamp: a Nx1 vector of uint32 (in units of 1e-4 seconds), where N is
%       the number of events
%     samples: WxCxN 3-dimensional array, where W equals samples per second and
%       C is the number of channels, expressed in microvolts. Note that these
%       values must not be sign-inverted, or else the results will be wrong.
%     thresholds: C-element row vector of thresholds for the channels,
%       expressed in microvolts. Note that these are absolute values (positive
%       values) that correspond to the amplitude of negative voltage
%       deflection.
%     Fs: sampling rate of voltage samples in each trigger-captured window
%   This struct can be produced from a NSpike *.tt file using READ_SPIKE.
%
%   SESSION must be a struct array with the following fields:
%         subject: a string
%             day: an integer
%           epoch: a string, name of the session
%     environment: a string, descriptor of the environment
%          tstart: a time string of the form HH:MM:SS[.xxxx]
%            tend: a time string of the form HH:MM:SS[.xxxx]
%   Error will be raised if the .subject and .day fields are not identical for
%   all elements of SESSION. Also, error will be raised if any of the
%   tstart, tend intervals overlap.
%
%   The output CLUST is a struct whose fields are structs: clustattrib,
%   clustdata, graphattrib. 
%
%Depends on:
%   IS_SESSION (written by smk)
%   IS_SPIKE (written by smk)
%   TS2STR (written by smk)
%   FASTBITSET (written by mpk)
%   MEASURE_SPIKE_SHAPE (written by smk)
%
%Written by smk, 2009 February 14.
%

if (exist('is_session') ~= 2)
  error('MAKE_MATCLUST_STRUCT depends on m-file IS_SESSION (written by smk)');
end
if (exist('is_spike') ~= 2)
  error('MAKE_MATCLUST_STRUCT depends on m-file IS_SPIKE (written by smk)');
end
if (exist('ts2str') ~= 2)
  error('MAKE_MATCLUST_STRUCT depends on m-file TS2STR (written by smk)');
end
if (exist('fastbitset') ~= 3)
  error(['MAKE_MATCLUST_STRUCT depends on mex-file FASTBITSET ' ...
      '(written by mpk)']);
end
if (exist('measure_spike_shape') ~= 2)
  error('MAKE_MATCLUST_STRUCT depends on m-file MEASURE_SPIKE_SHAPE (written by smk)');
end

% How finely do we interpolate spike waveforms to estimate shape features?
INTERP_FACTOR = 10; 

% Matclust-specific constants. Do not tamper without thinking.
MAX_NUM_CLUST = 100;
MAX_NUM_FILTERS = 32;
TS_PER_SEC = 1e4;

% Sessions will be divided into time filters of at least this duration
TIME_CHUNK = 60*5; % 5-minute intervals

% Verify that filename specifies a real file and includes path
if (exist(filename) ~= 2)
  error('filename %s does not specify valid file on search path', ...
    filename);
elseif ~isdir(fileparts(filename))
  error('filename %s does not include path',filename);
end

% Check that SESSION is a valid session data struct array
if ~is_session(session)
  error('SESSION does not appear to be a valid sessions data struct array');
end
% Check whether all elements of SESSION have matching subject and day
if (numel(unique({session(:).subject})) > 1)
  error('SESSION argument has elements whose subject fields differ');
end
if (numel(unique([session(:).day])) > 1)
  error('SESSION argument has elements whose day fields differ');
end

% Load waveforms from the specified file. Note that we pass a filename, instead
% of simply supplying the actual spike data struct, because we want to (1) store
% the filename in the matclust meta-data, and (2) verify that this file loads
% correctly.
try
  tmp = load(filename);
catch
  error('could not load file %s',filename)
end
% Check that spike is valid spike data struct scalar
if ~isfield(tmp,'spike') || ~is_spike(tmp.spike) || ...
      ~isscalar(tmp.spike)
  error('%s does not load a scalar spike data struct', filename);
end
spike = tmp.spike;
clear('tmp');
% Check whether spike matches subject/day
if (spike.day ~= session(1).day)
  error('spike struct day %d does not match session day %d', ...
      spike.day,session(1).day);
end
if ~strcmp(spike.subject,session(1).subject)
  error('spike struct subject %s does not match session subject %s', ...
      spike.subject,session(1).subject);
end

% Compute parameters 
try
  [params, paramnames] = calculate_params(spike);
  assert(size(params,2) == numel(paramnames));
catch
  error('Error in generating cluster params');
end

matclust_filename = [regexp(filename,'[^/]+(?=_spike.mat$)','match','once') ...
    '_matclust.mat'];
clust.clustattrib = struct( ...
    'clusters'              , {[]}                                    , ...
    'filterindex'           , {[]}                                    , ...
    'takenpolys'            , {[]}                                    , ...
    'eventeditindex'        , {[]}                                    , ...
    'clustersOn'            , {[]}                                    , ...
    'currentfilepath'       , {''}                                    , ...
    'currentfilename'       , {matclust_filename}                     , ...
    'currentparamfilename'  , {[]}                                    , ...
    'nodata'                , {0}                                     , ...
    'datafile'              , {filename}                              , ...
    'lastaction'            , {''}                                    , ...
    'newchanges'            , {0}                                     , ...
    'states'                , {[]}                                    , ...
    'currstate'             , {1}                                     , ...
    'pointexclude'          , {zeros([size(params,1) 1],'int32')}     , ...
    'pointinclude'          , {zeros([size(params,1) 1],'int32')}     , ...
    'cluster0attrib'        , {struct('color',{[1 1 1]},'show',{1})}  , ...
    'hiddenclusters'        , {zeros([MAX_NUM_CLUST 2])}              , ...
    'currclust'             , {1}                                     , ...
    'dependencies'          , {false([MAX_NUM_CLUST MAX_NUM_CLUST])}  );
% all points belong to the 'all' time filter, which is indexed by 1
% note that 'all' is a *superset* of the union of the time filters defined
% below, because it also includes the gaps which are not included in any of the
% time filters
timefilters = ones([size(params,1) 1],'int32');
timefilterranges = zeros(0,2);
timefilternames = cell([MAX_NUM_FILTERS 1]);
j = 1; % j is index into clustdata.timefilterranges and clustdata.timefilternames
timefilters = fastbitset(timefilters,j,logical(1));
timefilternames{j} = 'all';
timefilterranges(j,1) = -1 + min( ...
    params(1,1),double(min(arrayfun(@(s) s.timerange(1),session))));
timefilterranges(j,2) = +1 + max( ...
    params(end,1),double(max(arrayfun(@(s) s.timerange(2),session))));
% start populating timefilterranges from the second row onwards; at this point,
% we need j to equal 1
assert(j == 1)
for i = 1:length(session)
  % subdivide each session into contiguous timeranges that are approximately
  % TIME_CHUNK in duration
  subdivisions = uint32(round(linspace( ...
      double(session(i).timerange(1)),double(session(i).timerange(end)), ...
      1 + max(1, ...
      round(double(diff(session(i).timerange))/TS_PER_SEC/TIME_CHUNK)) )));
      % number of subdivision is obtained by rounding division of total session
      % length by TIME_CHUNK, with a minimum of 1 subdivision
  if (numel(subdivisions) > 2)
    subdivisions = [subdivisions(1:end-1)', subdivisions(2:end)'];
  end
  assert(is_timerange(subdivisions));
  assert(subdivisions(1,1) == session(i).timerange(1));    
  assert(subdivisions(end,2) == session(i).timerange(end));    
  for k = 1:size(subdivisions,1)
    % j is index into clustdata.timefilterranges and clustdata.timefilternames
    j = j + 1; 
    timefilterranges(j,:) = double(subdivisions(k,:));
    timefilternames{j,1} = sprintf('%s (%s - %s)',session(i).epoch, ...
        ts2str(subdivisions(k,1)),ts2str(subdivisions(k,2)));
    % for each time filter that we define, find the threshold-trigger events
    % that fall within that time filter and apply bitmask to the
    % corresponding element of clustdata.timefilters
    timefilters = fastbitset(timefilters,j, ...
        (spike.timestamp >= subdivisions(k,1)) & ...
        (spike.timestamp < subdivisions(k,2)));
  end
end
if size(timefilterranges,1) > MAX_NUM_FILTERS
  error('too many time filters are defined');
else  
  timefiltermemmap = zeros([MAX_NUM_FILTERS 1]);
  timefiltermemmap(1:size(timefilterranges,1)) = 1:size(timefilterranges,1);
end
timefiltersOn = zeros([1 MAX_NUM_FILTERS]);
timefiltersOn(1) = 1;
otherfilternames = cell([MAX_NUM_FILTERS 1]);
for i = 1:length(otherfilternames)
  if i > 1
    otherfilternames{i} = num2str(i);
  else
    otherfilternames{1} = '1 All points';
  end
end
otherfiltermemmap = zeros([MAX_NUM_FILTERS 1]);
otherfiltermemmap(1) = 1;
otherfiltersOn = zeros([1 MAX_NUM_FILTERS]);
otherfiltersOn(1) = 1;
datarange = zeros([2 size(params,2)]);
for i = 1:size(params,2)
  datarange(1,i) = min(params(:,i)) - 0.05*range(params(:,i));
  datarange(2,i) = max(params(:,i)) + 0.05*range(params(:,i));
end
customvar = struct('session',session);
clust.clustdata = struct( ...
    'filledparam'           , {ones([1 size(params,2)])}              , ...
    'params'                , {params}                                , ...
    'origparams'            , {size(params,2)}                        , ...
    'names'                 , {paramnames}                            , ...
    'timefilterranges'      , {timefilterranges}                      , ...
    'timefilters'           , {timefilters}                           , ...
    'timefiltermemmap'      , timefiltermemmap                        , ...
    'timefiltersOn'         , {timefiltersOn}                         , ...
    'otherfilters'          , {ones([size(params,1) 1],'int32')}      , ...
    'otherfiltermemmap'     , otherfiltermemmap                       , ...
    'otherfiltersOn'        , {otherfiltersOn}                        , ...
    'filtermemmap'          , {[timefiltermemmap; otherfiltermemmap]} , ...
    'filteredpoints'        , {true([size(params,1) 1])}              , ...
    'datarange'             , {datarange}                             , ...
    'UnitsPerSec'           , {TS_PER_SEC}                            , ... 
    'timefilternames'       , {timefilternames}                       , ...
    'otherfilternames'      , {otherfilternames}                      , ...
    'customvar'             , {customvar}                             );

viewbox = repmat([zeros([size(params,2) 1]) ones([size(params,2) 1])], ...
    [1 2 size(params,2)]);
clust.graphattrib = struct( ...
    'polydraw'              , {0}                                     , ...
    'drawsquare'            , {0}                                     , ...
    'magdrawsquare'         , {0}                                     , ...
    'magsquare'             , {[]}                                    , ...
    'handdown'              , {0}                                     , ...
    'handinfo'              , {[]}                                    , ...
    'polyg'                 , {[]}                                    , ...
    'currentpolyhighlight'  , {[]}                                    , ...
    'squarehighlighton'     , {0}                                     , ...
    'pointmoved'            , {0}                                     , ...
    'plothighlightclear'    , {1}                                     , ...
    'polypress'             , {0}                                     , ...
    'relativepolypress'     , {[]}                                    , ...
    'shiftclick'            , {0}                                     , ...
    'viewbox'               , {viewbox}                               , ...
    'oldviewbox'            , {viewbox}                               , ...
    'backgroundcolor'       , {[0 0 0]}                               , ...
    'resolutionfactor'      , {[1 1]}                                 , ...
    'nonhiddenpoints'       , {true([size(params,1) 1])}              , ...
    'blackedOutClusters'    , {[]}                                    );

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% utility function for estimating peak and trough
function [params, paramnames] = calculate_params(spike)
  MEMLIMIT = 1e8;
  INTERP_FACTOR = 10; % 10x spline interpolation
  num_channels = size(spike.samples,2);
  params = nan([size(spike.samples,3), 1+3*num_channels],'double');
  paramnames = {};
  try
    shape = measure_spike_shape(spike,INTERP_FACTOR,MEMLIMIT);
  catch
    error('error in MEASURE_SPIKE_SHAPE');
  end
  params(:,1) = double(spike.timestamp);
  paramnames{1,1} = 'timestamp';
  for c = 1:num_channels
    params(:,1 + c) = double((-1) * shape.trough_amplitude(:,c));
    paramnames{:,1 + c} = sprintf('trough amplitude %d',c);
    params(:,1 + c + num_channels) = double(shape.peak_amplitude(:,c));
    paramnames{:,1 + c + num_channels} = sprintf('peak amplitude %d',c);
    params(:,1 + c + 2*num_channels) = double(shape.peak_time(:,c) - ...
        shape.trough_time(:,c));
    paramnames{:,1 + c + 2*num_channels} = sprintf('trough-peak width %d',c);
  end
end

