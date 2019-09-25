function rawpos = read_nspike_fixpos(posfilename,mpegfilename,mpegoffsetfilename,timestampsfilename,session)
%READ_NSPIKE_FIXPOS Import raw position marker records outputted by nspike_fixpos.
%   RAWPOS = READ_NSPIKE_FIXPOS(POS_FILENAME,MPEG_FILENAME,MPEGOFFSET_FILENAME,
%   TIMESTAMPS_FILENAME,SESSION) reads video-frame timestamps and raw pixel
%   coordinates from the packed binary data format in POSFILENAME and returns a
%   RAWPOS struct array.
%
%   POSFILENAME must specify a data file of the type outputted by nspike_fixpos.
%   MPEGFILENAME, MPEGINDEXFILENAME and MPEGTIMESTAMPFILENAME must specify
%   appropriate arguments for nspike_fixpos.
%
%   SESSION must be a struct array with the following fields:
%         subject: a string
%             day: an integer
%           epoch: a string, descriptor of the epoch
%     environment: a string, descriptor of the environment
%       timerange: a 1x2 vector of uint32 timestamps
%          tstart: a time string of the form HH:MM:SS[.xxxx]
%            tend: a time string of the form HH:MM:SS[.xxxx]
%
%   RAWPOS is a struct array with the following fields:
%     subject: inherited from SESSION
%     day: inherited from SESSION
%     session: inherited from SESSION
%     environment: inherited from SESSION
%     timestamp: Nx1 vector of uint32 timestamps
%     xfront: Nx1 vector of single float coordinates
%     yfront: Nx1 vector of single float coordinates
%     xback: Nx1 vector of single float coordinates -- may be empty
%     yback: Nx1 vector of single float coordinates -- may be empty
%     units: the string 'pixels'
%     sources: cell array containing the strings POS_FILENAME, MPEG_FILENAME,
%     MPEGOFFSET_FILENAME, TIMESTAMPS_FILENAME, SESSION.source
%
%   See also IS_RAWPOS, written by smk.
%
%Depends on:
%   IS_SESSION (written by smk)
%   READ_BINARY_RECORDS (written by smk)
%   IS_RAWPOS (written by smk)
%   VIDEO_READER (MATLAB class written by smk)
%
% Written by SMK, 2009 February 09.
%

% How much fluctuation of timestamps is acceptable? This is expressed as
% fraction of the median timestamp difference between consecutive frames
TIMESTAMP_TOLERANCE = 0.1; % We tolerate 10% fluctuation from median

if (exist('is_session') ~= 2)
  error('READ_NSPIKE_FIXPOS depends on m-function IS_SESSION (written by smk)');
end
if (exist('read_binary_records') ~= 2)
  error(['READ_NSPIKE_FIXPOS depends on m-function READ_BINARY_RECORDS ' ...
    '(written by smk)']);
end
if (exist('video_reader') ~= 2)
  error('READ_NSPIKE_FIXPOS depends on VIDEO_READER class (written by smk)');
end
if (exist('is_rawpos') ~= 2)
  error('READ_NSPIKE_FIXPOS depends on m-fuction IS_RAWPOS (written by smk)');
end

% Verify that posfilename, mpegfilename, mpegoffsetfilename, timestampsfilename
% all specify real files and includes path
if ~ischar(posfilename)
  error('posfilename must be a string');
elseif (exist(posfilename) ~= 2)
  error('posfilename %s does not refer to a file on search path', ...
      posfilename);
elseif ~isdir(regexp(posfilename,'.*/+(?=[^/]+$)','match','once'))
  error('posfilename %s does not include path',posfilename);
end
if ~ischar(mpegfilename)
  error('mpegfilename must be a string');
elseif (exist(mpegfilename) ~= 2)
  error('mpegfilename %s does not refer to a file on search path', ...
      mpegfilename);
elseif ~isdir(regexp(mpegfilename,'.*/+(?=[^/]+$)','match','once'))
  error('mpegfilename %s does not include path',mpegfilename);
end
if ~ischar(mpegoffsetfilename)
  error('mpegoffsetfilename must be a string');
elseif (exist(mpegoffsetfilename) ~= 2)
  error('mpegoffsetfilename %s does not refer to a file on search path', ...
      mpegoffsetfilename);
elseif ~isdir(regexp(mpegoffsetfilename,'.*/+(?=[^/]+$)','match','once'))
  error('mpegoffsetfilename %s does not include path',mpegoffsetfilename);
end
if ~ischar(timestampsfilename)
  error('timestampsfilename must be a string');
elseif (exist(timestampsfilename) ~= 2)
  error('timestampsfilename %s does not refer to a file on search path', ...
      timestampsfilename);
elseif ~isdir(regexp(timestampsfilename,'.*/+(?=[^/]+$)','match','once'))
  error('timestampsfilename %s does not include path',timestampsfilename);
end

% Check whether mpegfilename, mpegoffsetfilename, timestampsfilename can be used
% to construct a video_reader object
try
  video_reader_obj = video_reader(mpegfilename, mpegoffsetfilename, ...
      timestampsfilename);
  delete(video_reader_obj);
catch
  error(['Could not initialize a VIDEO_READER object with arguments ' ...
      '%s %s %s'],mpegfilename,mpegoffsetfilename,timestampsfilename);
end

% Check whether SESSION is a valid session struct array
if ~is_session(session)
  error('SESSION does not appear to be a valid sessions data struct array');
end
% Check whether all elements of sessions have matching subject and day
if (numel(unique({session(:).subject})) > 1)
  error('SESSION argument has elements whose subject fields differ');
end
if (numel(unique([session(:).day])) > 1)
  error('SESSION argument has elements whose day fields differ');
end

% read data from posfilename using READ_BINARY_RECORDS
posformat = cell2struct({ ...
    'timestamp', 'uint32', 1; ...
    'xfront'   , 'int16' , 1; ...
    'yfront'   , 'int16' , 1; ...
    'xback'    , 'int16' , 1; ...
    'yback'    , 'int16' , 1 }, ...
    {'name','type','count'},2);
eoh = '%%ENDHEADER\n';
% Read the entire *.pos file
try
  pos = read_binary_records(posfilename,eoh,posformat,Inf);
catch 
  error('Could not read records from pos file %s',posfilename);
end
% Convert position fields from int16 to single values
pos.xfront = single(pos.xfront);
pos.yfront = single(pos.yfront);
pos.xback = single(pos.xback);
pos.yback = single(pos.yback);
if isempty(pos.timestamp)
  error('pos file %s contains zero records',posfilename);
end
% If the first timestamp is zero, delete the first element of each field (this
% is the result of a bug in nspike_fixpos)
if (pos.timestamp(1) == 0) 
  if (numel(pos.timestamp) > 1)
    pos.timestamp = pos.timestamp(2:end);
    pos.xfront = pos.xfront(2:end);
    pos.yfront = pos.yfront(2:end);
    pos.xback = pos.xback(2:end);
    pos.yback = pos.yback(2:end);
  else
    error(['first timestamp in %s is zero (due to a bug in nspike_fixpos) ' ...
      'but there is only one sample'],posfilename);
  end
end
% Check whether timestamps are out of order
if any(diff(pos.timestamp) <= 0)
  error('Timestamps in pos file %s are not monotonically increasing', ...
      posfilename);
end

% Read the entire timestamps file and check that the timestamps from the *.pos
% file comprise a slice of this vector of timestamps
timestampsformat = cell2struct({ ...
    'timestamp', 'uint32', 1}, ...
    {'name','type','count'},2);
try
  ts = read_binary_records(timestampsfilename,eoh,timestampsformat,Inf);
catch
  error('Could not read records from timestamps file %s',timestampsfilename);
end
if isempty(ts.timestamp)
  error('Timestamps file %s contains zero records',timestampsfilename);
end
if any(diff(ts.timestamp) <= 0)
  error('Timestamps in timestamps file %s are not monotonically increasing',...
      timestampsfilename);
end
if ~all(ismember(pos.timestamp,ts.timestamp)) || ...
    ~isequal(pos.timestamp,ts.timestamp( ...
    find(pos.timestamp(1)==ts.timestamp):...
    find(pos.timestamp(end)==ts.timestamp)))
  error(['Timestamps in pos file %s do not comprise a slice of the ' ...
       'timestamps file %s'],posfilename,timestampsfilename);
end

% Replace samples that have position values assigned to zero (sentinel value used
% by nspike_fixpos). Note that this should occur *after* the consistency check
% with the timestamps file
zero_idx = find((pos.xfront == 0) | (pos.yfront == 0) | ...
    (pos.xback == 0) | (pos.yback == 0));
if ~isempty(zero_idx)
  warning(['%d samples with zero position values were found. These ' ...
      'samples will be censored'],numel(zero_idx));
  pos.xfront(zero_idx) = NaN;
  pos.yfront(zero_idx) = NaN;
  pos.xback(zero_idx) = NaN;
  pos.yback(zero_idx) = NaN;
end

% Initialize output struct, inherited meta-data from SESSION. Each element of
% RAWPOS corresponds to an element of SESSION.
rawpos = rmfield(session,{'tstart','tend','timerange','source'});
for f = 1:length(posformat)
  name = posformat(f).name;
  [rawpos(:).(name)] = deal({});
end
[rawpos(:).units] = deal('pixels');
[rawpos(:).sources] = deal({mpegfilename mpegoffsetfilename ...
    timestampsfilename posfilename});
for i = 1:numel(session)
  % append session source to sources field cell array
  rawpos(i).sources = {rawpos(i).sources{:} session(i).source};

  % Select samples that span the session timerange (make sure to include
  % flanking samples whose timestamps are just outside of timerange)
  start_idx = find(pos.timestamp < session(i).timerange(1),1,'last');
  if isempty(start_idx)
    error(['session tstart (%d) is earlier than the first available ' ...
        'position sample (%d) in %s'],session(i).timerange(1), ...
        pos.timestamp(1),posfilename);
  end
  end_idx = find(pos.timestamp > session(i).timerange(end),1,'first');
  if isempty(end_idx)
    error(['session tend (%d) is later than the last available ' ...
        'position sample (%d) in %s'],session(i).timerange(end), ...
        pos.timestamp(end),posfilename);
  end
  % now populate the structure with the selected samples
  if isempty(start_idx:end_idx)
    warning(['no valid position samples within specified time range: ' ...
        '%s to %s. rawpos struct will contain empty data fields'], ...
        session(i).tstart,session(i).tend);
  end
  for j = 1:length(fieldnames(pos))
    name = subsref(fieldnames(pos),substruct('{}',{j}));
    rawpos(i).(name) = pos.(name)(start_idx:end_idx,:);
  end
  % if timestamps are unusual, report diagnostic information to the user
  timestamp_diffs = diff(double(rawpos(i).timestamp));
  median_ts_diff = median(timestamp_diffs);
  idx = find((timestamp_diffs < median_ts_diff*(1-TIMESTAMP_TOLERANCE)) | ...
      (timestamp_diffs > median_ts_diff*(1+TIMESTAMP_TOLERANCE)));
  for j = 1:numel(idx)
    warning('anomalous timestamp difference between position samles: %d', ...
        timestamp_diffs(idx(j)));
  end
  % check for redundancy: are xfront and xback fields identical? are yfront and
  % yback fields identical?
  if isequal(rawpos(i).xback,rawpos(i).xfront) && ...
      isequal(rawpos(i).yback,rawpos(i).yfront)
    warning(['front and back tracking coordinates are identical: ' ...
        'xback and yback fields will be made empty to indicate this fact']);
    rawpos(i).xback = single([]);
    rawpos(i).yback = single([]);
  end
end

%{
if ~is_rawpos(rawpos)
  error('There is a bug in either READ_NSPIKE_FIXPOS or IS_RAWPOS');
end
%}

