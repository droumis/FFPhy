function dio = read_dio(diofilename,configfilename,session)
%READ_DIO Import digital I/O events from NSPike dio file.
%   DIO = READ_DIO(DIOFILENAME,CONFIGFILENAME,SESSION) reads
%   timestamps and digital I/O status codes from the packed binary format in
%   DIOFILENAME and returns a DIO struct. CONFIGFILENAME is read to determine
%   the configuration of the digital I/O ports
%
%   DIOFILENAME must specify an NSpike digital I/O event file of the type
%   outputted by nspike_extract.
%
%   CONFIGFILENAME must specify an NSpike config file of the type outputted by
%   nspike_extract, which contains configuration information for the main
%   NSpike host that controls digital I/O.
%
%   SESSION must be a struct array with the following fields:
%         subject: a string
%             day: an integer
%           epoch: a string, descriptor of the session
%     environment: a string, descriptor of the environment
%          tstart: a time string of the form HH:MM:SS[.xxxx]
%            tend: a time string of the form HH:MM:SS[.xxxx]
%
%   DIO is a struct array with the following fields:
%     subject: inherited from SESSION
%     day: inherited from SESSION
%     epoch: inherited from SESSION
%     environment: inherited from SESSION
%     timestamp: Nx1 vector of uint32 timestamps
%     status: Nxnum_ports vector of uint16 codes (num_ports is 4 in the
%       current version of NSpike).
%     ports: cell array of length num_ports containing strings, either 'input'
%       or 'output' to indicate the function of the port
%     source: DIOFILENAME
%
%   See also READ_DIO, written by smk.
%
%Depends on:
%   READ_BINARY_RECORDS (written by smk)
%   IS_SESSION (written by smk)
%   IS_DIO (written by smk)
%
%
%Written by SMK, 2009 April 28.
%

if (exist('is_session') ~= 2)
  error('READ_DIO depends on m-file IS_SESSION (written by smk)');
end
if (exist('read_binary_records') ~= 2)
  error('READ_DIO depends on m-file READ_BINARY_RECORDS (written by smk)');
end
if (exist('is_dio') ~= 2)
  error('READ_DIO depends on m-file IS_DIO (written by smk)');
end

% Verify that diofilename specifies a real file and includes path
if ~ischar(diofilename)
  error('diofilename must be a string');
elseif (exist(diofilename) ~= 2)
  error('diofilename %s does not refer to a valid file on search path',diofilename);
elseif ~isdir(fileparts(diofilename))
  error('diofilename %s does not include path',diofilename);
elseif isempty(regexp(diofilename,'^.+(?=\.dio$)','match','once'))
  error('file %s does not have *.dio suffix',diofilename);
end

% Read config file to determine the number and types of the digital I/O ports
if ~ischar(configfilename)
  error('configfilename must be a string');
elseif (exist(configfilename) ~= 2)
  error('configfilename %s does not refer to a valid file on search path',configfilename);
elseif ~isdir(fileparts(configfilename))
  error('configfilename %s does not include path',configfilename);
elseif isempty(regexp(configfilename,'^.+(?=\.config$)','match','once'))
  error('file %s does not have *.config suffix',configfilename);
end
d = dir(configfilename);
if isempty(d)
  error('no files found that match file specifier %s',configfilename);
end
if (length(d) > 1)
  error('more than one file matches file specifier %s; be more specific', ...
      configfilename);
end
% guard against memory crash
FILE_TOO_BIG = 1e5; 
if (d.bytes > FILE_TOO_BIG)
  error(['file %s is too big to be a config file; please provide ' ...
      'correct config file names'],configfilename);
end
fid = fopen(configfilename,'r');
if (fid == -1)
  error('Could not open file %s',configfilename);
end
while isempty(regexp(fgets(fid),'%BEGINCONFIG\n'))
  if feof(fid)
    error('File %s does not have %BEGINCONFIG character sequence', ...
        configfilename);
  end
end
txt = fread(fid,Inf,'*char')';
datatype_matches = regexp(txt,'datatype[[ \t]+\w+]*[ \t]+DIGITALIO');
if isempty(datatype_matches)
  error('config file %s does not contain DIGITALIO datatype');
end
num_ports = str2double(regexp(txt, ...
    '(?<=\sndioports[ \t]+)\d+','match'))';
if numel(num_ports) ~= 1
  error('config file %s does not contain proper ndioports line');
end
for p = 1:num_ports
  port_type = regexp(txt, ...
      sprintf('(?<=\\sdioport[ \\t]+%d[ \\t]+)\\w+',p-1),'match','once');
  if strcmp(port_type,'input') || strcmp(port_type,'output')
    ports{p} = port_type;
  else
    error('port type %s is not recognized',port_type);
  end
end

% Check that SESSION is a valid session struct array
if ~is_session(session)
  error('SESSION does not appear to be a valid sessions data struct array');
end
% Check whether all elements of session have matching subject and day
if (numel(unique({session(:).subject})) > 1)
  error('SESSION argument has elements whose subject fields differ');
end
if (numel(unique([session(:).day])) > 1)
  error('SESSION argument has elements whose day fields differ');
end

% read data from *.dio file using READ_BINARY_RECORDS
recformat = cell2struct({ ...
    'timestamp', 'uint32', 1; ...
    'status'   , 'uint16', num_ports }, ...
    {'name','type','count'}, 2);
eoh = ''; % no header
try
  all_records = read_binary_records(diofilename,eoh,recformat,Inf);
catch 
  error('Could not read records from file %s',diofilename);
end
if isempty(all_records.timestamp)
  error('File %s contains zero records',diofilename);
end
% check whether timestamps are out of order
if any(diff(all_records.timestamp) < 0)
  error('Timestamps in %s are not monotonically increasing', ...
      diofilename);
end

% Inherit meta-data from SESSION
dio = rmfield(session,{'tstart','tend'});
% Note which ports are input/output
[dio(:).ports] = deal(ports);
% Keep record of the original .dio file
[dio(:).source] = deal(diofilename);

for i = 1:numel(dio)
  % select samples that span the time interval from TSTART to TEND
  start_idx = find(all_records.timestamp >= session(i).timerange(1),1, ...
      'first');
  end_idx = find(all_records.timestamp < session(i).timerange(2),1,'last');
  select_idx = (start_idx:end_idx)';
  if isempty(select_idx)
    warning('no dio events within session %s (%s to %s)', ...
        session(i).epoch,session(i).tstart,session(i).tend);
  end
  for j = 1:length(fieldnames(all_records))
    name = subsref(fieldnames(all_records),substruct('{}',{j}));
    dio(i).(name) = all_records.(name)(select_idx,:);
  end
  % Report diagnostic information to the user
  timestamp_diffs = diff(dio(i).timestamp);
  if any(timestamp_diffs < 0)
    error('timestamps are out of order');
  end
  % ppend to output 
  dio(i,1) = orderfields(dio(i),fieldnames(dio));
end

if ~is_dio(dio)
  error('There is a bug in either READ_DIO or IS_DIO');
end

