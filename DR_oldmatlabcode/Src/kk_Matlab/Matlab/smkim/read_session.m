function session = read_session(filename);
%READ_SESSION Read session data from flat text into MATLAB data structures.
%   SESSION = READ_SESSION(FILENAME) reads session meta-data from the
%   flat text file specified by FILENAME (which must include the full
%   filesystem path and end with a *.TXT suffix) into a struct array
%   SESSION. The flat text file must contain lines of the following format:
%
%   [%s subject name] day[%d day num] [%s epoch name] [%s environment]
%   [hh:mm:ss.xxxx start time string] [hh:mm:ss.xxxx end time string]
%
%   Allowed characters in the %s fields are [a-zA-Z0-9]. For example:
%
%   ratA  day2  run1  Wtrack  02:34:56  02:50:01
%
%   Empty lines and lines that begin with the comment escape characer '#' are
%   ignored.
%
%   SESSION is a struct array with the following fields:
%         subject: a string
%             day: an integer
%         session: a string, descriptor of the session
%     environment: a string, descriptor of the environment
%          tstart: a time string of the form HH:MM:SS[.xxxx]
%            tend: a time string of the form HH:MM:SS[.xxxx]
%       timerange: a 1x2 vector of uint32 timestamps, equal to 
%                  [str2ts(tstart) str2ts(tend)]
%          source: FILENAME
%
%   See also IS_SESSION, written by smk.
%
%Depends on:
%   STR2TS (written by smk)
%   IS_SESSION (written by smk)
%
%Written by smk, 2009 February 12.
%     

if (exist('str2ts') ~= 2)
  error('READ_SESSION depends on m-file STR2TS (written by smk)');
end
if (exist('is_session') ~= 2)
  error('READ_SESSION depends on m-file IS_SESSION (written by smk)');
end

session = struct( ...
    'subject'    , {}, ...
    'day'        , {}, ...
    'epoch'      , {}, ...
    'environment', {}, ...
    'tstart'     , {}, ...
    'tend'       , {}, ...
    'timerange'  , {}, ...
    'source'     , {});

% Verify that filename specifies a real file and includes path
if ~ischar(filename)
  error('filename must be a string');
elseif (exist(filename) ~= 2)
  error('filename %s does not refer to a valid file on search path',filename);
elseif ~isdir(fileparts(filename))
  error('filename %s does not include path',filename);
elseif isempty(regexp(filename,'^.+(?=\.txt$)','match','once'))
  error('file %s does not have *.txt suffix',filename);
end
d = dir(filename);
if isempty(d)
  error('no files found that match file specifier %s',filename);
end
if (length(d) > 1)
  error('more than one file matches file specifier %s; be more specific', ...
      filename);
end
% guard against memory crash
FILE_TOO_BIG = 1e5; 
if (d.bytes > FILE_TOO_BIG)
  error(['file %s is too big to be a sessions flat text file; please ' ...
      'correct file name'],filename);
end
f = fopen(filename,'r');
if (f == -1)
  error('could not open file %s',filename);
end
while 1
  try
    l = fgetl(f);
  catch
    error('could not read line from file %s: %s',filename,ferror(f));
  end
  if (l == -1)
    break
  else
    % Strip away any trailing comments
    l = regexp(l,'^[^#]*','match','once');
    if isempty(l)
      continue;
    end
    % Find stuff between whitespace
    w = regexp(l,'\S+','match');
    if (length(w) == 0)
      continue;
    end
    if (length(w) ~= 6)
      error('line does not contain correct number of tokens: %s',l);
    end
    % Parse subject string
    if ~isempty(regexp(w{1},'[^a-zA-Z0-9]'))
      error('string %s is not a valid subject name (line: %s)',w{1},l);
    else
      subject = w{1};
    end
    % Parse day string
    tmp = str2double(regexp(w{2},'(?<=^day)\d+(?=$)','match','once'));
    % note that str2double returns NaN for invalid strings
    if (round(tmp) ~= tmp) || (tmp < 0)
      error('string %s is not a valid day specifier (line: %s)',w{2},l);
    else
      day = tmp;
    end
    % Parse epoch string
    if ~isempty(regexp(w{3},'[^a-zA-Z0-9]'))
      error('string %s is not a valid epoch name (line: %s)',w{3},l);
    else
      epoch = w{3};
    end
    % Parse environment string
    if ~isempty(regexp(w{4},'[^a-zA-Z0-9]'))
      error('string %s is not a valid environment name (line: %s)',w{4},l);
    else
      environment = w{4};
    end
    % Parse tstart and tend strings
    try
      if str2ts(w{6}) <= str2ts(w{5})
        error('time strings tstart and tend are out of order: %s',l);
      else
        tstart = w{5};
        tend = w{6};
        timerange = cellfun(@str2ts,{tstart, tend});
      end
    catch
      error('could not parse time string in this line: %s',l);
    end
    % Sanity checks: For a given animal on a given day, each epoch
    % must be unique, and no timeranges may overlap.
    same_subject_same_day = find(strcmp({session(:).subject},subject) & ...
        ([session(:).day] == day));
    if ~isempty(same_subject_same_day)
      if any(strcmp({session(same_subject_same_day).epoch},epoch))
        error('duplicate epoch fields for the same subject on the same day');
      end
      other_timeranges = vertcat(session(same_subject_same_day).timerange);
      if any((timerange(1) >= other_timeranges(:,1)) & ...
          (timerange(1) <= other_timeranges(:,end))) || ...
          any((timerange(end) >= other_timeranges(:,1)) & ...
          (timerange(end) <= other_timeranges(:,end)))
        error('overlapping session times for the same subject on the same day');
      end
    end
    % if we made it all the way to the end without throwing any errors, then we
    % can safely append this entry onto the output struct array
    session(end+1,1) = struct( ...
        'subject'    , subject    , ...
        'day'        , day        , ...
        'epoch'      , epoch       , ...
        'environment', environment, ...
        'tstart'     , tstart     , ...
        'tend'       , tend       , ...
        'timerange'  , timerange  , ...
        'source'     , filename   );
  end
end
fclose(f);
if isempty(session)
  error('no session information found in %s',filename);
end

if ~is_session(session)
  error('there is bug in either READ_SESSION or IS_SESSION');
end


