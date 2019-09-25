function recording_sites = read_recording_sites(filename);
%READ_RECORDING_SITES Read electrode recording sites from flat text into MATLAB data structures.
%   RECORDING_SITES = READ_RECORDING_SITES(FILENAME) reads electrode-location
%   meta-data from the flat text file specified by FILENAME (which must include
%   the full filesystem path and end with a *.TXT suffix) into a struct array
%   RECORDING_SITES. The flat text file must contain lines of the following
%   format:
%
%   [%s subject name] electrode[%d electrode number] day[%d day number]
%   [%s hemisphere] [%s region - this field may include spaces]
%
%   For example:
%
%   ratA  electrode09  day2  right  CA3
%
%   The subject, electrode, and hemisphere fields may only contain
%   [a-zA-Z0-9_]. The day field must be of the form day%1d. The region field
%   may contain any characters except newline. Everything after the whitespace
%   following the hemisphere field is interpreted as region.
%
%   Empty lines and lines that begin with the comment escape characer '#' are
%   ignored.
%
%   RECORDING_SITES is a struct array with the following fields:
%         subject: a string
%       electrode: integer index
%             day: an integer
%      hemisphere: a string
%          region: descriptor string
%          source: FILENAME
%
%   See also IS_RECORDING_SITES, written by smk.
%
%Depends on:
%   IS_RECORDING_SITES (written by smk)
%
%Written by smk, 2009 February 12.
%     

if (exist('is_recording_sites') ~= 2)
  error(['READ_RECORDING_SITES depends on m-file IS_RECORDING_SITES ' ...
      '(written by smk)']);
end

recording_sites = struct( ...
    'subject'    , {}, ...
    'electrode'  , {}, ...
    'day'        , {}, ...
    'hemisphere' , {}, ...
    'region'     , {}, ...
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
  error(['file %s is too big to be a recording sites flat text file; ' ...
      'please provide correct file name'],filename);
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
      continue
    end
    % Find stuff between whitespace
    w = regexp(l,'\S+','match');
    if (length(w) < 5)
      error('line does not contain correct number of tokens: %s',l);
    end
    % Parse subject string
    if ~isempty(regexp(w{1},'[^a-zA-Z0-9]'))
      error('string %s is not a valid subject name (line: %s)',w{1},l);
    else
      subject = w{1};
    end
    % Parse electrode string
    tmp = str2double(regexp(w{2},'(?<=^electrode)\d+(?=$)','match','once'));
    % note that str2double returns NaN for invalid strings
    if (round(tmp) ~= tmp) || (tmp < 0)
      error('string %s is not a valid electrode specifier (line: %s)',w{2},l);
    else
      electrode = tmp;
    end
    % Parse day string
    tmp = str2double(regexp(w{3},'(?<=^day)\d+(?=$)','match','once'));
    % note that str2double returns NaN for invalid strings
    if (round(tmp) ~= tmp) || (tmp < 0)
      error('string %s is not a valid day specifier (line: %s)',w{2},l);
    else
      day = tmp;
    end
    % Parse hemisphere string
    if ~isempty(regexp(w{4},'[^_a-zA-Z0-9]'))
      error('string %s is not a valid hemisphere descriptor (line: %s)',w{4},l);
    else
      hemisphere = w{4};
    end
    % Parse region string
    if isempty(strcat(w{5:end}))
      error('line %s does not contain a non-empty region descriptor',l);
    else
      region = '';
      for i = 1:length(w(5:end))
        region = [region, w{5+i-1}, ' '];
      end
      region = region(1:end-1);
    end
    % Sanity checks: For a given animal on a given day, electrode numbers must be unique
    same_subject_same_day = find(strcmp({recording_sites(:).subject},subject) & ...
        ([recording_sites(:).day] == day));
    if ~isempty(same_subject_same_day)
      if any([recording_sites(same_subject_same_day).electrode] == electrode)
        error('duplicate entries for same electrode');
      end
    end
    % if we made it all the way to the end without throwing any errors, then we
    % can safely append this info to the collection
    recording_sites(end+1,1) = struct( ...
        'subject'    , subject    , ...
        'electrode'  , electrode  , ...
        'day'        , day        , ...
        'hemisphere' , hemisphere , ...
        'region'     , region     , ...
        'source'     , filename   );
  end
end
fclose(f);
if isempty(recording_sites)
  error('no recording sites information found in %s',filename);
end

if ~is_recording_sites(recording_sites)
  error(['There is a bug in either READ_RECORDING_SITES or ' ...
      'IS_RECORDING_SITES']);
end

