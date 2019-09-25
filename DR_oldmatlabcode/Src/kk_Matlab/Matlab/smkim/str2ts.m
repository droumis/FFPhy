function timestamp = str2ts(timestring)
%STR2TS Convert a time string to an integer timestamp (1e4 units/second).
%   TIMESTAMP = STR2TS(STRING) returns the uint32 timestamp (expressed in units
%   of 0.1 milliseconds) that corresponds to STRING. STRING must be of the form
%   'HHH:MM:SS.XXXX'. The trailing sub-second characters "XXXX" may be given to
%   any precision up to 4 decimal places or omitted altogether, and the leading
%   hour characters "HHH" must include 1-3 digits. Examples of valid
%   timestrings:
%
%       '00:00:00.0000' (smallest allowed uint32 timestamp)
%       '0:20:14'
%       '00:20:14.'
%       '000:20:14.0000'
%       '119:18:16.7295' (largest allowed uint32 timestamp)
%
%   See also TS2STR, IS_TIMESTRING written by smk.
%
%Written by smk, 5 February 2009.
%

if ~ischar(timestring) || (size(timestring,1) ~= 1)
  error('timestring must be a single-row char argument');
end

% DO NOT TAMPER WITH THIS CONSTANT CONVERSION FACTOR!
TS_PER_SEC = 1e4;
% DO NOT TAMPER WITH THIS REGEX! 
re = '(?<h>^\d{1,3}):(?<m>\d{2}):(?<s>\d{2})(?<r>\.\d{0,4}$|$)';
m = regexp(timestring,re,'names');
if isempty(m)
  error('%s is not a well-formed timestring',timestring);
else
  hour = str2double(m.h);
  minute = str2double(m.m);
  second = str2double(m.s);
  if isempty(m.r) | strcmp(m.r,'.')
    subsecond = 0;
  else
    subsecond = str2double(m.r);
  end
  if (hour > 119) || (hour < 0)
    error('Invalid timestring: hours must be between 0 and 99');
  end
  if (minute > 59) || (minute < 0)
    error('Invalid timestring: minutes must be between 0 and 59');
  end
  if (second > 59) || (second < 0)
    error('Invalid timestring: seconds must be between 0 and 59');
  end
  if (subsecond > 999) || (subsecond < 0)
    error('This is an unfortunate bug');
  end
end
double_value = round(TS_PER_SEC*(3600*hour + 60*minute + second + subsecond));
timestamp = uint32(double_value);
if (double(timestamp) ~= double_value)
  error('Invalid timestring: exceeds range of uint32');
end

