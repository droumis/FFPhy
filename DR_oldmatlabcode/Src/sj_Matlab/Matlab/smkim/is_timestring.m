function bool = is_timestring(timestring)
%IS_TIMESTRING Validate whether the input is a valid timestring.
%
%   IS_TIMESTRING(TIMESTRING) returns true if TIMESTRING is a string of the
%   format "HH:MM:SS[.XXXX]" which resolves to a valid uint32 timestamp.
%   
%   See also STR2TS, TS2STR written by smk.
%
%Depends on:
%   STR2TS (written by smk)
%
%Written by SMK 2009 June 22.
%

if (exist('str2ts') ~= 2)
  error('IS_TIMESTRING depends on m-file STR2TS');
end

try
  timestamp = str2ts(timestring);
  if (isa(timestamp,'uint32') && isreal(timestamp))
    bool = true;
  else
    error('something is wrong with STR2TS');  
  end
catch
  bool = false;
end

