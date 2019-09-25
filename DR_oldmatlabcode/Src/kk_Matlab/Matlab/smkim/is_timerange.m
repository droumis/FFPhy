function bool = is_timerange(timerange)
%IS_TIMERANGE Validate whether the input is a valid uint32 timestamp interval.
%
%   IS_TIMERANGE(TIMERANGE) returns true if TIMERANGE is an Nx2 array of uint32
%   timestamp intervals.
%   
%   See also IS_TIMESTAMP, IS_INTERVALS.
%
%Depends on:
%   IS_INTERVAL (written by SMK)
%
%Written by SMK 2009 June 22.
%

if (exist('is_intervals') ~= 2)
  error('IS_TIMERANGE depends on m-file IS_INTERVALS (written by SMK)');
end

if is_intervals(timerange) && isa(timerange,'uint32')
  bool = true;
else
  bool = false;
end

