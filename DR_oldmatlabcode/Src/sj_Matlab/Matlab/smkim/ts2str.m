function timestring = ts2str(timestamp)
%TS2STR Convert an integer timestamp (1e4 units/second) to a time string.
%   STRING = TS2STR(TIMESTAMP) returns a string of the form "HH:MM:SS.xxxx" that
%   corresponds to the uint32 timestamp TIMESTAMP (expressed in units of 1e-4
%   seconds). TIMESTAMP can not be so large that "HH" exceeds 99.
%
%   See also STR2TS, IS_TIMESTRING written by smk.
%
%Written by smk, 5 February 2009.
%

% DO NOT TAMPER WITH THIS CONSTANT CONVERSION FACTOR!
TS_PER_SEC = 1e4;

if ~isa(timestamp,'uint32') || ~isreal(timestamp) || ~isscalar(timestamp)
  error('TIMESTAMP argument must be real uint32 scalar');
end

% perform calculations in double precision, because MATLAB behaves strangely
% when doing integer arithmetic (e.g. does not always truncate in the expected
% ways)
timestamp = double(timestamp);
xxxx = rem(timestamp,TS_PER_SEC);
ss = rem(floor(timestamp/TS_PER_SEC),60);
mm = rem(floor(timestamp/TS_PER_SEC/60),60);
hhh = floor(timestamp/TS_PER_SEC/3600);

timestring = sprintf('%03d:%02d:%02d.%04d',hhh,mm,ss,xxxx);

