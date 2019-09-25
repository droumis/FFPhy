% TRUNCATE: Truncates (fixes) a matrix to a specified number of decimal places.
%
%     Usage: res = truncate(x,dp)
%
%           x =   matrix of real values.
%           dp =  number of decimal positions to be preserved.
%           --------------------------------------------------
%           res = corresponding 'de-fuzzed' matrix.
%

function res = truncate(x,dp)
  res = fix(x * 10.^dp) * 10.^(-dp);

  return;
