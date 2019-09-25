% IsOdd:  Boolean function indicating whether the integers in a matrix are (=1)
%         or are not (=0) odd.  Non-integer values are rounded.
%
%     Usage: b = isodd(X)
%
%         X =  input matrix.
%         ---------------------------------
%         b = corresponding boolean matrix.
%

% RE Strauss, 5/31/02

function b = isodd(X)
  b = mod(round(X),2);
  return;
  