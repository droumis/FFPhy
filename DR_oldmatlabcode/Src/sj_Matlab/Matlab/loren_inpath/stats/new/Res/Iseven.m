% IsEven: Boolean function indicating whether the integers in a matrix are (=1)
%         or are not (=0) even.  Non-integer values are rounded.
%
%     Usage: b = iseven(X)
%
%         X =  input matrix.
%         ---------------------------------
%         b = corresponding boolean matrix.
%

% RE Strauss, 5/31/02

function b = iseven(X)
  b = abs(mod(round(X),2)-1);
  return;
  