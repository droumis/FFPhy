% PadRows: Given two matrices, pads the shorter one with rows of NaN's
%            (or other value).
%
%     Syntax:  [A,B] = padrows(A,B,{val})
%
%         A,B = matrices to be compared. Either may be of any size.
%         val = optional value to be padded [default = NaN].
%         ---------------------------------------------------------
%         A,B = adjusted matrices.
%

% RE Strauss, 6/8/93
%   2/28/03 - pad with NaN's rather than zeros.
%   3/9/03 -  make padded value optional.

function [A,B] = padrows(A,B,val)
  if (~nargin) help padrows; return; end;
  
  if (nargin < 3) val = []; end;
  if (isempty(val)) val = NaN; end;

  [ra,ca] = size(A);
  [rb,cb] = size(B);

  if (ra < rb)            % A shorter than B
    A = [A; val*ones((rb-ra),ca)];
  end;

  if (ra > rb)            % B shorter than A
    B = [B; val*ones((ra-rb),cb)];
  end;
  return;
