% PadCols: Given two matrices, pads the narrower one with columns of NaN's
%            (or other value).
%
%     Syntax:  [A,B] = padcols(A,B,{val})
%
%         A,B = matrices to be compared. Either may be of any size.
%         val = optional value to be padded [default = NaN].
%         ---------------------------------------------------------
%         A,B = adjusted matrices.
%

% RE Strauss, 6/8/93
%   3/9/03 -  make padded value optional; default to NaN rather than zero.

function [A,B] = padcols(A,B)
  if (~nargin) help padcols; return; end;
  
  if (nargin < 3) val = []; end;
  if (isempty(val)) val = NaN; end;

   [ra,ca] = size(A);
   [rb,cb] = size(B);

   if (ca < cb)            % A narrower than B
      A = [A, val*ones(ra,(cb-ca))];
   end;

   if (ca > cb)            % B narrower than A
      B = [B, val*ones(rb,(ca-cb))];
   end;
   return;
