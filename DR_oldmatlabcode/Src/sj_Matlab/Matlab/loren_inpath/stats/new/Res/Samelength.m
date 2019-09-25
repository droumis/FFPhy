% SAMELENGTH: Determines whether a series of matrices have the same number of 
%             rows.  Optionally expands scalars to column vectors and row 
%             vectors to matrices having the common number of rows.
%
%     Usage: [ok,x1,...,x9] = samelength(x1,{x2},...,{x9})
%
%         x1,{x2},...,{x9} = up to 9 matrices
%         ----------------------------------------------------------------------
%         ok =        boolean value indicating whether input matrices have the 
%                       same number of rows.  If any corresponding output 
%                       matrices are also requested, expands scalars to column 
%                       vectors and row vectors to matrices before making 
%                       this determination; otherwise only the sizes of the 
%                       input matrices are considered.
%         x1,...,x9 = output vectors corresponding to input matrices, with 
%                       scalars and row vectors expanded.  If b is false, 
%                       original matrices are returned.
%

% RE Strauss, 1/3/01
%   2/13/01 - changed 'b' to 'ok'

function [ok,x1,x2,x3,x4,x5,x6,x7,x8,x9] = samelength(x1,x2,x3,x4,x5,x6,x7,x8,x9)
  len = zeros(1,nargin);

  if (nargin>=1) len(1) = size(x1,1); end;
  if (nargin>=2) len(2) = size(x2,1); end;
  if (nargin>=3) len(3) = size(x3,1); end;
  if (nargin>=4) len(4) = size(x4,1); end;
  if (nargin>=5) len(5) = size(x5,1); end;
  if (nargin>=6) len(6) = size(x6,1); end;
  if (nargin>=7) len(7) = size(x7,1); end;
  if (nargin>=8) len(8) = size(x8,1); end;
  if (nargin>=9) len(9) = size(x9,1); end;

  ulen = uniquef(len,1);
  nulen = length(ulen);

  ok = 0;
  if (nulen>2)
    return;
  end;
  if (nulen==1)
    ok = 1;
    return;
  end;
  if (ulen(1)~=1 | nargout==1)
    return;
  end;

  ok = 1;
  s = ones(max(len),1);

  if (nargin>=1), if (len(1)==1) x1 = s*x1; end; end;
  if (nargin>=2), if (len(2)==1) x2 = s*x2; end; end;
  if (nargin>=3), if (len(3)==1) x3 = s*x3; end; end;
  if (nargin>=4), if (len(4)==1) x4 = s*x4; end; end;
  if (nargin>=5), if (len(5)==1) x5 = s*x5; end; end;
  if (nargin>=6), if (len(6)==1) x6 = s*x6; end; end;
  if (nargin>=7), if (len(7)==1) x7 = s*x7; end; end;
  if (nargin>=8), if (len(8)==1) x8 = s*x8; end; end;
  if (nargin>=9), if (len(9)==1) x9 = s*x9; end; end;

  return;
