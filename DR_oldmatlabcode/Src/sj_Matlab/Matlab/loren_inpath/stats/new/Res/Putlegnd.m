% PUTLEGND: Writes a series of string labels onto the current plot.
%
%     Usage: putlegnd(x,y,strmat,{fontsize})
%
%           x,y =      proportional positions of x,y coordinates, in relation 
%                        to current axis ranges.
%           strmat =   [n x m] string matrix, containing the strings to be printed 
%                        (by row).
%           fontsize = optional font size for text strings [default=12].
%

% RE Strauss, 7/1/98

function putlegnd(x,y,strmat,fontsize)
  if (nargin < 4) fontsize = []; end;

  if (isempty(fontsize))
    fontsize = 12;
  end;

  [r,c] = size(strmat);

  y_incr = fontsize/200;
  for i = 1:r
    puttext(x,y-(i-1)*y_incr,strmat(i,:),[],fontsize);
  end;
  
  return;
