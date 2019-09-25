% Range: calculates range as max(x)-min(x), by column or overall.
%
%     Usage: r = range(x,{total})
%
%         x =     [n x p] matrix.
%         total = optional boolean vector indicating, if true, that a single range
%                   is to be calculated across the entire matrix rather than
%                   by column [default = 0].
%         ------------------------------------------------------------------------
%         r =     [1 x p] vector of column ranges (if total==0) or a scalar range
%                   if (total==1).
%

% RE Strauss, 1/8/03

function r = range(x,total)
  if (nargin < 1) help range; return; end;
  
  if (nargin < 2) total = []; end;
  
  if (isempty(total)) total = 0; end;
  
  if (total)
    x = x(:);
  end;
  
  r = max(x) - min(x);
  
  return;
  
