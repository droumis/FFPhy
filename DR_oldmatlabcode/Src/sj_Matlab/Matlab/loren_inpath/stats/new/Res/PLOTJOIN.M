% PLOTJOIN: Given two sets of X,Y vectors specifying two matched
%           scatterplots, plots one set as O's and the other as X's and
%           joins comparable observations with line segments.
%
%     Usage: plotjoin(X1,Y1,X2,Y2,{nosymbols})
%
%           X1,Y1 =     vectors specifying first set of points.
%           X2,Y2 =     vectors specifying second set of points.
%           nosymbols = optional boolean flag indicating, if true, that 
%                         connecting lines are to be drawn without symbols.
%

% RE Strauss, 5/3/97
%   8/20/99 - buffer axis boundaries; changed plot colors for Matlab v5.
%   7/20/00 - added 'nosymbols' option.

function plotjoin(X1,Y1,X2,Y2,nosymbols)
  if (nargin < 5) nosymbols = []; end;

  if (isempty(nosymbols))
    nosymbols = 0;
  end;

  X1 = X1(:);
  X2 = X2(:);
  Y1 = Y1(:);
  Y2 = Y2(:);

  n1 = length(X1);
  n2 = length(Y2);
  m1 = length(X2);
  m2 = length(Y2);
  if (any([n1~=n2,n2~=m1,m1~=m2]))
    error('  PLOTJOIN: Four vectors must be of same length');
  end;

%  set(gca,'box','on');
  if (~nosymbols)
    plot(X1,Y1,'ok',X2,Y2,'xk');
  end;
  hold on;
  for i=1:n1
    plot([X1(i),X2(i)],[Y1(i),Y2(i)],'k');
  end;
  putbnd([X1;X2],[Y1;Y2]);
  hold off;

  return;

