% LINEJOIN: Given two sets of X,Y vectors specifying two matched
%           scatterplots, plots one set as O's and the other as X's and
%           joins comparable observations with line segments.
%
%     Usage: linejoin(X1,Y1,X2,Y2)
%
%           X1,Y1 - vectors specifying first set of points
%           X2,Y2 - vectors specifying second set of points
%

% RE Strauss, 3/5/95
%   8/20/99 - changed plot colors for Matlab v5.

function linejoin(X1,Y1,X2,Y2)
  n1 = length(X1);
  n2 = length(Y2);
  m1 = length(X2);
  m2 = length(Y2);
  if (any([n1~=n2,n2~=m1,m1~=m2]))
    err('  Four vectors must be of same length');
  end;

  hold on;
  set(gca,'box','on');
  plot(X1,Y1,'ok',X2,Y2,'xk');
  for i=1:n1
    plot([X1(i),X2(i)],[Y1(i),Y2(i)],'k');
  end;
  hold off;

  return;

