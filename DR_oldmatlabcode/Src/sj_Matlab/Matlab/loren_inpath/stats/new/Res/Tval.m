% TVAL: 2-sample t-statistic value given a data vector x and a group-id vector g
%
%   Syntax: t = tval(x,g)
%

% RE Strauss, 5/20/96

function t = tval(x,g)
  grp_id = uniquef(g,1);
  if (length(grp_id)>2)
    error('  TVAL: more than 2 groups specified');
  end;

  x1 = x(g==grp_id(1));
  x2 = x(g==grp_id(2));

  n1 = length(x1);
  n2 = length(x2);

  if (n1<2 | n2<2)
    error('  TVAL: need at least 2 observations per group');
  end;

  m1 = mean(x1);
  m2 = mean(x2);

  v1 = var(x1);
  v2 = var(x2);

  t = (m2-m1)./sqrt((((n1-1).*v1+(n2-1).*v2)./(n1+n2-2)).*((n1+n2)./(n1.*n2)));

  return;
