% UNIDSVAR: Calculate variance of the uniform interval [p,q].
%
%     Usage: s2 = unidsvar(p,q)
%

% RE Strauss, 10/8/96

function s2 = unidsvar(p,q)
  s2 = var(linspace(p,q,5000));
  return;
