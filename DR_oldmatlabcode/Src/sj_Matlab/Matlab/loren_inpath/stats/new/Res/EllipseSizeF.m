% EllipseSizeF: Objective function for numerical integration in EllipseSize().
%
%     Usage: fk = ellipsesizef(theta,k)
%

% RE Strauss, 12/7/01

function fk = ellipsesizef(theta,k)
  s = sin(theta);
  fk = sqrt(1 - (k.*k).*(s.*s));
  
  return;