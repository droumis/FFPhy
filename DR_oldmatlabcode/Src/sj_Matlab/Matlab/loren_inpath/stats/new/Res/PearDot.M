% PEARDOT:  Objective function for PearsonFitMoments.
%             c = [a c0 c1 c2]
%

% RE Strauss, 4/20/95
%   111/19/02 - changed name of calling function from pearcurv().

function dpdx = peardot(x,p,c)
  a = c(1);
  c0 = c(2);
  c1 = c(3);
  c2 = c(4);

  d = c0 + c1*x + c2*x*x;
  dpdx = sign(d) * (-p) * (a + x) / max(eps,abs(d));
  return;
