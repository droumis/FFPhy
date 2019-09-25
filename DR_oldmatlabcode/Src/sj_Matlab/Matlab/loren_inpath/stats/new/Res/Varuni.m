% VARUNI: Calculate variance of the interval [p,q] on the uniform distribution [0,1].
%

function s2 = varuni(p,q)
  nint = 1000;
  s2 = var([p:((q-p)/nint):q]);
  return;
