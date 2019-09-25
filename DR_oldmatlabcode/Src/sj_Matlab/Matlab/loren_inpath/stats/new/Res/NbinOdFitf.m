% NBINODFITF: Objective function for NbinOdFit().
%
%     Usage: dev = nbinodfitf(param,fx,stat)
%
%         param =   [m,v].
%         fx =      absolute frequencies (counts) corresponding to 
%                     X=0:(length(fx)-1).
%         stat =    value indicating the optimization statistic to be 
%                     minimized in the fitting:
%                       0 : sum of squared deviations [default].
%                       1 : chi-squared deviations.
%                       2 : sum of absolute deviations.
%         --------------------------------------------------------------
%         dev =     value of optimization statistic.
%

% Hilborn, R. & M. Mangel. 1997. The ecological detective: confronting
%   models with data.  Princeton Univ. Press.

% RE Strauss, 2/19/03
%   3/4/03 - added chi-squared and absolute deviations.

function dev = nbinodfitf(param,fx,stat)
  m = param(1);
  v = param(2);
  
  x = 0:length(fx)-1;
  y = nbinodcdf(x,m,v).*sum(fx);
  
  switch (stat)
    case 0,
      dev = sum((fx(:)-y(:)).^2);
    case 1,
      dev = sum((fx(:)-y(:)).^2./y(:));
    case 2,
      dev = sum(abs(fx(:)-y(:)));
  end;
  
  return;
  
