% KSProb: Calculates Q_KS(lambda) = pr(D > Dmax), the significance function for the 
%         Kolmogorov-Smirnov test, as given by Press et al. (1992:624).
%
%     Usage: pr = ksprob(lambda)
%
%         lambda =  expression that is a function of D, sample size, and (for the 2D KS test)
%                     the correlation coefficient.
%         -----------------------------------------------------------------------------------
%         pr =      significance level of the test.
%

% Press, WH, SA Teukolsky, WT Vetterling, and BP Flannery. 1992.  Numerical recipes in C: the
%   art of scientific programming.  Cambridge University Press.

% RE Strauss, 5/29/03

function pr = ksprob(lambda)
  eps1 = 1e-3;
  eps2 = 1e-8;
  fac = 2;
  pr = 0;
  termbf = 0;
  a2 = -2*lambda*lambda;
  
  for j = 1:100
    term = fac*exp(a2*j*j);
    pr = pr + term;
    if (abs(term)<=eps1*termbf | abs(term)<=eps2*pr)
      return;
    end;
    fac = -fac;
    termbf = abs(term);
  end;
  pr = 1;

  return;
  