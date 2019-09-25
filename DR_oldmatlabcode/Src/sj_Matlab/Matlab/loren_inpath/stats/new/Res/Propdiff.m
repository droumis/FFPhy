% PROPDIFF: Tests directly for the difference in two proportions (as opposed to 
%           indirectly via a 2x2 contingency test), and provides a confidence 
%           interval on the difference and the probability under the null 
%           hypothesis that the true difference = 0.  The confidence interval 
%           and probability can be bootstrapped; if not, asymptotic normal 
%           approximations are used and a flag is returned indicating whether 
%           the normal approximation is reasonable.  
%
%     Usage: [d,s2d,ci,prob,holds] = 
%                               propdiff(num1,den1,num2,den2,{iter},{ci_level})
%
%         num1,den1 = numerator and denominator of first proportion.
%         num2,den2 = numerator and denominator of second proportion.
%         iter =      optional number of iterations for bootstrapped confidence 
%                       interval and probability [default = 0].
%         ci_level =  optional confidence level for interval [default = 0.95].
%         ----------------------------------------------------------------------
%         d =         observed difference in proportions (p1-p2).
%         s2d =       sampling variance of the difference.
%         ci =        [1x2] vector of lower and upper confidence bounds.
%         prob =      2-tailed probability under the null hypothesis of delta=0.
%         holds =     [3x1] vector of information concerning the normal approx:
%                       boolean flag indicating whether the normal approximation 
%                         can reasonably be expected to hold;
%                       minimum sample size for the first proportion;
%                       minimum sample size for the second proportion.
%

% McPherson,G. 1990. Statistics in scientific investigation: its basis, 
%   application, and interpretation.  Springer-Verlag.

% RE Strauss, 1/17/00
%   7/20/01 - allow for proportions of zero;
%               use abs(diff) rather than diff.

function [d,s2d,ci,prob,holds] = propdiff(num1,den1,num2,den2,iter,ci_level)
  if (nargin < 5) iter = []; end;
  if (nargin < 6) ci_level = []; end;

  get_prob = 0;
  get_holds = 0;
  if (nargout >= 3)
    get_prob = 1;
  end;
  if (nargout >= 5)
    get_holds = 1;
  end;

  if (isempty(iter))
    iter = 0;
  end;
  if (isempty(ci_level))
    ci_level = 0.95;
  end;

  if (ci_level > 1)
    ci_level = ci_level/100;
  end;
  alpha = 1-ci_level;

  if (num1>den1 | num2>den2)
    error('  PROPDIFF: numerators must be <= denominators.');
  end;

  p1 = num1/den1;
  p2 = num2/den2;
  n1 = den1;
  n2 = den2;
  d = abs(p1-p2);

  if (p1 < eps)
    p1 = eps;
  end;
  if (p2 < eps)
    p2 = eps;
  end;
  if (p1 > 1-eps)
    p1 = 1-eps;
  end;
  if (p2 > 1-eps)
    p2 = 1-eps;
  end;

  s2d = p1*(1-p1)/n1 + p2*(1-p2)/n2;
  sd = sqrt(s2d);
  holds = [];
  ci = [];
  prob = [];

  if (sd<eps)
    get_prob = 0;
  end;

  if (get_prob)
    if (iter)                               % Bootstrapped results
      x = makegrps([1 0 1 0],[num1 den1-num1 num2 den2-num2]);
      g = makegrps([1 2],[n1 n2]);
      [ci,prob] = bootstrp('propdifff',[1 1 0 0],iter,alpha,x,g);
      if (d<0)
        ci = -ci;
      end;

    else                                    % Normal approximation
      if (get_holds)
        p1 = max([p1,eps]);
        p2 = max([p2,eps]);
        p1 = min([p1,1-eps]);
        p2 = min([p2,1-eps]);
        n1min = 2/(3*(p1*(1-p1)).^2);
        n2min = 2/(3*(p2*(1-p2)).^2);
        if (n1>=n1min & n2>=n2min)
          holds = [1; n1min; n2min];
        else
          holds = [0; n1min; n2min];
        end;
      end;

      ci = norminv([alpha/2,1-(alpha/2)],d,sd);
      prob = 2*(1-normcdf(d/sd));
    end;
  end;

  return;

