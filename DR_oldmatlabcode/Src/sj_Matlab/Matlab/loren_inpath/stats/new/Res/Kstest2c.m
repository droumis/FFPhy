% KSTEST2C: Kolmogorov-Smirnov two-sample test of identity of two data
%           distributions, for continuous data, with optional randomization.
%           Calculates either the KS-statistic (max difference between distributions)
%           or an MSE-statistic (mean squared difference).  Returns either the 
%           asymptotic probability (for the KS-statistic only) or a randomized 
%           probability.
%
%     Syntax: [pr,stathat,power] = kstest2c(X,grps,stat,iter,alpha,doplot)
%
%         X1 =      data vector (length n1+n2=N).
%         grps =    corresponding group-membership vector.
%         stat  =   flag indicating the statistic to be calculated:
%                     0 = KS statistic [default],
%                     1 = MSE statistic.
%         iter  =   number of iterations; if iter=0 [default], then only the
%                     observed statistic value is returned.
%         alpha =   expected probability of Type I error [default=0.05].
%         doplot =  boolean flag indicating that plot of cumulative distributions
%                     is to be produced [default=0].
%         -----------------------------------------------------------------------
%         pr =      estimated significance level.
%         stathat = observed statistic value.
%         power  =  estimated power level at given alpha.
%

% Two-tailed probability level of Dmax is estimated using the first three terms 
%   of the Smirnov (1948) formula.  See SPSS 7.5 Statistical Algorithms.
% Smirnov, NV. 1948. Table for estimating the goodness of fit of empirical 
%   distributions.  Ann. Math. Stat. 19:279-281.

% RE Strauss, 9/26/96
%   5/10/99 - miscellaneous improvements (modified from KSTEST2).
%   5/28/03 - addition of "if (~nargin)" statement.

function [pr,stathat,power] = kstest2c(X,grps,stat,iter,alpha,doplot)
  if (~nargin) help kstest2c; return; end;

  if (nargin < 3) stat = []; end;
  if (nargin < 4) iter = []; end;
  if (nargin < 5) alpha = []; end;
  if (nargin < 6) doplot = []; end;

  if (isempty(stat))   stat = 0; end;
  if (isempty(iter))   iter = 0; end;
  if (isempty(alpha))  alpha = 0.05; end;
  if (isempty(doplot)) doplot = 0; end;

  if (alpha > 1)
    alpha = 0.01*alpha;
  end;

  get_power = 0;                            % Get power by randomization only
  if (nargout>2 & iter>0)
    get_power = 1;
  end;

  [g,f] = uniquef(grps,1);
  if (length(g) > 2)
    error('KSTEST2C: this is a 2-sample test.');
  end;

  stathat = kstst2cf(X,grps,[],[],stat,doplot);    % Get observed statistic value

  pr = 0;
  if (iter==0 & stat==1)                    % Get probability of Dmax
    z = stathat * sqrt(prod(f)/sum(f));
    if (z < 0.27)
      pr = 1;
    elseif (z < 1)
      Q = exp(-1.233701*z^(-2));
      pr = 1 - 2.506628*(Q+Q^9+Q^25)/z;
    elseif (z < 3.1)
      Q = exp(-1.233701*z^(-2));
      pr = 2*(Q-Q^4+Q^9-Q^16);
    else
      pr = 0;
    end;
  end;

  if (iter>0)                       % Randomize
    procs = [0 1 get_power];
    [ci,pr,power] = bootstrp('kststcf',procs,iter,alpha,X,grps,0,stat);
  end;

  return;
