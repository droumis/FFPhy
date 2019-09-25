% KSTEST2: Randomized Kolmogorov-Smirnov test of identity of two data
%          distributions, for continuous data.
%          Calculates either the KS-statistic (max difference between distributions)
%          or an MSE-statistic (mean squared difference).
%
%     Syntax: [stathat,signif,power] = kstest2(X1,X2,{stat},{iter},{alpha})
%
%         X1 =    vector (length n1) of scores for group 1.
%         X2 =    vector (length n2) of scores for group 2.
%         stat  = flag indicating the statistic to be calculated:
%                   0 = KS statistic [default],
%                   1 = MSE statistic.
%         iter  = number of iterations; if iter=0 [default], then only the
%                   observed statistic value is returned.
%         alpha = expected probability of Type I error [default = 0.05].
%         ----------------------------------------------------------------
%         stathat = observed statistic value.
%         signif =  estimated significance level.
%         power  =  estimated power level at given alpha.
%

% RE Strauss, 9/26/96
%   9/19/99 - updated handling of default input arguments.

function [stathat,signif,power] = kstest2(X1,X2,stat,iter,alpha)
  if (nargin < 3) stat = []; end;
  if (nargin < 4) iter = []; end;
  if (nargin < 5) alpha = []; end;

  if (size(X1,2)>1)                       % Convert row vectors to col vectors
    X1 = X1';
  end;
  if (size(X2,2)>1)
    X2 = X2';
  end;

  if (isempty(iter))                         % Check input/output arguments
    iter = 0;
  end;
  if (isempty(alpha))
    alpha = 0.05;
  end;

  get_signif = 0;
  get_power = 0;
  if (nargout > 1)
    get_signif = 1;
  end;
  if (nargout > 2)
    get_power = 1;
  end;

  if (isempty(stat))
    stat = 0;
  end;
  if (stat < 0 | stat > 1)
    error('  Error: invalid statistic flag');
  end;
  if (iter == 0)
    get_signif = 0;
    get_power = 0;
  end;

  n1 = length(X1);
  n2 = length(X2);

  X = [X1;X2];                          % Concatenate data
  grps = [zeros(n1,1);ones(n2,1)];      % Identify groups

  stathat = kstest2f(X,grps,[],[],stat);   % Get observed statistic value

  if (get_signif)                       % Randomize
    procs = [0, get_signif, get_power];
    [ci,signif,power] = bootstrp('kstest2f',procs,iter,alpha,X,grps,0,stat);
  end;

  return;
