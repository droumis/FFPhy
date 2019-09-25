% KSTEST1C: Randomized 1-sample Kolmogorov-Smirnov test of the goodness of fit
%           of a data distribution to a continuous probability distribution.
%           The probability distribution is given as a set of x,y coordinates
%           tracing the pdf, specified as [x=(min:incr:max),f(x)], which is
%           integrated by cubic spline interpolation.  The range of
%           the abscissa of the test is given by the range of the data; thus
%           the pdf must include at least the sample range.
%           Calculates either an MSE (mean squared difference) goodness-of-fit 
%           statistic or the KS statistic (max difference between distributions).
%
%     Usage: [stathat,pr,power] = kstest1c(X,tdist,{stat},{plotdist},{iter},{alpha})
%
%         X =       [n x 1] vector of data scores, from which cumulative data
%                     distribution is created.
%         tdist =   [m x 2] matrix of pdf coordinates of probability distribution, from 
%                     which cumulative reference distribution is created.
%         stat  =   optional flag indicating the statistic to be calculated:
%                     0 = MSE statistic [default]
%                     1 = KS statistic
%         plotdist = flag indicating whether (TRUE, =1) or not (FALSE, =0) to
%                     plot the observed and theoretical cumulative functions
%                     [default = TRUE].
%         iter  =   number of iterations; if iter=0 [default], then only the
%                   observed statistic value is returned.
%         alpha =   expected probability of Type I error [default = 0.05].
%         ----------------------------------------------------------------------
%         stathat = observed statistic value.
%         pr =      estimated pricance level.
%         power  =  estimated power level at given alpha.
%

% RE Strauss, 9/25/96
%   9/15/98 - permit return of chi-square goodness-of-fit statistic.
%   8/20/99 - changed plot colors for Matlab v5.

function [stathat,pr,power] = kstest1c(X,tdist,stat,plotdist,iter,alpha)
  if (nargin < 3) stat = []; end;
  if (nargin < 4) plotdist = []; end;
  if (nargin < 5) iter = []; end;
  if (nargin < 6) alpha = []; end;

  if (size(X,2)>1)                        % Convert row vectors to col vectors
    X = X';
  end;
  if (size(tdist,2)>2)
    tdist = tdist';
  end;

  [rt,ct] = size(tdist);
  if (ct~=2 | rt<2)
    error('  Invalid theoretical distribution');
  end;

  if (isempty(stat))                      % Default input-argument values
    stat = 0;
  end;
  if (isempty(plotdist))                  
    plotdist = 1;
  end;
  if (isempty(iter))
    iter = 0;
  end;
  if (isempty(alpha))
    alpha = 0.05;
  end;

  if (stat < 0 | stat > 2)
    error('  KSTEST1C: invalid statistic flag');
  end;

  get_pr = 0;
  get_power = 0;
  if (nargout>1 & iter>0)
    get_pr = 1;
  end;
  if (nargout>2 & iter>0)
    get_power = 1;
  end;

  X = sort(X);                            % Sort data values and get range
  rnge = [X(1) X(length(X))];

  % Convert the reference pdf to cubic-spine piecewise polynomials (pp-form),
  % and from there to the pp-form of the integral (row vector).
  % Hanselman & Littlefield, "Mastering Matlab", pp. 142-148.

  pp = spline(tdist(:,1),tdist(:,2));

  stathat = ks1cf(X,[],[],0,pp,stat,rnge); % Get observed statistic value

  if (plotdist)
    y = [0:1/(length(X)-1):1]';               % Cumulative increments

    Xt = linspace(rnge(1),rnge(2));
    t = spintgrl(pp,Xt);                      % Evaluate integral at data points
    t = t - t(1);                             % Anchor to 0 at left
    t = t / t(length(t));                     % Anchor to 1 at right

    leny = length(y);
    leny2 = 2*leny-1;
    Xy = zeros(leny2,1);
    Xx = Xy;

    Xx([1:2:leny2]) = X([1:leny]);            % Expand to staircase form
    Xx([2:2:(leny2-1)]) = X([2:leny]);
    Xy([1:2:leny2]) = y([1:leny]);
    Xy([2:2:(leny2-1)]) = y([1:(leny-1)]);

    figure;                                   % Plot cum functions
    plot(Xt,t,':k',Xx,Xy,'-k');
    putbnd([Xt Xx'],[t Xy']);
    legend('Null','Empirical');
    putxlab('X');
    putylab('Cumulative relative frequency');
  end;

  if (get_pr | get_power)                     % Randomize
    procs = [0, get_pr, get_power];
    user_null = 1;
    [ci,pr,power] = bootstrp('ks1cf',procs,iter,alpha,X,[],user_null,pp,stat,rnge);
  end;

  return;
