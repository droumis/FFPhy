% KSTEST1D: Randomized 1-sample Kolmogorov-Smirnov test of the goodness of fit
%           of a data distribution to a discrete probability distribution.
%           The data and reference distributions are given as correspondng vectors of 
%           frequencies representing the discrete 'bins' of the distributions.
%           Calculates either an MSE-statistic (mean squared difference),
%           the KS-statistic (max difference between distributions), or the 
%           Fisher X2 statistic.
%
%     Syntax: [stathat,signif,power] = kstest1d(X,Xf,Rf,{stat},{plotdist},{iter},{alpha})
%
%         X =       [n x 1] vector of abscissa values of discrete distributions.
%         Xf =      [n x 1] vector of data frequencies (absolute counts), one per abscissa.
%         Rf =      [n x 2] matrix of corresponding freqencies (absolute or relative) for 
%                     reference distribution.
%         stat  =   flag indicating the statistic to be calculated:
%                     0 = MSE statistic [default],
%                     1 = KS statistic,
%                     2 = X2 statistic.
%         plotdist = flag indicating whether (TRUE = 1) or not (FALSE = 0) to 
%                     plot the observed and theoretical cumulative functions 
%                     [default = TRUE].
%         iter  =   number of iterations; if iter=0 [default], then only the
%                   observed statistic value is returned.
%         alpha =   expected probability of Type I error [default = 0.05].
%         ----------------------------------------------------------------------
%         stathat = observed statistic value.
%         pr =      estimated significance level.
%         power  =  estimated power level at given alpha.
%

% RE Strauss, 9/26/96

function [stathat,pr,power] = kstest1d(X,Xf,Rf,stat,plotdist,iter,alpha)
  TRUE = 1; FALSE = 0;

  if (nargin < 4)
    stat = [];
  end;
  if (nargin < 5)
    plotdist = [];
  end;
  if (nargin < 6)
    iter = [];
  end;
  if (nargin < 7)
    alpha = [];
  end;

  if (size(X,2)>1)                        % Convert row vectors to col vectors
    X = X';
  end;
  if (size(tdist,2)>2)
    tdist = tdist';
  end;

  if (isempty(stat))                      % Default input-argument values
    stat = 0;
  end;
  if (isempty(plotdist))                  
    plotdist = TRUE;
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

  get_pr = FALSE;
  get_power = FALSE;
  if (nargout>1 & iter>0)
    get_pr = TRUE;
  end;
  if (nargout>2 & iter>0)
    get_power = TRUE;
  end;

  
  Xval = makegrps(X,Xf);                      % Expand freqs into observatons
  if (max(Rf)<1)                              % If reference freqs are relative, make absolute


  stathat = ks1df(X,[],[],tdist,stat,rnge);   % Get observed statistic value

  if (plotdist)
    X = sort(X);                              % X = sorted data sample
    y = [0:1/(length(X)-1):1]';               % Cumulative increments

    leny = length(y);
    leny2 = 2*leny-1;
    Xy = zeros(leny2,1);
    Xx = Xy;

    Xx([1:2:leny2]) = X([1:leny]);            % Expand to staircase form
    Xx([2:2:(leny2-1)]) = X([2:leny]);
    Xy([1:2:leny2]) = y([1:leny]);
    Xy([2:2:(leny2-1)]) = y([1:(leny-1)]);

    t = tdist(:,1);                           % Theor cdf
    ft = tdist(:,2);                          % Increments

    lent = length(t);
    lent2 = 2*lent-1;
    Xt = zeros(lent2,1);
    Xft = Xt;

    Xt([1:2:leny2]) = t([1:leny]);            % Expand to staircase form
    Xt([2:2:(leny2-1)]) = t([2:leny]);
    Xft([1:2:leny2]) = ft([1:leny]);
    Xft([2:2:(leny2-1)]) = ft([1:(leny-1)]);

    clf;                                      % Plot cum functions
    hold on;
    plot(Xx,Xy);
    plot(Xt,Xft);
    putxlab('X');
    putylab('Cumulative relative frequency');
    hold off;
  end;

  if (get_pr)                       % Randomize
    procs = [0, get_pr, get_power];
    [ci,pr,power] = bootstrp('ks1df',procs,iter,alpha,X,[],tdist,stat,rnge);
  end;

  return;
