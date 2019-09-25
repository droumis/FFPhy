% PROBZERO: Finds the curve describing the probablity of observing 0 of N
%           'successes', for the case in which no instances of a particular
%           kind of observation have been observed.
%
%     Syntax: [prob,maxprop] = probzero(N,{alpha},{lenprob},{noplot})
%
%           N =       sample size.
%           alpha =   desired probability cutoff-level at which the binomial
%                     parameter P is to be estimated [default = 0.05].
%           lenprob = number of function sampling points (=rows of output matrix
%                     'prob') to be reported for the probability curve
%                     [default = 50].
%           noplot =  flag to suppress plot [default = FALSE].
%           --------------------------------------------------------------------
%           prob =    2-column matrix of probabilities of observing zero
%                     'successes' (col 2) as a decreasing function of possible
%                     values of binomial parameter P (col 1).
%           maxprop = maximum proportion of 'successes' in population
%                     (binomial parameter P), for probability 1-alpha.
%                     Interpretation: if the true proportion of 'successes' in 
%                     the population were maxprop, then the observed results 
%                     would occur more than alpha-percent of the time.  The 
%                     interval [0, maxprop] is a (1-alpha)% confidence interval 
%                     in this sense.
%

% RE Strauss, 5/24/96
%   9/3/99 - changed plot colors for Matlab v5.

% Note: The estimation of minimum sample sizes for this problem is problematic: as
% more observations are sampled beyond N, the minimum sample size increases rather
% than stabilizes.

function [prob,maxprop] = probzero(N,alpha,lenprob,noplot)
  if (nargin < 2) alpha = []; end;
  if (nargin < 3) lenprob = []; end;
  if (nargin < 4) noplot = []; end;

  tol =  1e-6;                    % Estimation tolerance
  minp = 1e-4;                    % Minimum proportion P for function
  minprob = .005;                 % Minimum probability detected

  if (isempty(alpha))             % Default values
    alpha = 0.05;
  end;
  if (isempty(lenprob))
    lenprob = 50;
  end;
  if (isempty(noplot))
    noplot = 0;
  end;

  find_maxprop = 0;               % Output arguments
  if (nargout > 1)
    find_maxprop = 1;
  end;

  if (exp(log(eps).*N) < minprob)  % Maximum extent of probability curve
    maxp = 0;
    incr = 0.1;
    while (incr > tol)
      while (exp(log(1-(maxp+incr)).*N) > minprob)
        maxp = maxp + incr;
      end;
      incr = 0.1*incr;
    end;
  else
    maxp = 1;
  end;

  prob = zeros(lenprob,2);        % Allocate output matrix
  incr = maxp./lenprob;
  i = 0;

  for p = minp:incr:maxp          % Get probability curve
    i = i+1;
    prob(i,:) = [p, exp(log(1-p).*N)];
  end;

  if (find_maxprop)               % Find max proportion at prob alpha
    i = find(prob(:,2)>alpha);
    maxprop = prob(max(i),1);
    incr = prob(max(i)+1,1) - maxprop;

    while (incr > tol)
      while (exp(log(1-(maxprop+incr)).*N) > alpha)
        maxprop = maxprop + incr;
      end;
      incr = 0.1*incr;
    end;
    maxprop
  end;

  if (~noplot)
    clf;
    plot(prob(:,1),prob(:,2),'k');
    xlabel('Proportion of successes in population');
    ylabel('Probability of observing zero successes');
    puttitle(sprintf('N = %d',N));
    if (find_maxprop)
      hold on;
      plot([0,maxprop],[alpha,alpha],'--k');
      plot([maxprop,maxprop],[0,alpha],'--k');
      hold off;
    end;
  end;

  return;
