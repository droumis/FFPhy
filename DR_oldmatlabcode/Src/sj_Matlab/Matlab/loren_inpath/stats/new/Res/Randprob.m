% RANDPROB: Finds the right-tailed probability level of an observed
%           statistic value from a sorted distribution of randomized
%           statistics.  If the observed statistic value is a scalar,
%           the distribution is expected to be a vector of length
%           'rand_iter'.  If a vector of p observed statistic values
%           is passed, the corresponding distributions are expected as
%           the columns of a (rand_iter x p) matrix.
%
%           If the observed value is less than the miminum randomized value,
%           the probability returned is 1/length(distribution), which is an
%           upper bound.
%           If the observed value is greater than the miminum randomized value,
%           the probability returned is 1-(1/length(distribution)), which is a
%           lower bound.
%
%     Usage: prob = randprob(s,distrib,{doplot})
%
%           s =       [1 x p] row vector of observed statistic values.
%           distrib = [rand_iter x p] matrix of randomized statistics.
%           --------------------------------------------------------------------
%           prob =    [1 x p] row vector of probabilities.
%           doplot =  optional boolean flag indicating, if true, that plot(s) of 
%                       distribution and critical value is to be produced 
%                       [default = 0].
%

% RE Strauss, 3/16/97
%   2/27/01 - added check for sorted cols; added plot option.
%   3/3/01 -  changed histgram from absolute to relative frequencies.

function prob = randprob(s,distrib,doplot)
  if (nargin < 3) doplot = []; end;

  if (isempty(doplot))
    doplot = 0;
  end;

  p = length(s);
  [iter,q] = size(distrib);
  if (p ~= q)
    disp('  RANDPROB: matrix sizes incompatible.');
    return;
  end;

  prob = zeros(1,p);                    % Allocate return vector

  for v=1:p;                            % Cycle thru variables
    S = s(v);                           % Extract statistic & distrib
    D = distrib(:,v);

    notsorted = sum(D(2:iter)>=D(1:iter-1));  % Sort if not sorted
    if (notsorted)
      D = sort(D);
    end;

    if (S <= D(1))
      Pr = 1-(1/iter);
    elseif (S > D(iter))
      Pr = 1/iter;
    elseif (S == D(iter))
      i = iter;
      while(S == D(i-1))
        i = i-1;
      end;
      Pr = 1-((i-1)/iter);
    else
      lower = 0;                          % Initialize lower & upper limits
      upper = iter+1;
      while ((upper-lower) > 1)           % Bisection search
        mid = floor((upper+lower)/2);
        if (S < D(mid))
          upper = mid;
        else
          lower = mid;
        end;
      end;

      while (lower>1 & (D(lower)==D(upper) | S==D(lower))) % Slide to left of series of
        lower = lower-1;                                   %   equal values
        upper = upper-1;
      end;

      f = (D(upper)-S)/(D(upper)-D(lower)); % Linearly interpolate
      Pr = (iter-upper+f+1)/iter;
    end;
    prob(v) = Pr;

    if (doplot)
      figure;
      [npts,bars] = histgram(D,[],[],[],[],[],'rel');
      
      hold on;
      maxbar = max(bars(:,2));
      plot([S S],[0 1.05*maxbar],'k');
      hold off;
    end;
  end;

  return;
