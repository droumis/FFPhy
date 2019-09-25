% BOOTPROB: Finds the interpolated right-tailed probability level of an
%           observed statistic value from a sorted distribution of
%           bootstrapped statistics.  If the observed statistic value is
%           a scalar, the distribution is expected to be a vector of length
%           'boot_iter'.  If a vector of 'p' observed statistic values
%           is passed, the corresponding distributions are expected as
%           the columns of a (boot_iter x p) matrix.
%
%     Usage: prob = bootprob(s,distrib)
%
%           s =       [1 x p] row vector of observed statistic values.
%           distrib = [boot_iter x p] matrix of bootstrapped statistics.
%           ---------------------------------------------------------------
%           prob =    [p x 1] column vector of probabilities.
%

% RE Strauss, 1/30/96
%   12/3/00 - sort distributions if necessary.

function prob = bootprob(s,distrib)

  p = length(s);
  [boot_iter,q] = size(distrib);
  if (p ~= q)
    disp('  BOOTPROB: matrix sizes incompatible');
    return;
  end;

  prob = zeros(p,1);                    % Allocate return vector
  bsort = issorted(distrib);            % Determine whether distributions are sorted

  for var=1:p;                          % Cycle thru variables
    S = s(var);                         % Extract statistic & distrib
    D = distrib(:,var);

    if (~bsort(var))
      D = sort(D);
    end;

    if (S <= D(1))
      Pr = 1;
    elseif (S > D(boot_iter))
      Pr = 0;
    else
      lower = 0;                          % Initialize lower & upper limits
      upper = boot_iter+1;
      while ((upper-lower) > 1)           % Bisection search
        mid = floor((upper+lower)/2);
        if (S > D(mid))
          lower = mid;
        else
          upper = mid;
        end;
      end;
      f = (D(upper)-S)/(D(upper)-D(lower)); % Linearly interpolate
      Pr = (boot_iter-upper+f+1)/boot_iter;
    end;
    prob(var) = Pr;
  end;
  return;
