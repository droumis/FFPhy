% BINOFITF: Objective function for binofit(): the MSE statistic value for the 
%           goodness of fit of a data distribution Xto a binomial distribution
%           specified by parameters N,p.
%
%     Syntax: stathat = binofitf(p,N,freq)
%
%         p,N =     parameters of a binomial distribution.
%         freq =    vector of absolute frequencies corresponding to integral
%                     values of X (of N)
%         -------------------------------------------------------------------
%         stathat = corresponding test-statistic value.
%

% RE Strauss, 6/7/96

function stathat = binofitf(p,N,freq)
  d = (N+1) - length(freq);             % Check whether observed vector is too short
  if (d > 0)                            % If so, pad with zeros
    freq = [freq; zeros(d,1)];
  end;
  X = [0:N]';

  n = sum(freq);                        % Total sample size
  e = n * binopdf(X,N,p);               % Expected binomial counts

  stathat = (freq-e)'*(freq-e);

  return;

