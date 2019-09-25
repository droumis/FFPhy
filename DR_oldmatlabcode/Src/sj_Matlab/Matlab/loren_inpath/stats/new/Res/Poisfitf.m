% POISFITF: Objective function for POISFIT.  Returns the MSE between the 
%           observed frequencies and the Poisson distribution specified by 
%           lambda.
%
%     Syntax: stathat = poisfitf(lambda,freq)
%
%         lambda =  parameters of a Poisson distribution.
%         freq =    vector of absolute frequencies corresponding to integral 
%                     values of X (of N)
%         -------------------------------------------------------------------
%         mse =     corresponding probability (test-statistic) value.
%

% RE Strauss, 2/9/97

function mse = poisfitf(lambda,freq)
  X = [0:(length(freq)-1)]';

  n = sum(freq);                        % Total sample size
  e = n * poisspdf(X,lambda);           % Expected poisson counts

  mse = (freq-e)'*(freq-e);

  return;

