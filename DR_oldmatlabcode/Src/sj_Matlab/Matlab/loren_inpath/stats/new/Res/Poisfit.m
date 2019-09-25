% POISFIT: Fits a Poisson distribution P(lambda) to a discrete
%          distribution (X=0,1,2,...,N) by minimizing the probability of a
%          chi-square goodness-of-fit (power deviation) statistic.
%
%     Syntax: [lambda,expect,X2,df,pr] = poisfit(obs,X)
%
%         obs =     vector of absolute frequencies.
%         X =       optional vector of integral values of X (of N) 
%                     represented by the distribution [default: 0:1:N].
%         -------------------------------------------------------------
%         lambda =  parameter of best-fitting Poisson distribution.
%         stathat = corresponding test-statistic value.
%         exp =     expected values.
%         X2 =      chi-square value at best fit.
%         df =      degrees of freedom.
%         pr =      probability of the observed test-statistic.
%

% RE Strauss, 6/7/96
%   7/1/99 - modify to use GOODFIT, the power deviation chi-square statistic.

function [lambda,exp,X2,df,pr] = poisfit(obs,X)
  if (nargin < 2)                 % Default X vector
    X = [0:(length(obs)-1)]';
  elseif (size(X,2)>1)            % Convert input to column vectors
    X = X';
  end;
  if (size(obs,2)>1)
    obs = obs';
  end;

  lenx = length(X);
  if (length(obs) ~= lenx)
    error('  Error: input vectors must be same size');
  end;

  if (min(X) > 0)                 % Pad X & obs with lower values if necessary
    obs = [zeros(min(X),1); obs];
  end;

  lambda = fmin('poisfitf',0.1,max(X)-0.1,[],obs);   % Optimize lambda

  exp = sum(obs) * poisspdf(X,lambda);
  [X2,df,pr] = goodfit(obs,exp);

  return;
