% BINOFIT: Fits a binomial distribution B(N,p) to a discrete
%          distribution (X=0,1,2,...,N) by minimizing the mean-squared 
%          difference between data and binomial distribution.
%
%     Syntax: [N,p,expect,mse] = binofit(freq,{X})
%
%         freq =    vector of absolute frequencies.
%         X =       optional vector of integral values of X (of N)
%                     represented by the distribution [default: 0:1:N].
%         -------------------------------------------------------------
%         N,p =     parameters of best-fitting binomial distribution.
%         expect =  row vector of expected frequencies.
%         mse =     mean squared difference between data and best-fitting 
%         binomial distribution.
%

% RE Strauss, 6/7/96
%   2/13/01 - general improvements and change in output arguments.

function [N,p,expect,mse] = binofit(freq,X)
  if (nargin < 2) X = []; end;

  freq = freq(:);
  if (isempty(X))                       % Default X vector
    X = [0:(length(freq)-1)]';
  end;

  lenx = length(X);
  if (length(freq) ~= lenx)
    error('  BINOFIT: input vectors must be same size');
  end;

  if (min(X) > 0)                       % Pad X & freq with lower values 
    freq = [zeros(min(X),1); freq];     %   if necessary
  end;

  N = length(freq)-1;                   % Initial N at max observed X
  p = fmin('binofitf',0,1,[],N,freq);   % Optimize p
  mse = binofitf(p,N,freq);             % Initial MSE statistic

  mse_save = mse + 1;
  while (mse < mse_save)                % Until fit degrades,
    N_save = N;                         % Save current values
    p_save = p;
    mse_save = mse;

    N = N + 1;                          % Next N, next values
    p = fmin('binofitf',0,1,[],N,freq); % Optimize p
    mse = binofitf(p,N,freq);           % MSE statistic
  end;

  N = N_save;
  p = p_save;
  mse = mse_save;

  expect = sum(freq) * binopdf(X,N,p);  % Expected binomial counts

  return;
