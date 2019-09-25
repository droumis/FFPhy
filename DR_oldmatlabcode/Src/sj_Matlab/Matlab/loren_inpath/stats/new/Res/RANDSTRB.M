% RANDSTRB: Selects a random sample from a normal distribution and then
%           adjusts the sample to have values for the mean, standard
%           deviation, skewness, and kurtosis within sampling error of the
%           target values.
%
%     Usage: [X,m,s,skew,kurt,lambda1,lambda2,sse]
%                       = randstrb([N|X],mu,sigma,g1,g2,popl)
%
%           N -     sample size to be drawn from population (scalar)
%             or
%           X -     sample already drawn from population (column vector)
%           mu -    population mean
%           sigma - population standard deviation
%           g1 -    population skewness
%           g2 -    population kurtosis
%           popl -  boolean flag indicating whether to use the population
%                     parameters (=TRUE) or to sample estimates of them (=FALSE)
%                     [default=FALSE]
%           X -     random sample
%           m -     sample mean
%           s -     sample standard deviation
%           skew -  sample skewness
%           kurt -  sample kurtosis
%           sse -   final sse from fit
%

%     Finds the optimal value of the 'inverse' Box-Cox parameter lambda,
%     which specifies the transformation (in this case) from a normal
%     distribution to a skewed one.

% RE Strauss, 5/26/95
%   1/4/00 -  changed fminu() to fmins().

function [X,m,s,skew,kurt,lambda1,lambda2,sse] = randstrb(N,mu,sigma,g1,g2,popl)
  if (nargin < 6) popl = []; end;

  if (isempty(popl))
    popl = 0;
  end;

  if (isscalar(N))                    % If scalar passed,
    X = normrnd(100,2,N,1);           %   initial random-normal sample,
  else                                %     avoiding negative values
    X = N;                            % Else stash vector X
    [r,c] = size(X);
    if (c==1)                         % Transpose if X is row vector
      N = r;
    else
      X = X';
      N = c;
    end;
    if (min(X) <= 0)
      X = X-min(X)+1;                 % Shift positive if necessary
    end;
  end;

  if (popl)                           % Use population parameters
    samp_mean = mu;
    samp_std =  sigma;
    samp_g1 =   g1;
    samp_g2 =   g2;
  else                                % Else use sampling distribs
    se_mean = sigma / sqrt(N);       % Standard errors
    se_std =  sigma / sqrt(2*N);
    se_g1 =   sqrt((6*N*(N-1))/((N-2)*(N+1)*(N+3)));
    se_g2 =   sqrt(((24*N*(N-1)*(N-1))/((N-3)*(N-2)*(N+3)*(N+5))));

    samp_mean = normrnd(mu,se_mean); % Random statistic values
    samp_std =  normrnd(sigma,se_std);
    while (samp_std < 1)
      samp_std = normrnd(sigma,se_std);
    end;
    samp_g1 =   normrnd(g1,se_g1);
    samp_g2 =   normrnd(g2,se_g2);
    while (samp_g2 < 1)
      samp_g2 = normrnd(g2,se_g2);
    end;
  end;

  figure(1),histgram(X);
  title('Initial Distribution');

  init_m = mean(X);
  init_s = std(X);
  init_skew = skewness(X);
  init_kurt = kurtosis(X);

  options = foptions;                 % Optimization options
%  options(1) = 1;                      % Tabular display of results
%  options(3) = 1e-6;                    % Termination tolerance on F(X)
  options(9) = 0;                       % No user-defined gradients
%  options(14) = 1000;                   % Number of iterations

  sampg = [samp_g1 samp_g2 samp_std];
  init = [1 1];
  lambda = fmins('momopt',init,options,[],sampg,X);  % Optimize lambdas
  sse = momopt(lambda,sampg,X);       % Final error

  lambda1 = lambda(1);
  lambda2 = lambda(2);

%  X = skewtran(lambda1,X);
  X = kurttran(lambda2,X);
  X = skewtran(lambda1,X);

  s = std(X);                         % Transformed standard deviation
  X = (samp_std/s)*X;                 % Final rescale
  m = mean(X);                        % Transformed mean
  X = X-m+samp_mean;                  % Final rescale

  m = mean(X);                        % Statistics of transformed sample
  s = std(X);
  skew = skewness(X);
  kurt = kurtosis(X);

  figure(2), histgram(X);
  title('Final Distribution');

  [init_m init_s init_skew init_kurt]
  [m s skew kurt sse]
  
  return;

