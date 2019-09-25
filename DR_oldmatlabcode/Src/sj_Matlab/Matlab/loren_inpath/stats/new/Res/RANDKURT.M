% RANDKURT: Selects a random sample from a normal distribution and then
%           adjusts the sample to have values for the mean, standard
%           deviation, and kurtosis within sampling error of the
%           target values.
%
%     Usage: [X,m,s,kurt,sse] = randkurt([N|X],mu,sigma,g2,popl)
%
%           N -     sample size to be drawn from population (scalar)
%             or
%           X -     sample already drawn from population (column vector)
%           mu -    population mean
%           sigma - population standard deviation
%           g2 -    population kurtosis
%           popl -  boolean flag indicating whether to use the population
%                     parameters (=TRUE) or to sample estimates of them (=FALSE)
%                     [default=FALSE]
%           X -     random sample
%           m -     sample mean
%           s -     sample standard deviation
%           kurt -  sample kurtosis
%           sse -   final sse from fit
%

%     Finds the optimal value of the 'inverse' Box-Cox parameter lambda,
%     which specifies the transformation (in this case) from a normal
%     distribution to a skewed one.

% RE Strauss, 4/14/95
%   8/20/99 - changed plot colors for Matlab v5.

function [X,m,s,kurt,sse] = randkurt(N,mu,sigma,g2,popl)
  if (nargin < 6) 
    popl = 0;
  end;

  if (length(N)==1)                   % If scalar passed,
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
  end;

  if (popl)                           % Use population parameters
    samp_mean = mu;
    samp_stdev = sigma;
    samp_g2 = g2
  else                                % Else use sampling distribs
    se_mean =  sigma / sqrt(N);       % Standard errors
    se_stdev = sigma / sqrt(2*N);
    se_g2 =    sqrt(((24*N*(N-1)*(N-1))/((N-3)*(N-2)*(N+3)*(N+5))));

    samp_mean =  normrnd(mu,se_mean); % Random statistic values
    samp_stdev = normrnd(sigma,se_stdev);
    samp_g2 =    normrnd(g2,se_g2);
  end;

  figure(1),histgram(X);
  title('Initial Distribution');
kurt = kurtosis(X)

  lambda_range = [-1:.01:4];

  sse = zeros(length(lambda_range),2);
  i = 0;
  for l1 = lambda_range;
    i = i+1;
    sse(i,1) = l1;
    sse(i,2) = kurtopt(l1,g2,X);
  end;

  figure(2),
  plot(sse(:,1),sse(:,2),'k');

  options = foptions;
%  options(1) = 1;                    % Tabular display of results
  options(3) = 1e-6;                 % Termination tolerance on F(X)
  options(9) = 0;                    % No user-defined gradients
  options(14) = 1000;                % Number of iterations
%  lambda = fmins('skurtopt',init_lambda,options,[],samp_g2,X);  % Optimize

%  sse = kurtopt(lambda,sampg,X);

%  m = mean(X);
%  d = X-m;
%  dirdev = sign(d);                   % Directions of deviations from mean
%  if (abs(lambda) < eps)             % Transform kurtosis
%    X = m + dirdev * log(abs(d));
%  else
%    X = m + dirdev * (((abs(d)).^lambda)/lambda);
%  end;

%  s = std(X);                          % Transformed standard deviation
%  X = (samp_stdev/s)*X;                % Final rescale
%  m = mean(X);                         % Transformed mean
%  X = X-m+samp_mean;                   % Final rescale

%  figure(3),
%  histgram(X);
%  title('Final Distribution');

%  m = mean(X);                         % Statistics of transformed sample
%  s = std(X);
%  skew = skewness(X)

  return;

