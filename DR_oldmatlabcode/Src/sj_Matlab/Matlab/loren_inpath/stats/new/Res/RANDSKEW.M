% RANDSTRB: Selects a random sample from a normal distribution and then
%           adjusts the sample to have values for the mean, standard
%           deviation, and skewness within sampling error of the
%           target values.
%
%     Usage: [X,m,s,skew,sse] = randstrb([N|X],mu,sigma,g1,popl)
%
%           N -     sample size to be drawn from population (scalar)
%             or
%           X -     sample already drawn from population (column vector)
%           mu -    population mean
%           sigma - population standard deviation
%           g1 -    population skewness
%           popl -  boolean flag indicating whether to use the population
%                     parameters (=TRUE) or to sample estimates of them (=FALSE)
%                     [default=FALSE]
%           X -     random sample
%           m -     sample mean
%           s -     sample standard deviation
%           skew -  sample skewness
%           sse -   final sse from fit
%

%     Finds the optimal value of the 'inverse' Box-Cox parameter lambda,
%     which specifies the transformation (in this case) from a normal
%     distribution to a skewed one.

% RE Strauss, 4/12/95
%   8/20/99 - changed plot colors for Matlab v5.

function [X,m,s,skew,sse] = randstrb(N,mu,sigma,g1,popl)
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
    samp_g1 = g1
  else                                % Else use sampling distribs
    se_mean =  sigma / sqrt(N);       % Standard errors
    se_stdev = sigma / sqrt(2*N);
    se_g1 =    sqrt((6*N*(N-1))/((N-2)*(N+1)*(N+3)));

    samp_mean =  normrnd(mu,se_mean); % Random statistic values
    samp_stdev = normrnd(sigma,se_stdev);
    samp_g1 =    normrnd(g1,se_g1);
  end;

  figure(1),histgram(X);
  title('Initial Distribution');
skew = skewness(X)

  if (abs(samp_g1-skew) < eps)             % Initial lambda
    init_lambda = 1;
  elseif (samp_g1-skew > 0)
    init_lambda = 1.5;
  else
    init_lambda = 0.5;
  end;

  lambda_range = [-100:1:100];

  sse = zeros(length(lambda_range),2);
  i = 0;
  for l1 = lambda_range;
    i = i+1;
    sse(i,1) = l1;
    sse(i,2) = skewopt(l1,g1,X);
  end;

  figure(2),
  plot(sse(:,1),sse(:,2),'k');

  options = foptions;
%  options(1) = 1;                    % Tabular display of results
  options(3) = 1e-6;                 % Termination tolerance on F(X)
  options(9) = 0;                    % No user-defined gradients
  options(14) = 1000;                % Number of iterations
%  lambda = fmins('skewopt',init_lambda,options,[],samp_g1,X);  % Optimize

%  sse = skewopt(lambda,sampg,X);

%  m = mean(X);
%  if (abs(lambda) < eps)             % Transform skewness
%    X = log(X);
%  else
%    X = sign(lambda)*(X.^lambda)/(abs(m).^lambda);
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

