% RANDSPLT: Selects a random sample from a normal distribution and then
%           adjusts the sample to have values for the mean, standard
%           deviation, skewness, and kurtosis within sampling error of the
%           target values.
%
%     Usage: [lambda1,lambda2,sse] = randsplt([N|X],mu,sigma,g1,g2,popl)
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
%           lambda1 - vector of lambda values for skewness
%           lambda2 - vector of lambda values for kurtosis
%           sse -   lambda1*lambda2 matrix of sse values
%

%     Finds the optimal value of the 'inverse' Box-Cox parameter lambda,
%     which specifies the transformation (in this case) from a normal
%     distribution to a skewed one.

% RE Strauss, 4/18/95

function [lambda1,lambda2,sse] = randsplt(N,mu,sigma,g1,g2,popl)
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
    samp_g1 =   normrnd(g1,se_g1);
    samp_g2 =   normrnd(g2,se_g2);
  end;

%  figure(1),histgram(X);
%  title('Initial Distribution');
skew = skewness(X);
kurt = kurtosis(X);

lambda1 = [0:2:100];
lambda2 = [0.1:0.06:3];

%  sse = zeros(length(lambda2),length(lambda1));
  sse = zeros(length(lambda1)*length(lambda2),3);
%  i = 0;
  k = 0;
  for l1 = lambda1;
%    i = i+1;
%    j = 0;
    for l2 = lambda2;
%      j = j+1;
      k = k+1;
%      sse(j,i) = min(20, momopt([l1,l2],[g1,g2,samp_std],X) );
      sse(k,[1:2]) = [l1 l2];
      sse(k,3) = min(20, momopt([l1,l2],[g1,g2,samp_std],X) );
    end;
  end;

%  figure, surf(lambda1,lambda2,sse);
%  axis([min(lambda1) max(lambda1) min(lambda2) max(lambda2) 0 max(max(sse))]);
%  view(-110,40);
%  figure, contour(lambda1,lambda2,sse,[.1 .2 .5 1 2 5 10]);
%  axis([min(lambda1) max(lambda1) min(lambda2) max(lambda2)]);

  return;

