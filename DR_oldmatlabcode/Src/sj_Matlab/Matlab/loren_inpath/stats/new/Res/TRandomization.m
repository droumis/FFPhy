% TRandomization: randomized t-test.
%
%     Usage: [that,n,pvalue] = trandomization(x,g,{iter})
%
%         x =      [n x 1] vector of data for a single variable.
%         g =      corresponding group-membership vector.
%         iter =   number of randomization iterations [default = 1000].
%         -----------------------------------------------------------
%         that =   value of the t-statistic from the data.
%         n =      vector of sample sizes.
%         pvalue = significance level of the test.
%

function [that,n,pvalue] = trandomization(x,g,iter)
  if (nargin<3)
    iter = [];
  end;
  
  if (isempty(iter))
    iter = 1000;
  end;

  [grpid,n] = uniquef(g);
  if (length(grpid)~=2)
    error('  TRandomization: must be exactly 2 groups.');
  end;
  
  that = tstatistic(x,g);         % T-hat from the data
  N = length(x);                  % Total sample size
  
  tnull = zeros(iter,1);          % Allocated the vector for the null distribution
  pvalue = 0;
  for it = 1:iter                 % Randomization loop
    r = randperm(N);                % Random permutation
    g = g(r);                       % Rearrange the grouping vector
    t = tstatistic(x,g);            % Randomized t statistic
    tnull(it) = t;                  % Stash the value of t
  end;
  
  figure;
  histgram(tnull);

  that = abs(that);                       % Use absolute value for right tail
  pvalue_right = sum(tnull>=that)/iter;   % Right-tail probability
  pvalue_left = sum(tnull<=that)/iter;    % Left-tail probability
  pvalue = pvalue_left + pvalue_right;    % Total probability

  return;
  
