% GOODFIT:  Power divergence goodness-of-fit test (Cressie & Read 1984), 
%           indexed by the parameter lambda.  Pearson's chi-squared test 
%           (lambda=1), the likelihood-ratio G test (lambda=0), and the 
%           Freeman-Tukey test (lambda=-0.5) are all members of the power 
%           divergence family of test statistics.  Cressie & Read (1984) 
%           recommend lambda=2/3 for best approximation to a chi-squared 
%           sampling distribution.
%
%     Usage: [X2,df,pr,X2cell] = goodfit(obs,exp,{df},{lambda})
%
%           obs =     matrix of observed counts.
%           exp =     matrix of correpsonding expected counts.
%           df =      optional degrees of freedom [default = nbins-1].
%           lambda =  optional power-divergence parameter [default = 2/3].
%           --------------------------------------------------------------
%           X2 =      total test statistic value.
%           df =      degrees of freedom.
%           pr =      right-tailed chi-square probability.
%           X2cell =  test-statistic values for individual cells.
%

% Cressie, N & TRC Read. 1984. Multinomial goodness-of-fit tests.  
%   J Roy Stat Soc 46:440-464.
% Young, LJ & JH Young. 1998. Statistical Ecology.  Kluwer Academic Publs.

% RE Strauss, 7/1/99

function [X2,df,pr,X2cell] = goodfit(obs,exp,df,lambda)
  if (nargin < 3) df = []; end;
  if (nargin < 4) lambda = []; end;

  if (size(obs)~=size(exp))
    error('GOODFIT: matrices of obs & exp counts not compatible');
  end;

  if (isempty(lambda))
    lambda = 2/3;
  end;

  [r,c] = size(obs);                % Vector or matrix?
  isvector = 0;
  if (min([r,c]) == 1)                % If vector,
    isvector = 1;                     %   set flag
    if (r==1)                         %   transpose row to column vectors
      obs = obs';
      exp = exp';
      k = c;
    else
      k = r;
    end;
  end;

  if (isempty(df))                  % Determine df
    if (isvector)
      df = k-1;
    else
      df = (r-1)*(c-1);
    end;
  end;

  i = find(exp<5);                  % Check for small expected cells
  if (~isempty(i))
    if (isvector)
      while (~isempty(i) & df>1)      % Combine cells
        [exp,obs] = sortmat(-exp,obs);
        exp = -exp;
        ln = length(exp);
        exp(ln-1) = exp(ln-1)+exp(ln);
        exp(ln) = [];
        obs(ln-1) = obs(ln-1)+obs(ln);
        obs(ln) = [];
        df = df-1;
        i = find(exp<5);
      end;
%      if (~isempty(i))
%        disp('GOODFIT warning: small expected cell counts in matrix');
%      end;
    else
      disp('GOODFIT warning: small expected cell counts in matrix');
    end;
  end;

  X2cell = (2/(lambda*(lambda+1))) * obs.*((obs./exp).^lambda-1);
  X2 = sum(sum(X2cell));
  pr = 1 - chi2cdf(X2,df);

  return;
