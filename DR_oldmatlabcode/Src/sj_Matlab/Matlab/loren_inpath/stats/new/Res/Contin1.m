% CONTIN1: One-way discrete goodness-of-fit test, either asymptotic or randomized, 
%          based on specified observed and expected frequencies.  If expected 
%          frequencies are not given, the observed values are assumed to be 
%          equiprobable.
%            See also GOODFIT and GOODFITP.
%
%     Syntax: [pr,tot_x2,df,power,cell_x2,cell_pr] =
%                contin1(obs,{exp},{useG},{iter},{estdata},{alpha})
%
%          obs =      [r x c] matrix of observed counts.
%          exp =      [r x c] matrix of expected proportions or counts.
%          useG =     optional flag indicating whether or not
%                       the G-statistic is to be used in place of the conventional 
%                       X^2 [default = FALSE = 0].
%          iter =     optional number of randomization iterations [default=0].
%          estdata =  optional flag specifying whether or not
%                       the expected frequencies were estimated from the data;
%                       if true, the degrees of freedom for the asymptotic
%                       probability are decremented [default=0].
%          alpha =    Pr(Type I error), used for estimating power [default=0.05].
%          -----------------------------------------------------------------------
%          pr =       overall randomized significance level of observed table, 
%                       either asymptotic (if iter=0) or randomized.
%          totx2 =    overall observed value of statistic: chi-square or G
%                       (William's correction).
%          df =       degrees of freedom for asymptotic test.
%          power =    overall power of the test.
%          cell_x2 =  [r x c] matrix of observed chi-square contributions, by cell.
%          cell_pr =  [r x c] matrix of randomized chi-square probabilities, 
%                       by cell. 
%

% RE Strauss, 4/15/96
%   5/9/99 -   miscellaneous improvements, changes in defaults.
%   10/31/99 - removed Yates' correction for chi-square test 
%                (Overall, JE. 1980.  Psychol. Bull. 87:132-135).

function [pr,totx2,df,power,cell_x2,cell_pr] = contin1(obs,exp,useG,iter,estdata,alpha)
  if (nargin < 2) exp = []; end;          % Optional arguments
  if (nargin < 3) useG = []; end;
  if (nargin < 4) iter = []; end;
  if (nargin < 5) estdata = []; end;
  if (nargin < 6) alpha = []; end;

  if (isempty(exp))                       % Default input-argument values
    exp = (sum(sum(obs))/length(obs))*ones(size(obs));
  end;
  if (isempty(iter))
    iter = 0;
  end;
  if (isempty(useG))
    useG = 0;
  end;
  if (isempty(estdata))
    estdata = 0;
  end;
  if (isempty(alpha))
    alpha = 0.05;
  end;

  if (alpha > 1)
    alpha = alpha/100;
  end;

  find_power = 0;                     % Flags and default values for computations
  find_cellvals = 0;

  totx2 = [];
  df = [];
  power = [];
  cell_x2 = [];
  cell_pr = [];

  if (nargout>3 & iter>0)
    find_power = 1;
  end;
  if (nargout > 4)
    find_cellvals = 1;
  end;

  [r,c] = size(obs);                      % Check vector sizes
  k = max([r,c]);                         % Number of classes
  N = sum(obs);                           % Sample size

  if (all(exp) < 1)                       % Convert expected proportions to counts
    exp = exp*N;
  end;

  if (min([r,c]) > 1)
    error('  CONTIN1: use CONTIN2 or CONTIN3 for multi-way contingency tables');
  end;
  if (size(exp) ~= [r,c])
    error('  CONTIN1: observed and expected vectors must be same size');
  end;
  if (any(exp<=0))
    error('  CONTIN1: expected values must be positive');
  end;
  if (abs(sum(exp)-N)>1e-6)
    error('  CONTIN1: totals of expected and observed counts must be identical');
  end;

  if (r==k)                               % Transpose to row vectors
    obs = obs';
    exp = exp';
  end;

  if (useG)                               % Observed cell G or chi-square values
    cell_x2 = -2.*exp.*log(obs./exp);       % Max-likelihood statistic
  else
    cell_x2 = (obs-exp).^2./exp;            % Pearson chi-square statistic
  end;

  totx2 = sum(cell_x2);                   % Observed total chi-square
  if (useG)                               % William's correction for overall G
    q = 1 + ((k+1))/(6*N);
    totx2 = totx2 ./ q;
  end;

  Rtotx2 = zeros(iter,1);                 % Allocate randomization matrices
  if (find_cellvals)
    Rcell_x2 = zeros(iter,k);
  end;

  if (iter==0)                            % Asymptotic probabilities
    df = k-1;
    if (estdata)
      df= df-1;
    end;
    pr = 1-chi2cdf(totx2,df);

  else                                    % Randomized probabilities and power
    for it = 1:iter                         % Null distribution
      Robs = prbcount(exp./N,N);              % Randomize expected counts, 
                                              %   disallowing zero cells
      if (useG)                               % Convert counts to chi-square vals
        Robs = -2.*exp.*log(Robs./exp);
      else
        Robs = (Robs-exp).^2./exp;        
      end;

      tRobs = sum(Robs);
      if (useG)                               % William's correction for G
        tRobs = tRobs ./ q;
      end;

      Rtotx2(it) = tRobs;                     % Stash in randomization matrices
      if (find_cellvals)
        Rcell_x2(it,:) = Robs;
      end;
    end;

    Rtotx2 = sort(Rtotx2);
    pr = randprob(totx2,Rtotx2);            % Overall probability
    if (find_cellvals)
      Rcell_x2 = sort(Rcell_x2);
      cell_pr = randprob(cell_x2,Rcell_x2);   % Cell probabilities
    end;

    if (find_power)                         % Estimate power
      indx = iter*(1-alpha);                  % Critical value from null distribution
      low =  min(max(floor(indx),1),iter);
      high = max(min(ceil(indx),iter),1);

      if (high-low > 0)
        delta = (indx - low)/(high - low);
        critval = Rtotx2(low) + delta*(Rtotx2(high)-Rtotx2(low));
      else
        critval = Rtotx2(low);
      end;

      power = 0;
      incr = 1/iter;

      for it = 1:iter                         % Sampling distribution
        Robs = prbcount(obs./N,N);              % Randomize observed counts, 
                                                %   disallowing zero cells
        if (useG)                               % Convert counts to chi-square vals
          Robs = -2.*exp.*log(Robs./exp);
        else
          Robs = (Robs-exp).^2./exp;        
        end;

        tRobs = sum(Robs);
        if (useG)                               % William's correction for G
          tRobs = tRobs ./ q;
        end;

        if (tRobs >= critval)                % Power = proportion of sampled
          power = power + incr;                 %   chi-square values >= critical value
        end;                                    %   from null distribution
      end;
    end;
  end;

  return;
