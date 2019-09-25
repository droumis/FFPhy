% VARCOMP: Variance components from a two-level nested AVOVA with unequal
%          sample sizes (Sokal & Rohlf 1981, pp. 294-299).  
%          Optionally performs a randomization test of significance and/or 
%          calculates bootstrapped confidence intervals, both holding constant 
%          the numbers of individuals/group and replicates/individual.
%
%     Syntax: [table,F,high_table,low_table] 
%                 = varcomp(grps,data,{iter},{boot_iter},{CI_level})
%
%           grps = [N x 1] vector of group-membership identifiers
%                     for N individuals.
%           data = [N x max(NR)] matrix of data for N individuals and a
%                     maximum of NR replicates per individual;
%                     zeros = missing data for individuals having fewer than
%                     max(NR) replicates.
%           iter = number of randomization iterations to estimate the null 
%                     distributions (and thus significance levels) of the 
%                     F-statistics [default = 0].
%           boot_iter = number of bootstrap iterations (bootstrapping 
%                     replicate observations within groups [default = 0].
%           CI_level = confidence level for bootstrapped confidence intervals
%                     [default = 95].
%           -------------------------------------------------------------------
%           table = [varG dfG ssG msG
%                    varI dfI ssI msI
%                    varR dfR ssR msR
%                    varT dfT ssT msT]
%
%                   where  G = among-group effect
%                          I = among-individual, within-group effect
%                          R = among-replicate, within-individual effect
%                          T = total effect
%
%                          var = variance component
%                          df =  degrees of freedom
%                          ss =  sum of squares
%                          ms =  mean square
%
%           F =     [fGI dfG dfI prGI
%                    fIR dfI dfR prIR]
%
%                   where  f =  F-statistic using Satterthwaite's approx of msI
%                          df = numerator and denominator degrees of freedom
%                          pr = significance level, asymptotic or randomized
%
%           high_table = bootstrapped upper confidence limits of table
%           low_table =  bootstrapped lower confidence limits of table
%

% RE Strauss, 2/8/96

function [table,F,high_table,low_table] = varcomp(grps,data,iter,boot_iter,CI_level)
  TRUE = 1; FALSE = 0;

  % Set execution flags based on input & output arguments

  get_signif = FALSE;
  if (nargout > 1)
    get_signif = TRUE;
  end;

  randomize = FALSE;
  if (nargin > 2)
    if (iter > 0 & get_signif)
      randomize = TRUE;
    end;
  end;

  bootstrap = FALSE;
  if (nargin > 3 & nargout > 2)
    bootstrap = TRUE;
  end;

  if (nargin < 4)
    CI_level = 0.95;
  end;

  % Delete individuals for which all replicate observations missing

  nobs =  size(data,1);               % Number of individuals
  indx = [];
  for i = 1:nobs
    if (~any(data(i,:)))
      indx = [indx; i];
    end;
  end;
  if (length(indx)>0)
    data(indx,:) = [];
    grps(indx) = [];
    nobs = nobs - length(indx);
  end;

  % ANOVA and variance components

  grpid = uniquef(grps);              % Group identifiers
  ngrps = length(grpid);              % Number of groups

  G = design(grps);                   % Design matrix for groups

  d = (data>0);                       % Non-zero (non-missing) values
  rep = sum(d')';                     % Number of replicates per individual
  nrep = sum(rep);                    % Total replicates across individuals

  dfR = 0;
  dfI = 0;
  dfG = ngrps - 1;
  dfT = -1;

  ssy2 = 0; ssy3 = 0;
  N3 =   0;   N4 = 0;

  for k = 1:ngrps                     % For all groups,
    indx = find(G(:,k));                % Locate data for current group
    Xk = data(indx,:);                  % Stash in separate matrix
    nk = size(Xk,1);                    % Sample size for group
    rXk = rep(indx);                    % Reps/individual for group
    grprXk = sum(rXk);                  % Reps/group

    N3 = N3 + grprXk*grprXk;            % MSE coefficients
    N4 = N4 + (rXk'*rXk)/grprXk;

    dfR = dfR + sum(rXk-1);             % Increment degrees of freedom
    dfI = dfI + nk-1;
    dfT = dfT + sum(rXk);

    sXk = sum(Xk')';                    % Sum of reps/individual
    ssy2 = ssy2 + sum((sXk.*sXk)./rXk); % Sum of within-individual mean-squared reps
    sy = sum(sXk);                      % Sum of within-individual reps
    ssy3 = ssy3 + sy*sy/grprXk;         % Sum within-group mean-squared reps
  end; % for k=1:ngrps

  ssy1 = sum(sum(data.^2));
  CT = (sum(sum(data))^2)/(dfT+1);

  ssR = ssy1 - ssy2;                    % Sums of squares
  ssI = ssy2 - ssy3;
  ssG = ssy3 - CT;
  ssT = ssy1 - CT;

  msR = ssR / dfR;                      % Mean squares
  msI = ssI / dfI;
  msG = ssG / dfG;
  msT = ssT / dfT;

  N1 = dfT + 1;                         % MSE coefficients
  N2 = rep' * rep;

  np0 = (N4-N2/N1)/dfG;
  n0 =  (N1-N4)/dfI;
  nb0 = (N1-N3/N1)/dfG;

  varR = msR;                           % Variance components
  varI = (msI-msR)/n0;
  varG = (msG-msR-np0*varI)/nb0;
  varT = varR + varI + varG;

  table = [varG dfG ssG msG
           varI dfI ssI msI
           varR dfR ssR msR
           varT dfT ssT msT];

  F = [];
  high_table = [];
  low_table = [];

  % Asymptotic F-tests

  if (get_signif)                       % Asymptotic F-tests
    w2 = np0 / n0;                        % Satterthwaite's msI correction
    w1 = 1 - w2;
    msIp = w1*msR + w2*msI;
    dfIp = floor((msIp*msIp)/((w1*w1*msR*msR/dfR)+(w2*w2*msIp*msI/dfI)));

    fGI = msG / msIp;                     % F-statistics
    fIR = msIp / msR;

    if (~randomize)
      prGI = 1-fcdf(fGI,dfG,dfIp);        % Probabilities
      prIR = 1-fcdf(fIR,dfIp,dfR);

      F = [fGI dfG dfIp prGI
           fIR dfIp dfR prIR];
    end;
  end;  % Asymptotic F-tests

  % Randomized F-tests

  if (randomize)                        % Randomized F-tests
    BF = zeros(iter,2);

    for it = 1:iter
if (rem(it,10)==0)
  it
end;
      indx = data>0;                      % Non-zero replicate values
      vals = data(indx);
      p = randperm(length(vals));         % Random permutation
      X = data;
      X(indx) = vals(p);                  % Permute into data matrix

      ssy2 = 0; ssy3 = 0;
      N4 = 0;

      for k = 1:ngrps                     % For all groups,
        indx = find(G(:,k));                % Locate data for current group
        Xk = X(indx,:);                  % Stash in separate matrix
        nk = size(Xk,1);                    % Sample size for group
        rXk = rep(indx);                    % Reps/individual for group
        grprXk = sum(rXk);                  % Reps/group
        N4 = N4 + (rXk'*rXk)/grprXk;        % MSE coefficients
        sXk = sum(Xk')';                    % Sum of reps/individual
        ssy2 = ssy2 + sum((sXk.*sXk)./rXk); % Sum of within-individual mean-squared reps
        sy = sum(sXk);                      % Sum of within-individual reps
        ssy3 = ssy3 + sy*sy/grprXk;         % Sum within-group mean-squared reps
      end; % for k=1:ngrps

      msR = (ssy1 - ssy2) / dfR;            % Mean squares
      msI = (ssy2 - ssy3) / dfI;
      msG = (ssy3 - CT) / dfG;

      np0 = (N4-N2/N1)/dfG;                 % MSE coefficients
      n0 =  (N1-N4)/dfI;
      w2 = np0 / n0;                        % Satterthwaite's msI correction
      w1 = 1 - w2;
      msIp = w1*msR + w2*msI;
      dfIp = floor((msIp*msIp)/((w1*w1*msR*msR/dfR)+(w2*w2*msIp*msI/dfI)));

      fGIp = msG / msIp;                     % F-statistics
      fIRp = msIp / msR;
      BF(it,:) = [fGIp fIRp];
    end;

    BF = sort(BF);
    pr = randprob([fGI fIR],BF);

    F = [fGI dfG dfI pr(1)
         fIR dfI dfR pr(2)];
  end;  % Randomized F-tests

  % Bootstrapped confidence intervals

  if (bootstrap)                        % Bootstrapped confidence intervals
    for it = 1:boot_iter

    end;
  end;  % Bootstrapped confidence intervals

  return;
