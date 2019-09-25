% MORPHINT: Measures and assesses the significance of morphological integration 
%           (character suites).  Uses Leamy's (1994) modification of a measure of 
%           integration proposed by Cheverud, Rutledge & Atchley (1983:897) 
%           (adjusted variance among the eigenvalues), but defaults to the 
%           covariance matrix rather than the correlation matrix.  Eigenvalues 
%           are scaled to a sum of one before integration is estimated, 
%           consistent with sets of eigenvalues estimated from correlation 
%           matrices.  'Size-free' patterns of integration are assessed by 
%           ignoring the eigenvector of PC1, under the assumption that all size 
%           variation is described by PC1, and scaling the remaining eigenvalues 
%           to a mean of one.  The index varies from zero (no correlations among 
%           characters) to unity (singular covariance or correlation matrix).
%              Significance levels of individual character suites (H0: I=0) are 
%           assessed via a permutation test in which values of variables are 
%           independently randomized under the null hypothesis, and the 
%           integration index I is used as the test statistic.  
%              The test for heterogeneity among character suite index values 
%           (H0: I_1 = ... = I_k) is done via the same permutation test.
%
%     Usage: [I,Ici,pr_w,signif_w,pr_a,evals,pc1] = ...
%                         morphint(suites,X,{iter},{sizefree},{usecorr},{alpha})
%
%           suites =    grouping vector (length p) indicating s hypothetical 
%                         character suites.
%           X =         [n x p] data matrix for single group of specimens.
%           iter =      number of iterations for random-permutation test 
%                         [default = 0].
%           sizefree =  boolean variable indicating whether (=1) or not (=0) to 
%                         ignore PC1 in estimation of integration indices 
%                         [default = 0].
%           usecorr =   boolean variable indicating (=1) that principal 
%                         component analyses are to done on the correlation 
%                         matrices of characters [default = covariance matrices].
%           alpha =     critical significance level for sequential Bonferroni 
%                         tests [default = 0.05].
%           ---------------------------------------------------------------------
%           I =         [1 x s] vector of integration indices (in numerical 
%                         sequence of suite identifiers).
%           Ici =       [2 x s] matrix of confidence intervals for integration 
%                         indices, based on bootstrapping observations.
%           pr_w =      corresponding [1 x s] vector of statistical levels for 
%                         individual character suites, based on character-
%                         permutation test.
%           signif_w =  corresponding [1 x s] vector of boolean values indicating 
%                         whether character suites are statistically significant 
%                         from one another, based on a sequential Bonferroni test.
%           pr_a =      statistical level for test of homogeneity among character 
%                         suites.
%           evals =     [p x s] matrix of eigenvalues for each suite.
%           pc1 =       [p x s] matrix of PC1 loadings for each suite (as vector 
%                         correlations).
%

% RE Strauss, 12/11/98
%   12/16/98 - added sequential Bonferroni test for individual character suites, 
%                and permutation test for heterogeneity among suites.
%   6/2/01 -   changed measure of integration to I3; added bootstrapped confidence 
%                intervals.
%   7/12/01 -  export eigenvalues and PC1 loadings.

function [I,Ici,pr_w,signif_w,pr_a,evals,pc1] = ...
                                morphint(suites,X,iter,sizefree,usecorr,alpha)

  if (nargin < 3) iter = []; end;
  if (nargin < 4) sizefree = []; end;
  if (nargin < 5) usecorr = []; end;
  if (nargin < 6) alpha = []; end;

  get_pr_w = 0;                           % Check for output arguments,
  get_pr_a = 0;                           %   set flags
  get_signif_w = 0;
  get_ci = 0;

  if (nargout > 1)
    get_ci = 1;
  end;
  if (nargout > 2)
    get_pr_w = 1;
  end;
  if (nargout > 3)
    get_signif_w = 1;
  end;
  if (nargout > 4)
    get_pr_a = 1;
  end;

  if (isempty(sizefree))                  % Default argument values
    sizefree = 0;
  end;
  if (isempty(iter))
    iter = 0;
  end;
  if (isempty(usecorr))
    usecorr = 0;
  end;
  if (isempty(alpha))
    alpha = 0.05;
  end;

  [nobs,p] = size(X);
  nchars = length(suites);

  if (nchars ~= p)
    error('  MORPHINT: suite vector and data matrix not compatible');
  end;

  [suite_ids,npersuite] = uniquef(suites,1); % Sorted suite identifiers
  nsuites = length(suite_ids);

  I = zeros(1,nsuites);                   % Allocate output matrices
  Ici = [];
  pr_w = [];
  pr_a = [];
  signif_w = [];
  evals = zeros(max(npersuite),nsuites);
  pc1 = NaN*ones(p,nsuites);

  for s = 1:nsuites                       % Estimate integration indices
    v = find(suites==suite_ids(s));         % Vars for current suite
    x = X(:,v);                             % Extract suite from data matrix
    [I(s),evals(1:npersuite(s),s),pc1(v,s)] = ...
             morphinti(x,sizefree,usecorr); % Integration coefficient
  end;

  if (iter & get_ci)                      % Bootstrapped confidence intervals
    bI = zeros(iter,nsuites);          
    for it = 1:iter
      bX = bootsamp(X);
      for s = 1:nsuites
        bx = bX(:,suites==suite_ids(s));
        bI(it,s) = morphinti(bx,sizefree,usecorr);
      end;
    end;
    Ici = bootci(bI,I,[],1-alpha);
  end;

  if (iter & get_pr_w)                    % Probability via permutation of characters
    varI = var(I);                           % Variance among I index values

    Iperm = zeros(1,nsuites);
    pr_w = zeros(1,nsuites);
    pr_a = 0;
    incr = 1./iter;

    for it = 1:iter
      for nc = 1:nchars                       % Independently permute char values
        X(:,nc) = X(randperm(nobs),nc);
      end;

      for s = 1:nsuites
        x = X(:,suites==suite_ids(s));          % Extract suite of vars from data matrix
        Iperm(s) = morphinti(x,sizefree,usecorr);   % Integration coefficient
      end;
      
      i = find(Iperm >= I);                     % Within-suite test
      if (~isempty(i))
        pr_w(i) = pr_w(i) + incr;
      end;

      if (get_pr_a)
        varIperm = var(Iperm);
        if (varIperm >= varI)
          pr_a = pr_a + incr;
        end;
      end;
    end;

    if (get_signif_w)
      signif_w = seqbonf(pr_w',alpha)';
    end;
  end;

  return;

