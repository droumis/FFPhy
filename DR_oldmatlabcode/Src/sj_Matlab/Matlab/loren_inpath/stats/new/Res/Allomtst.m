% ALLOMTST: Tests for heterogeneity among within-group size vectors (PC1 
%           loadings) based on randomized vector correlations (randomly permuted 
%           group labels).  Assumes that the data have been log-transformed.
%           All groups are zero-centered before size vectors are estimated.
%             The test statistic for the overall test for homogeneity is the 
%           1s-complement of the mean squared vector correlations of all 
%           within-group size vectors from the pooled among-group size vector.  
%
%     Usage: [allomcf,rv,pr,prp,hgrps] = allomtst(logX,grps,{iter},{alpha})
%
%           logX =    [n x p] matrix of log-transformed morphometric data.
%           grps =    [n x 1] group-membership vector for k groups.
%           iter =    number of permuted iterations [default = 0].  If iter=0, 
%                       size vectors and vector correlations are provided but 
%                       statistical tests are not performed.
%           alpha =   optional overall alpha level for pairwise significance 
%                       tests [default = 0.05].
%           ------------------------------------------------------------------
%           allomcf = [p x k] matrix of within-group allometric coefficients.
%           rv =      [k x k] symmetric matrix of between-group vector 
%                       correlations for within-group PC1 vectors.
%           pr =      bootstrap probability for the null hypothesis that all 
%                       within-group size vectors are equal.
%           pr_pairs = [k x k] matrix  of pairwise probabilities (above diagonal) 
%                       and corresponding sequential Bonferroni significance 
%                       decisions (below diagonal).
%           hgrps =   [nhgrps x k] matrix specifying homogeneous subsets.
%

% Klingenberg, CP & R Froese. 1991. A multivariate comparison of allometric 
%     growth coefficients.  Syst Zool 40:410-419.

% RE Strauss, 1/15/00

function [allomcf,rv,pr,pr_pairs,hgrps] = allomtst(logX,grps,iter,alpha)
  if (nargin < 3) iter = []; end;
  if (nargin < 4) alpha = []; end;

  get_vectcorr = 0;
  get_probs = 0;
  get_hgrps = 0;

  if (nargout > 1)
    get_vectcorr = 1;
  end;
  if (nargout > 2)
    get_probs = 1;
  end;
  if (nargout  > 4)
    get_hgrps = 1;
  end;

  if (isempty(iter) | ~get_probs)
    iter = 0;
  end;
  if (~iter)
    get_probs = 0;
  end;
  if (isempty(alpha))
    alpha = 0.05;
  end;

  [N,P] = size(logX);

  [grpid,n] = uniquef(grps);
  if (any(n<P))
    disp('  ALLOMTST warning: At least one group size is less than');
    disp('                    number of variables');
  end;

  cX = grpcentr(logX,grps);               % Center data by group
  Sw = sizevect(cX,grps,1);               % Within-group size vectors
  allomcf = allom(Sw);                    % Allometric coefficients

  if (get_vectcorr)                       % Pairwise vector correlations
    rv = vectcorr(Sw);
  end;

  if (get_probs)
    r = trilow(rv);
    tsp_hat = 1 - r.*r;                   % Vector of pairwise test statistics
    prp = zeros(size(tsp_hat));             

    Sa = sizevect(logX,[],1);             % Pooled among-group size vector
    r = vectcorr(Sw,Sa);                  % Vect corrs of w-grp with a-grp vector
    ts_hat = 1 - mean(r.*r);              % Observed value of overall test statistic
    pr = 0;

    incr = 1./iter;                       % Probability increment
    for it = 1:iter                       % Iterate random permutations
      g = grps(randperm(N));                % Permute group labels
      cX = grpcentr(logX,g);                % Center data by group
      Sw = sizevect(cX,g,1);                % Within-group size vectors

      prv = vectcorr(Sw);                   % Pairwise tests
      r = trilow(prv);
      tsp = 1 - r.*r;
      i = find(tsp >= tsp_hat);
      if (~isempty(i))
        prp(i) = prp(i) + incr*ones(length(i),1);
      end;

      r = vectcorr(Sw,Sa);                  % Overall test
      ts = 1 - mean(r.*r);
      if (ts >= ts_hat)
        pr = pr + incr;
      end;
    end;

    [signif,critlim] = seqbonf(prp,alpha);
    pr_pairs = trisqmat([signif prp]);

    if (get_hgrps)
      hgrps = homosub(trisqmat(signif));
    end;
  end;

  return;
