% PAIRWISE: Evaluates pairwise post-hoc anovas, either parametric (ANOVA) or 
%           nonparametric (Kruskal-Wallace).
%
%     Usage: [pr,S,pr_pairs,S_pairs,hsets] = ...
%                 pairwise(x,grps,{KW},{iter},{hflag},{alpha})
%
%         x =         [n x 1] observations for a single variable.
%         grps =      [n x 1] classification variable for k groups.
%         KW =        optional boolean flag indicating that Kruskal-Wallace 
%                       tests are to be performed [default = 0 = anova].
%         iter =      optional number of randomization iterations [default=0].
%         hflag =     optional boolean flag indicating that homogeneous subsets 
%                       are to be evaluated at individual pairwise levels of 
%                       alpha, rather than by sequential Bonferroni critical 
%                       limits [default = 0 = Bonferroni limits].
%         alpha =     optional overall alpha level for pairwise significance 
%                       tests [default = max(0.05, pr)].
%         ----------------------------------------------------------------------
%         pr =        right-tailed probabilities.
%         S =         corresponding test-statistic value (F or H).
%         pr_pairs =  [k x k] matrix  of pairwise probabilities (above diagonal) 
%                       and corresponding sequential Bonferroni significance 
%                       decisions (below diagonal).
%         S_pairs =   corresponding pairwise test-statistic values (F or H).
%         hsets =     [nhsets x k] matrix specifying homogeneous subsets.
%

% RE Strauss, 9/12/00

function [pr,S,pr_pairs,S_pairs,hsets] = pairwise(x,grps,KW,iter,hflag,alpha)

  if (nargin < 3) KW = []; end;
  if (nargin < 4) iter = []; end;
  if (nargin < 5) hflag = []; end;
  if (nargin < 6) alpha = []; end;

  get_hsets = 0;
  if (nargout > 4)
    get_hsets = 1;
  end;

  default_alpha = 0.05;

  if (isempty(KW))
    KW = 0;
  end;
  if (isempty(iter))
    iter = 0;
  end;
  if (isempty(hflag))
    hflag = 0;
  end;
  if (isempty(alpha))
    alpha = default_alpha;
  end;

  if (min(size(x))>1)
    error('  PAIRWISE: data must comprise a vector');
  end;
  if (min(size(grps))>1)
    error('  PAIRWISE: classification variable must be a vector');
  end;
  
  N = length(x);
  if (length(grps)~=N)
    error('  PAIRWISE: data and classification vectors not of same length');
  end;

  x = x(:);                               % Convert input to col vectors
  grps = grps(:);
  if (KW)                                 % Convert to ranks for KW equivalent
    x = ranks(x);
  end;

  [grps,x] = sortmat(grps,x);             % Sort data by group

  [grp_id,freq] = uniquef(grps);          % Group identifiers & freqs
  ngrps = length(grp_id);                 % Number of groups
  freq = [0; cumsum(freq)];               % Cumulative frequencies

  if (ngrps<3)
    error('  PAIRWISE: need at least three groups for pairwise tests.');
  end;

  are_ranks = 0;
  if (KW)
    x = ranks(x);                           % Convert data to ranks
    are_ranks = 1;                          % Flag indicating ranks
  end;

  pr = NaN;
  pr_pairs = NaN*ones(ngrps,ngrps);  
  S_pairs = zeros(ngrps,ngrps);
  hsets = [];

  if (KW)                                   % Test statistic for all grps
    S = kwstat(x,grps,are_ranks);    
    [s,pr] = anova(x,grps);       
  else
    [S,pr] = anova(x,grps);
  end;

  for i = 1:(ngrps-1)                       % Pairwise statistic values
    gi = (freq(i)+1):freq(i+1);
    for j = (i+1):ngrps
      g = [gi (freq(j)+1):freq(j+1)];
      if (KW)
        S_pairs(i,j) = kwstat(x(g),grps(g));
        [s,pr_pairs(i,j)] = anova(x(g),grps(g));
      else
        [S_pairs(i,j),pr_pairs(i,j)] = anova(x(g),grps(g));
      end;
      S_pairs(j,i) = S_pairs(i,j);
      pr_pairs(j,i) = pr_pairs(i,j);
    end;
  end;

  if (iter)
    incr = 1./iter;
    pr = 0;
    pr_pairs = zeros(ngrps,ngrps);
    for it = 1:iter                         % Iterate the probabilities
      xr = x(randperm(N));                    % Randomly permute data
      if (KW)
        Sr = kwstat(xr,grps,are_ranks);         % Randomized KW statistic for all grps
      else
        Sr = anova(xr,grps);                    % Randomized anova for all grps
      end;
      if (Sr >= S)
        pr = pr+incr;
      end;

      for i = 1:(ngrps-1)
        gi = (freq(i)+1):freq(i+1);
        for j = (i+1):ngrps
          g = [gi (freq(j)+1):freq(j+1)];
          if (KW)
            Sr = kwstat(xr(g),grps(g));
          else
            Sr = anova(xr(g),grps(g));
          end;
          if (Sr >= S_pairs(i,j))
            pr_pairs(i,j) = pr_pairs(i,j)+incr;
            pr_pairs(j,i) = pr_pairs(i,j);
          end;
        end;
      end;
    end; % Iterations
  end;

  if (isempty(alpha))
    alpha = max([default_alpha, pr]);
  end;

  p = trilow(pr_pairs');                  
  if (hflag)
    signif = (p <= alpha);                % Pairwise significance levels
  else
    [signif,critlim] = seqbonf(p,alpha);  % Sequential Bonferroni adjustment
  end;
  pr_pairs = trisqmat([signif p]);

  if (get_hsets)
    hsets = homosub(trisqmat(signif));
  end;

  return;
