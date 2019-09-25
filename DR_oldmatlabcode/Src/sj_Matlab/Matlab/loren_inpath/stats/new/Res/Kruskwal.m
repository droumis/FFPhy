% KRUSKWAL: Randomized Kruskal-Wallis test for average differences among 
%           populations.  If iter>0, both the overall and the pairwise tests are 
%           randomized.   If iter=0, probabilities are estimated by an anova on 
%           ranks of data.  Finds all possible homogeneous subsets consistent 
%           with pairwise tests.
%
%     Usage: [pr,H,grp_medians,pr_pairs,H_pairs,hsets] = ...
%               kruskwal(x,grps,{iter},{pairwise},{hflag},{alpha})
%
%         x =           [n x 1] observations for a single variable.
%         grps =        [n x 1] classification variable for k groups.
%         iter =        optional number of randomization iterations [default=0].
%         pairwise =    optional boolean flag indicating, if true, that all 
%                         pairwise multiple comparisons are to be done and 
%                         homogeneous subsets are to be found [default=0].
%         hflag =       optional boolean flag indicating that homogeneous subsets 
%                         are to be evaluated at individual pairwise levels of 
%                         alpha, rather than by sequential Bonferroni critical 
%                         limits [default = 0 = Bonferroni limits].
%         alpha =       optional overall alpha level for pairwise significance 
%                         tests [default = max(0.05, pr)].
%         ------------------------------------------------------------------------
%         pr =          right-tailed KW probabilities.
%         H =           corresponding KW statistic value.
%         grp_medians = [k x 1] vector of medians, by group.
%         pr_pairs =    [k x k] matrix  of pairwise probabilities (above diagonal) 
%                         and corresponding sequential Bonferroni significance 
%                         decisions (below diagonal).
%         H_pairs =     corresponding pairwise KW statistic values.
%         hsets =       [nhsets x k] matrix specifying homogeneous subsets.
%

% RE Strauss, 11/7/99
%   12/21/99 - added homogeneous subsets based on pairwise differences.
%   12/25/99 - use observed pr for alpha if pr > 0.05.
%   8/30/00 -  allow for zero iterations.
%   9/12/00 -  added default value for alpha.
%   4/1/01 -   if not randomized, use anova of ranked data for probabilities;
%                call pairwise().

function [pr,H,grp_medians,pr_pairs,H_pairs,hsets] = ...
                                      kruskwal(x,grps,iter,pairwise,hflag,alpha)

  if (nargin < 3) iter = []; end;
  if (nargin < 4) pairwise = []; end;
  if (nargin < 5) hflag = []; end;
  if (nargin < 6) alpha = []; end;

  default_alpha = 0.05;

  if (isempty(iter))
    iter = 0;
  end;
  if (isempty(pairwise))
    pairwise = 0;
  end;
  if (isempty(hflag))
    hflag = 0;
  end;
  if (isempty(alpha))
    alpha = default_alpha;
  end;

  if (min(size(x))>1)
    error('  KRUSKWAL: data must be a vector');
  end;
  if (min(size(grps))>1)
    error('  KRUSKWAL: classification variable must be a vector');
  end;
  
  N = length(x);
  if (length(grps)~=N)
    error('  KRUSKWAL: data and classification vectors not of same length');
  end;

  x = x(:);                               % Convert input to col vectors
  grps = grps(:);

  [grps,x] = sortmat(grps,x);             % Sort data by group
  [grp_id,freq] = uniquef(grps);          % Group identifiers & freqs
  ngrps = length(grp_id);                 % Number of groups
  freq = [0; cumsum(freq)];               % Cumulative frequencies

  if (ngrps<3)
    if (ngrps<2)
      error('  KRUSKWAL: need at least two groups');
    else
      disp('  KRUSKAL warning: no pairwise tests for two groups');
      pairwise = 0;
    end;
  end;

  grp_medians = medians(x,grps);          % Get medians by group
  x = ranks(x);                           % Convert data to ranks
  are_ranks = 1;                          % Flag indicating ranks

  pr = [];
  pr_pairs = [];
  H_pairs = [];
  hsets = [];

  if (pairwise)
    [pr,H,pr_pairs,H_pairs,hsets] = pairwise(x,grps,1,iter,hflag,alpha);
    return;
  end;

  H = kwstat(x,grps,are_ranks);           % KW statistic for all grps

  if (iter)
    incr = 1./iter;
    pr = 0;
    for it = 1:iter                         % Iterate the probabilities
      xr = x(randperm(N));                    % Randomly permute data
      Hr = kwstat(xr,grps,are_ranks);         % Randomized KW statistic for all grps
      if (Hr >= H)
        pr = pr+incr;
      end;
    end;
  else
    x = ranks(x);
    [x,pr] = anova(x,grps);
  end;

  return;
