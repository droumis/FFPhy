% ANOVANST: Two-level nested unbalanced random-effects analysis of variance.  
%           Performed separately for each column of the data matrix.
%           Randomization is two-stage:
%             1) random permutation among subgroup, within group  
%                  (for among-subgroup F)
%             2) random permutation among and subgroup and group  
%                  (for among-group F)
%
%   Syntax: [F,pr,df,ss,ms,varcomp,varprop,hsets,nhsets] = 
%                                     anovanst(X,grps,subgrps,{iter},{alpha})
%
%         X =       [n x p] matrix of observations for n observations and p variables.
%         grps =    [n x 1] classification variable.
%         subgrps = [n x 1] classification variable.
%         iter =    optional number of iterations for random-permutation 
%                     estimate of statistical significance 
%                     [default = 0 = asymptotic estimates].
%         alpha =   optional significance level for homogeneous subsets 
%                     [default = 0.05].
%         ----------------------------------------------------------------------------
%         F =       [2 x p] matrix of F-statistics (among-group, among-subgroup).
%         pr =      [2 x p] matrix of corresponding significance levels, 
%                     asymptotic (if iter=0) or randomized (if iter>0).
%         df =      [4 x p] matrix of degrees of freedom (among-group, 
%                     among-subgroup, within-subgroup, total).
%         ss =      [4 x p] matrix of sums-of-squares (among-group,  
%                     among-subgroup, within-subgroup, total).
%         ms =      [3 x p] matrix of mean-squares (among-group, 
%                     among-subgroup, within-subgroup).
%         varcomp = [3 x p] matrix of variance-component estimates (among-group,
%                     among-subgroup, within-subgroup).
%         varprop = [3 x p] matrix of variance components as proportions of
%                     total (among-group, among-subgroup, within-subgroup).
%         hsets =   [nhsets x k] matrix specifying homogeneous subsets for 
%                     k groups; there might be several equally consistent sets for each 
%                     variable.  Subsets for different variables are concatenated.
%         nhsets =  [nhsets x 1] vector indicating the identify of the variable 
%                     to which the corresponding homogeneous subsets are associated.
%

% Sokal & Rohlf, 1981, pp. 294-299.

% RE Strauss, 3/23/97
%   1/22/00 - removed anova calculations as functions 'anovansf' and 'anovanstg'.
%   2/6/00 -  added randomization.
%   7/14/00 - added homogeneous subsets.
%   11/4/00 - corrected problem with unsorted permuted distributions.
%   12/3/00 - use 'homosub' facilities for nested anovas.
%   2/18/02 - added 'nhsets' output vector.

function [F,pr,df,ss,ms,varcomp,varprop,hsets,nhsets] = anovanst(X,grps,subgrps,iter,alpha)
  if (nargin < 4) iter = []; end;
  if (nargin < 5) alpha = []; end;
  
  getvc = 0;
  get_hsets = 0;
  if (nargout >= 6) getvc = 1; end;
  if (nargout >= 8) get_hsets = 1; end;

  [n,p] = size(X);

  if (~isvector(grps) | ~isvector(subgrps))
    error('  ANOVANST: group and subgroup identifiers must be vectors.');
  end;

  if (length(grps)~=n | length(subgrps)~=n)
    error('  ANOVANST: input matrices must have same number of observations.');
  end;

  if (isempty(iter))
    iter = 0;
  end;
  if (isempty(alpha))
    alpha = 0.05;
  end;

  getpr = 0;
  if (~iter)
    getpr = 1;
  end;

  F  = zeros(2,p);
  v1 = zeros(2,p);
  v2 = zeros(2,p);
  pr = zeros(2,p);
  df = zeros(4,p);
  ss = zeros(4,p);
  ms = zeros(3,p);
  varcomp = zeros(3,p);
  varprop = zeros(3,p);

  % Group and subgroup composition

  [subgrps,G,SG,ngrps,n4] = anovanstg(grps,subgrps);

  % ANOVAs

  [F,pr,df,ss,ms,varcomp,varprop] = ...
      anovanstf(X,grps,subgrps,G,SG,n4,getpr,getvc);


  if (iter)                             % Two-stage random permutation of data
    Fw = zeros(iter,p);                   % Within-group, among-subgroup F
    Fa = zeros(iter,p);                   % Among-group F

    [ug,fg] = uniquef(grps);              % Identity and sizes of subgrps

    for it = 1:iter                       % Iterate
      Xp = randpermg(X,grps);               % Permute subgrp data within grps
      Fr = anovanstf(Xp,grps,subgrps,G,SG,n4,0,0);
      Fw(it,:) = Fr(2,:);

      Xp = randpermg(X);                    % Permute all data
      Fr = anovanstf(Xp,grps,subgrps,G,SG,n4,0,0);
      Fa(it,:) = Fr(1,:);
    end;

    pr(1,:) = bootprob(F(1,:),Fa)';       % Find randomized probabilities
    pr(2,:) = bootprob(F(2,:),Fw)';

  end;

  if (get_hsets)                        % Heterogeneous subsets based on nested anovas
    hsets = [];
    nhsets = [];
    for ip = 1:p
      h = homosub(X(:,ip),grps,subgrps,alpha);
      hsets = [hsets; h];
      nhsets = [nhsets; ip*ones(size(h,1),1)];
    end;
  end;

  return;
