% RankDistMean: Ranks the observations based on absolute distance from the mean, where the
%               observation having rank 1 is closest to the mean.  Optionally ranks within
%               groups.
%
%     Usage: [rankobs,obsrank1] = rankdistmean(X,{g})
%
%         X =       [n x p] data matrix.
%         g =       optional [n x 1] group-membership vector for k groups.
%         --------------------------------------------------------------------
%         rankobs = [n x p] matrix for which each observation is replaced by its rank.
%         obsrank = [k x p] matrix in which the first row comprises the indices (positions) 
%                     of observations of rank 1 for the first group, the second row for group 2,
%                     etc.
%

% RE Strauss, 3/21/03

function [rankobs,obsrank1] = rankdistmean(X,g)
  if (~nargin) help rankdistfrom mean; return; end;
  
  if (nargin < 2) g = []; end;
  
  [N,P] = size(X);
  if (isempty(g)) g = ones(N,1); end;
    
  absdev = abs(grpcentr(X,g));
  rankobs = ranks(absdev,g);
  
  ug = unique(g);
  ngrps = length(ug);
  obsrank1 = zeros(ngrps,P);
  
  for ig = 1:ngrps
    r = rankobs(g==ug(ig),:);
    [rmin,i] = min(r);
    obsrank1(ig,:) = i;
  end;

  return;
  