% DISTDIFF: Given a pairwise distance or similarity matrix, and a boolean matrix 
%           indicating whether pairs of units are within the same group or 
%           in different groups, performs three tests for significant 
%           differences within- versus between-group: Mantel's test, t-test, and 
%           Mann-Whitney rank-sum test.  The latter two tests do not allow for 
%           the non-independence of pairwise measures.
%
%     Usage: [mantel_pr,ttest_pr,mannwhit_pr,n0,n1] = distdiff(distmat,grpmat)
%         
%           distmat =     [n x n] square symmetric matrix of pairwise distances
%                           or similarities.
%           grpmat  =     [n x n] corresponding boolean matrix indicating whether 
%                           pairs are in the same (=0) or different (=1) groups.
%           -------------------------------------------------------------------
%           mantel_pr =   probability for Mantel's test.
%           ttest_pr  =   probability for t-test.
%           mannwhit_pr = probability for Mann-Whitney test.
%           n0 =          number of pairs in same group.
%           n1 =          number of pairs in different groups.
%

% RE Strauss, 6/30/99

function [mantel_pr,ttest_pr,mannwhit_pr,n0,n1] = distdiff(distmat,grpmat)
  dm = trilow(distmat);               % Extract lower triangular matrices
  gm = trilow(grpmat);

  n0 = sum(gm==0);
  n1 = sum(gm==1);

  [r,Z,mantel_pr] = mantel(distmat,grpmat);
  [t,ttest_pr] = ttest(dm,gm);
  [mannwhit_pr,U,N,W] = mannwhit(dm,gm);

  return;