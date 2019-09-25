% GENEIGEN: Performs a generalized eigen analysis for two square matrices,
%           and returns the sorted eigenvectors and eigenvalues.
%
%     Syntax: [evects,evals] = geneigen(A,B)
%
%        A,B =    square matrices
%        ------------------------------------------------------
%        evects = matrix of eigenvectors (columns)
%        evals =  vector of eigenvalues in descending sequence
%

function [evects,evals] = geneigen(A,B)
  [evects,D] = eig(A,B);
  evects = real(evects);
  D = real(D);
  [evals,k] = sort(diag(-D));
  evals = -evals;
  evects = evects(:,k);

  return;

