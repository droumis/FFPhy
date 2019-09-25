% MISSRANDTEST: Assesses the distribution of number of missing values per row or 
%               column of a data matrix against an expected uniform marginal 
%               distribution.  Returns the relative difference between observed 
%               and expected values for rows, columns, and total (relative to their
%               maximum possible values).  Because the total realtive difference
%               assumes independence of rows and columns, a 
%               contingency msd statistic (comparable to a 2-way contingency 
%               chi-squared statistic) is also returned.  All of these statistics, 
%               which are measures of nonuniformness, must be compared against null 
%               distributions determined by randomization to determine 
%               significance levels.
%               
%
%   Usage:  [msd,stats] = missrandtest(X)
%
%       X =     [n x p] data matrix containing missing values.
%       -----------------------------------------------------------------------
%       msd =   [1 x 4] vector of mean-squared differences between observed and 
%                 expected marginal totals.
%                   1) 
%       stats = [1 x 3] vector of associated statistics for randomization:
%                   number of rows
%                   number of columns
%                   total number of missing values (N)
%                   
%

% RE Strauss, 3/16/01

function [msd,stats] = randmissmeas(X)
  [r,c] = size(X);

  X = ~isfinite(X);                     % Convert to binary matrix

  colobs = sort(sum(X));                     % Marginal totals
  rowobs = sort(rowsum(X));
  N = sum(colobs);

  stats = [N,r,c]
  
  rowexp = sort(binopdf());


  rowexp = (N/r)*ones(1,r);                 % Expected uniform row values
  trowhat = mean((rowmiss-rowexp).^2);
  msd(1) = trowhat;

  colexp = (N/c)*ones(1,c);                 % Expected uniform col values
  tcolhat = mean((colmiss-colexp).^2);
  msd(2) = tcolhat;

  


  return;
