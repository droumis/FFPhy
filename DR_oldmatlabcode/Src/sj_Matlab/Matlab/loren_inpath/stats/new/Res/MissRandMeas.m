% MISSRANDMEAS: Assesses the distribution of number of missing values per row or 
%               column of a data matrix against an expected uniform marginal 
%               distribution.  Returns the relative difference between observed 
%               and expected values for rows, columns, and total (relative to their
%               maximum possible values when all missing values are concentrated
%               within a single row or column).  Because the total relative difference
%               assumes independence of rows and columns, a 
%               contingency MAD statistic (comparable to a 2-way contingency 
%               chi-squared statistic) is also returned.  All of these statistics, 
%               which are measures of nonuniformity, must be compared against null 
%               distributions determined by randomization to determine 
%               significance levels.  If the matrix contains fewer than 2 missing
%               values, NaN's are returned for all measures of nonuniformity.
%               
%   Usage:  [nonuni,stats] = missrandmeas(X)
%
%       X =     [n x p] data matrix containing missing values.
%       -----------------------------------------------------------------------
%       nonuni =   [1 x 4] vector of mean-squared differences between observed and 
%                 expected marginal totals.
%                   1) relative difference for matrix;
%                   2) relative difference for rows;
%                   3) relative difference for columns;
%                   4) contingency MAD statistic;
%       stats = [1 x 4] vector of associated statistics for randomization:
%                   1) number of rows
%                   2) number of columns
%                   3) total number of missing values (N)
%                   4) proportion of missing values in matrix
%                   
%

% RE Strauss, 8/24/01

function [nonuni,stats] = missrandmeas(X)
  [r,c] = size(X);
  
  X = ~isfinite(X);                         % Convert to binary matrix
  if (sum(sum(X))>1)                        % If any missing data,
    [nonuni,stats] = missrandmeasf(X);      %   get stats
  else
    [r,c] = size(X);
    colobs = sum(X);                        % Marginal totals
    rowobs = sum(X');
    N = sum(colobs);
    p = N/(r*c);
    stats = [r,c,N,p];                      % Missing-data counts
    nonuni = [NaN NaN NaN NaN];
  end;
  
  return;
