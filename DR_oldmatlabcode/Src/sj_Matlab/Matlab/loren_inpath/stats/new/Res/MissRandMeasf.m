% MISSRANDMEASF:  Calculates measures of nonrandomness for missing data within a 
%                 data matrix.  Relative difference is the sum of absolute deviations
%                 of row or column totals from an expected uniform distribution,
%                 relative to the maximum possible when all missing values are
%                 concentrated within a single row or column.
%
%     Usage: [nonuni,stats] = missrandmeasf(X)
%
%       X =       [n x p] data matrix containing missing values.
%       -----------------------------------------------------------------------
%       nonuni =  [1 x 4] vector of relative differences between observed and 
%                 expected marginal totals.
%                   1) relative difference for matrix;
%                   2) relative difference for rows;
%                   3) relative difference for columns;
%                   4) contingency MAD statistic;
%       stats =   [1 x 4] vector of associated statistics for randomization:
%                   1) number of rows;
%                   2) number of columns;
%                   3) total number of missing values (N);
%                   4) proportion of missing values in matrix;
%

% RE Strauss, 8/24/01
%   10/30/02 - allow for all missing values to be in same row or column.

function [nonuni,stats] = missrandmeasf(X)
  [r,c] = size(X);

  colobs = -sort(-sum(X));                  % Marginal totals
  rowobs = -sort(-sum(X'));
  N = sum(colobs);

  p = N/(r*c);
  stats = [r,c,N,p];                        % Missing-data counts
  
  rowexp = (1/r)*ones(1,r);                 % Expected uniform row values
  rowexp = -sort(-prbcount(rowexp,N,[],1,1));
  rowdev = sum(abs(rowobs-rowexp));
  maxrowexp = [N zeros(1,r-1)];
  maxrowdev = sum(abs(rowobs-maxrowexp));
  if (maxrowdev==0)
    maxrowdev = 2*sum(maxrowexp);
  end;
  rowstat = rowdev/maxrowdev;
  
  colexp = (1/c)*ones(1,c);                 % Expected uniform col values
  colexp = -sort(-prbcount(colexp,N,[],1,1));
  coldev = sum(abs(colobs-colexp));
  maxcolexp = [N zeros(1,c-1)];
  maxcoldev = sum(abs(colobs-maxcolexp));
  if (maxcoldev==0)
    maxcoldev = 2*sum(maxcolexp);
  end;
  colstat = coldev/maxcoldev;
  
  totstat = (rowdev+coldev)/(maxrowdev+maxcoldev);
  
  cellexp = rowexp'*colexp/N;
  contstat = sum(sum(abs(X-cellexp)))/(r*c);
  
  nonuni = [rowstat colstat totstat contstat];
  
  return;
  
  
  
  