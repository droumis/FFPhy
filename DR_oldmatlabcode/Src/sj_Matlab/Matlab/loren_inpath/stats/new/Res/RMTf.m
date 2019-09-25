% MISSRANDTESTF:  Calculates measures of nonrandomness for missing data within a 
%                 data matrix.
%
%     Usage: [msd,stats] = missrandtest(X)
%
%       X =     [n x p] data matrix containing missing values.
%       -----------------------------------------------------------------------
%       msd =   [1 x 4] vector of relative differences between observed and 
%                 expected marginal totals.
%                   1) relative difference for matrix;
%                   2) relative difference for rows;
%                   3) relative difference for columns;
%                   4) 
%       stats = [1 x 3] vector of associated statistics for randomization:
%                   number of rows
%                   number of columns
%                   total number of missing values (N)
%

function [msd,stats] = missrandtest(X)
  [r,c] = size(X);

%  X = ~isfinite(X);                      % Convert to binary matrix

  colobs = sum(X)                         % Marginal totals
  rowobs = sum(X')
  N = sum(colobs)

  stats = [N,r,c]
  
  rowexp = (N/r)*ones(1,r)               % Expected uniform row values
  rowdev = sum(abs(rowobs-rowexp))        
  maxrowexp = [N zeros(1,r-1)]
  maxrowdev = sum(abs(rowobs-maxrowexp))
  rowstat = rowdev/maxrowdev
  
  colexp = (N/c)*ones(1,c)               % Expected uniform col values
  coldev = sum(abs(colobs-colexp))
  maxcolexp = [N zeros(1,c-1)]
  maxcoldev = sum(abs(colobs-maxcolexp))
  colstat = coldev/maxcoldev
  
  totstat = (rowdev+coldev)/(maxrowdev+maxcoldev)
  
  
  