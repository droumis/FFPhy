% RANDMISSTEST: Tests the distribution of number of missing values per row or 
%               column of a data matrix against an expected uniform distribution.  
%               Tests observations separately from variables (necessary because 
%               rows and columns of the data matrix are both considered to be 
%               completely permutable).  Uses a randomized goodness-of-fit test 
%               based on the mean squared deviation between observed and expected 
%               values.  The joint probability for rows and columns is based on 
%               the total mean squared deviation for means and columns, under the
%               assumption that rows and columns are independent.
%
%     Usage: [pval,mse] = randmisstest(X,{iter})
%
%         X =       [n x p] data matrix.
%         iter =    optional number of iterations for goodness-of-fit and 
%                     contingency tests [default = 1000].
%         --------------------------------------------------------------------
%         pval =    3-element vector of significance levels for missing values 
%                     per row and column, and joint significance level.
%         mse =     2-element vector of mean squared deviations between 
%                     observed and expected numbers of missing values 
%                     per row and column.
%

% RE Strauss, 3/16/01
%   7/18/01 - added joint probability for rows and columns, under the assumption 
%               that rows and columns are independent.

function [pval,mse] = randmisstest(X,iter)
  if (nargin < 2) iter = []; end;

  if (isempty(iter))
    iter = 1000;
  end;

  [n,p] = size(X);
  X = ~isfinite(X);                         % Replace with boolean matrix                 
  
  colmiss = sum(X);                         % Sums of missing data by columns
  rowmiss = sum(X');                        % Sums of missing data by rows
  N = sum(colmiss);                         % Total missing data

  pval = [1 1 1];
  mse = [0 0];
  if (N<3)                                  % Return if fewer than 3 missing data
    return;
  end;

  rowexp = (N/n)*ones(1,n);                 % Expected uniform row values
  trowhat = mean((rowmiss-rowexp).^2);
  mse(1) = trowhat;

  colexp = (N/p)*ones(1,p);                 % Expected uniform col values
  tcolhat = mean((colmiss-colexp).^2);
  mse(2) = tcolhat;

  trow = zeros(iter,1);
  tcol = zeros(iter,1);

  for it = 1:iter                           % Randomized probabilities
    BX = randbinm(n,p,N);                     % Randomize missing values
    colmiss = sum(BX);
    rowmiss = sum(BX');
    trow(it) = mean((rowmiss-rowexp).^2);
    tcol(it) = mean((colmiss-colexp).^2);
  end;

  pval(1) = randprob(trowhat,trow);         % Significance levels
  pval(2) = randprob(tcolhat,tcol);
  pval(3) = 1-((1-pval(1))*(1-pval(2)));      % Joint

  return;
