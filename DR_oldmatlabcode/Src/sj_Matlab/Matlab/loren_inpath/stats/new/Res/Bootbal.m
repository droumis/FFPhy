% BOOTBAL: Creates a matrix of randomized observation-indices for a balanced
%          bootstrap.
%
%     Syntax: indx = bootbal(nobs,iter)
%
%          nobs = number of observations
%          iter = number of desired bootstrap iterations
%          -----------------------------------------------------------------
%          indx = [iter x nobs] matrix of randomized indices, such that each
%                   observation is equally represented across the matrix
%

% Note: This procedure is rather inefficient, but is done 'in place' using a
%   double-permutation scheme to avoid memory problems for very large numbers
%   of observations and iterations.

% RE Strauss, 3/23/96

function indx = bootbal(nobs,iter)
  indx = zeros(iter,nobs);              % Allocate index matrix

  for r = 1:iter                        % Randomly permute obs for each iteration
    indx(r,:) = randperm(nobs);
  end;

  for c = 1:nobs                        % Randomly permute across iterations
    i = randperm(iter);
    indx(:,c) = indx(i,c);
  end;

  return;
