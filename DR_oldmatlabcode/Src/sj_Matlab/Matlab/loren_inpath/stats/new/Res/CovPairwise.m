% CovPairwise: Given a data matrix containing missing data, estimates means, variances and
%              covariances using elementwise deletion of missing values for all possible pairs
%              of variables.  Optionally tests whether the covariance matrix is positive
%              definite and, if not, adjusts it to be so.
%
%     Usage: [C,M,N] = covpairwise(X,{makeposdef})
%
%         X =          [n x p] data matrix.
%         makeposdef = optional boolean flag indicating, if true, that covariance matrix is
%                        to be adjusted to be positive definite [default = 0].
%         ---------------------------------------------------------------------------------
%         C =          [p x p] covariance matrix.
%         M =          [1 x p] vector of means.
%         N =          [p x p] matrix of pairwise sample sizes.
%

% RE Strauss, 7/1/02

function [C,M,N] = covpairwise(X,makeposdef)
  if (nargin < 2) makeposdef = []; end;
  
  if (isempty(makeposdef))
    makeposdef = 0;
  end;

  [n,p] = size(X);
  C = zeros(p,p);                     
  M = zeros(1,p);
  N = zeros(p,p);
  nc = 0;
  sumc = 0;

  for i = 1:p                         % For all variables,
    y = X(isfinite(X(:,i)),i);          % Eliminate missing values
    if (~isempty(y))
      M(i) =  mean(y);                    % Estimate means & variances
      C(i,i) = var(y);
      N(i,i) = length(y);
    else
      M(i) = NaN;
      C(i,i) = NaN;
    end;
  end;

  for i = 1:(p-1)                     % For all pairs of variables,
    for j = (i+1):p
      indx = find(isfinite(X(:,i)) & isfinite(X(:,j))); % Find pairwise available values
      if (length(indx)>1)
        y = X(indx,[i j]);
        y(:,1) = y(:,1) - M(i);           % Subtract respective means
        y(:,2) = y(:,2) - M(j);
        c = y'*y / (length(indx)-1);    % Adjusted covariance
        sumc = sumc + c(1,2);
        nc = nc + 1;
        C(i,j) = c(1,2);
        C(j,i) = c(1,2);
        n = size(y,1);
        N(i,j) = n;
        N(j,i) = n;
      end;
    end;
  end;

  [i,j] = find(abs(C)<eps);             % If any covariances are still zero,
  if (~isempty(i))                      %   set to mean covariance
    meanc = sumc / nc;
    for k = 1:length(i)
      C(i(k),j(k)) = meanc;
    end;
  end;

  if (makeposdef)
      [C,isposdef] = posdef(C);           % Make positive definate if not so
  end;
    
  return;
