% RANDMVN: Given either vectors of means and variances and a 
%          correlation matrix, or a vector of means and a covariance matrix, 
%          returns a matrix of random multivariate-normal scores.
%
%       Syntax: X = randmvn(N,mu,V,{R},{exact})
%
%           N =  number of MVN observations to be returned.
%           mu = vector (length P) of population means.
%           V =  vector (length P) of population variances, or [P x P] covariance 
%                  matrix.  If V is a scalar, it is used as the variance for all 
%                  variables.
%           R =  optional [P x P] correlation matrix (omitted if V is a covariance 
%                  matrix or if population correlations are all zero).  
%                  If R is a scalar, it is used as the correlation among all 
%                  pairs of variables.
%           exact = optional boolean flag indicating, if true, that covariances 
%                 among random observations are to be exactly those of the target 
%                 matrix.  Default (=0): covariances vary randomly.
%           ----------------------------------------------------------------------
%           X =  [N x P] matrix of random MVN values.
%

% Kaiser, HF & K Dickman.  1962.  Sample and population score matrices and 
%   sample correlation matrices from an arbitrary population correlation matrix.  
%   Psychometrika 27:179-182.

% RE Strauss, 1/2/97
%   8/29/99 -  small changes for Matlab v5.
%   10/5/99 -  convert from Cholesky decomposition to eigen factorization method.
%   11/20/99 - allow for target corrs to be exact.
%   12/8/99 -  added error message for too few observations generated.
%   1/15/00 -  allow for N < P if covars aren't to be exact.
%   5/4/00 -   adjusted for change in isvector(). 
%   5/6/00 -   convert imaginary values to real.
%   5/26/03 -  adjust for return by covcorr() of standard deviations rather than variances.

function X = randmvn(N,mu,V,R,exact)
  if (nargin < 4) R = []; end;
  if (nargin < 5) exact = []; end;

  [isvect,P] = isvector(mu);          % Check for vector mu
  if (~isvect)
    error('  RANDMVN: vector expected for mu');
  end;
  mu = mu(:);                         % Col vector

  if (N<=P & exact)
    disp('  RANDMVN warning: if number of observations is less than number of');
    disp('    variables, covariances will not be exact.');
    exact = 0;
  end;

  if (isempty(R))                     % Correlation matrix not supplied
    [isvect,lenV] = isvector(V);

    if (isvect)                         % Vector of variances
      if (lenV==1)                        % Scalar variance
        V = V * ones(P,1);                % Enlarge into vector
      elseif (lenv~=P)
        error('  RANDMVN: vectors mu and V not of same length');
      end;
      V = V(:);                           % Col vector
      V = sqrt(V);                        % Save as standard deviations
      R = eye(P);                         % Default correlation matrix

    else                                % Covariance matrix
      [r,c] = size(V);                    
      if (r~=c)
        error('  RANDMVN: invalid correlation matrix');
      else
        if (sum(abs(V-V'))>eps)
          error('  RANDMVN: invalid correlation matrix');
        end;
      end;

      [R,S] = covcorr(V);                 % Convert covariances to correlations
      V = S.^2;
    end;
  end;

  if (~isempty(R))                    % Correlation matrix supplied
    [isvect,lenV] = isvector(V);        % Check for covars vector
    if (~isvect & ~isscalar(V))         % Need vector of variances
      error('  RANDMVN: pass covariance or correlation matrix, not both.');
    end;

    if (lenV==1)                        % Scalar variance
      V = V * ones(P,1);                % Enlarge into vector
    end;
    V = V(:);                           % Col vector
    V = sqrt(V);                        % Save variances as standard deviations

    [r,c] = size(R);                    % Check structure of corr matrix
    if (r==1 & c==1)                    % If scalar, expand into matrix
      sr = R;
      R = ones(P,P);
      for i=1:(P-1)
        for j=(i+1):P
          R(i,j) = sr;
          R(j,i) = sr;
        end;
      end;
    end;

    if (~iscorr(R))
      error('  RANDMVN: Invalid correlation matrix.');
    end;
  end;

  % Convert mu & V to matching [N x P] matrices

  mu = ones(N,1) * mu';
  V =  ones(N,1) * V';
  
  % Generate random data matrix

  [E,lambda] = eig(R);                  % Eigen decomposition of corr matrix
  Z = randn(N,P);                       % Random eigenvector scores

  if (exact)                            % If target corrs are to be exact,
    Zcent = zscore(Z);                  %   scores must be standardized and
    [U,M] = eig(Z'*Zcent);              %   uncorrelated
    Z = zscore((sqrt(M)*U'*Zcent')');
  end;

  X = zscore(Z*sqrt(lambda)*E');        % Convert evect scores to observations
  X = X .* V + mu;                      % Scale obs to target means and variances
  X = real(X);                          % Convert any complex values to real

  return;

