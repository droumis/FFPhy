% FACTORP: Estimates the loadings of the primary factor of a covariance matrix,
%          given a corresponding boolean matrix specifying the variables to be
%          included in the secondary factors.  If no secondary-factor matrix is
%          specified, estimates the loadings of the general factor.
%          Also returns the residual covariance matrix.
%          Error handling is done by function WRIGHT.
%
%     Syntax: [f,rc] = factorp(c,{secnd},{primary})
%
%         c =   [n x n] covariance or correlation matrix.
%         secnd = optional [n x s] matrix of boolean values (T/F = 1/0)
%                 indicating by 1's the submatrices of variables (s of them)
%                 specifying the secondary factors.
%         primary = optional boolean flag indicating whether primary (=T)
%                 or general (=F) factor is to be estimated in the presence
%                 of secondary factors [default = T = primary].
%         ---------------------------------------------------------------------
%         f =   [n x 1] vector of factor loadings.
%         rc =  [n x n] matrix of residual covariances.
%

function [f,rc] = factorp(c,secnd,primary)
  TRUE = 1; FALSE = 0;

  if (nargin < 3)
    primary = TRUE;
  end;

  [nvar,p] = size(c);

  skip = eye(nvar);                   % Put 1's on diagonal of 'skip' matrix
  if (nargin>1)                       % If secondary factors specified,
    [n,nsec] = size(secnd);

    if (primary)                      % If primary factor is to be estimated,
      for s = 1:nsec                    % Iterate thru submatrix specifications
        secindx = find(secnd(:,s));         % Convert to vector of indices
        lensec = length(secindx);
        for i = 1:(lensec-1)              % Put 1's in sequestered off-diag submat
          for j = (i+1):lensec
            skip(secindx(i),secindx(j)) = TRUE;
            skip(secindx(j),secindx(i)) = TRUE;
          end;
        end;
      end;
    end;
  end;

  if (sum(sum(c<0)))                  % Any negative covariances?
    neg_covs = TRUE;
  else
    neg_covs = FALSE;
  end;

  if (nvar==2)                        % Exact: two variables
    f = ones(2,1) * sqrt(abs(c(1,2)));
    if (neg_covs)
      f(2) = -f(2);
    end;
  end;

  if (nvar==3 & ~neg_covs)            % Exact: three variables, pos covs
    f = zeros(3,1);
    f(1) = sqrt(c(1,2)*c(1,3)/max(c(2,3),eps));
    f(2) = sqrt(c(1,2)*c(2,3)/max(c(1,3),eps));
    f(3) = sqrt(c(1,3)*c(2,3)/max(c(1,2),eps));
  end;

  if (nvar>3 | neg_covs)              % Estimated solution
    % Estimate loadings from covariance matrix.  Initial guess is that
    % loadings are equal to the root-mean-square covariance
    % (including variance), by column.

%    init_f = sqrt(abs(mean(c)))';
    init_f = min(sqrt(abs(mean(c)))',sqrt(diag(c)));
    f = fmins('factorpf',init_f,[],[],c,skip);
  end;

  rc = c - f*f';                      % Residual covariances

  return;
