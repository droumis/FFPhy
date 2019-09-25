% POSDEF: If a correlation or covariance matrix is positive definite, returns 
%         the same matrix and a flag value of 1.  If not, returns a corresponding 
%         matrix that is positive definite by rescaling negative eigenvalues to a
%         value of 0, holding constant the total variance.  Finds the corresponding 
%         matrix by adjusting to be positive and then reconstructing matrix.  If the 
%         input matrix is not square-symmetric, returns a flag of 0 and a null matrix.
%
%     Usage: [Cpd,ispdef] = posdef(C)
%
%         C =      correlation or covariance matrix.
%         ------------------------------------------------------------------------
%         Cpd =    corresponding positive-definite correlation or covariance 
%                       matrix.
%         ispdef = boolean flag: 1 = input matrix is positive definite
%                                0 = input matrix is not positive definite
%                               -1 = input matrix is not positive definite but
%                                      cannot be adjusted to positive definite.
%

% RE Strauss, 11/20/99
%   11/2/01 - use issqsym() to determine whether is square-symmetric.
%   5/14/02 - replace nonpos evals with 10*eps (applied repeatedly) rather than eps.
%   5/17/02 - use while loop to ensure that all final eigenvalues are positive.
%   5/30/02 - add condition for (ispdef == -1)

function [C,ispdef] = posdef(C)
  tol = eps;

  ispdef = 0;
  if (~issqsym(C))                        % Check whether is square-symmetric
    C = [];
    return;
  end;

  [evect,eval] = eig(C);                  % Eigen decomposition
  eval = diag(eval);

  if (all(eval>0))                        % Matrix is pos def
    ispdef = 1;
  else                                    % Otherwise,
    e = eval;
% e    
    ineg = find(eval<eps);                  % Find nonpositive eigenvalues 
    eval(ineg) = tol*ones(length(ineg),1);  % Substitute tol for nonpos evals
    diagC = diag(C);
    totvar = sum(diagC);                    % Total variance
    Csave = C;

    while (any(e<eps) | any(~isreal(e)))    % Until all eigenvalues become non-negative and real
      eval(ineg) = 10*eval(ineg);             % Increase value for nonpos evals
      eval = eval*totvar/sum(eval);           % Adjust to maintain total variance

      C = evect*diag(eval)*evect';            % Reconstruct matrix

      tri = (tril(C)+triu(C)')./2;            % Adjust round-off asymmetrics
      C = trisqmat(trilow(tri),diag(C));
      e = eig(C);
% C
% eval
% e
      if (any(sum(C)<eps) | sum(sum(abs(C-Csave))/prod(size(C)))<1e-6) % If C hasn't changed, exit
        ispdef = -1;
        e = 1;
      end;
      Csave = C;
    end;
    R = covcorr(C);                         % Rescale to original diagonal elements
    C = covcorr(R,sqrt(diagC));
  end;

  return;
