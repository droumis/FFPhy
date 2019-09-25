% CONDFACTOR: Returns a modified version of the matrix-condition factor, given the 
%             corr/cov matrix.  Returns NAN if the input corr/cov matrix is not and 
%             cannot be made positive-definite.
%
%     Usage: [cf,evals] = condfactor(C,{condtype})
%
%         C =         correlation or covariance matrix.
%         condtype =  optional type of condition factor to be used:
%                       1 = condition factor, which measures stability of the 
%                           eigensolution, and is greatest when the ratio 
%                           (eigenvalue 1)/(eigenvalue 2) is large [default];
%                       0 = reciprocal condition factor, which measures 
%                           stability under matrix inversion, and is greatest 
%                           when the ratio (eigenvalue 1)/(eigenvalue 2) is  
%                           small.
%         -------------------------------------------------------------------
%         cf =        condition factor.
%         evals =     sorted eigenvalues.
%

% RE Strauss, 6/21/02
%   7/1/02 -  fix problem with sorting eigenvalues from large to small.
%   9/23/02 - return eigenvalues.

function [cf,evals] = condfactor(C,condtype)
  if (nargin < 2) condtype = []; end;
  
  if (isempty(condtype))
    condtype = 1;
  end;

  if (condtype~=0 & condtype~=1)
    error('  CONDFACT: invalid type of condition factor.');
  end;

  Ciscorr = 0;                            % Check if C is a correlation matrix
  if (iscorr(C))
    Ciscorr = 1;
  end;
  if (~Ciscorr)
    if (~iscov(C))
      error('  CONDFACT: invalid input correlation/covariance matrix.');
    end;
  end;

  [C,isposdf] = posdef(C);               % Make positive definate if not so
  if (Ciscorr)                            % Optionally rescale to corr matrix              
    C = covcorr(C);
  end;

  if (isposdf >= 0)
    evals = log(-sort(-eig(C)));  
  else
    evals = [NaN NaN];
  end;
  
  if (condtype)       
%    condfactorig = log(cond(C));         % Condition factor
    cf = evals(1)-evals(2);               % Modified: log(e1)-log(e2) = log(e1/e2)
  else
%    condfactorig = log(rcond(C));        % Reciprocal condition factor
    cf = evals(2)-evals(1);               % Modified: log(e2)-log(e1) = log(e2/e1)
  end;

  return;
  

