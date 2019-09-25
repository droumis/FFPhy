% BIMODAL: Calculates a coefficient of bimodality (original here, not in
%          literature to my knowledge), based on an inverse kurtosis, by column.
%          The index is centered on zero, the asymptotic value for a uniform
%          distribution (no modes).  Positive values are increasingly bimodal;
%          negative values are increasingly unimodal.  Index ranges from 
%          negative to positive infinity.
%
%     Syntax: [coeff,pr,ci,power] = bimodal(X,iter,alpha)
%
%          X =     [n x p] data matrix.
%          iter =  number of iterations for significance and confidence levels.
%          alpha = expected probability of Type I error.
%          --------------------------------------------------------------------
%          coeff = [1 x p] vector of bimodality coefficients.
%          pr =    [1 x p] vector of randomized significance levels, based on a 
%                    null uniform distribution.
%          ci =    [2 x p] matrix of randomized confidence intervals:
%                    row 1 - lower critical values
%                        2 - upper critical values
%          power = [1 x p] vector of power levels.
%

% RE Strauss, 10/21/97

function [coeff,pr,ci,power] = bimodal(X,iter,alpha)
  if (nargin < 2) iter = []; end;
  if (nargin < 3) alpha = []; end;

  if (isempty(iter))
    iter = 0;
  end;
  if (isempty(alpha))
    alpha = 0.05;
  end;

  coeff = bimodalf(X);

  if (iter > 0)
    if (nargout < 4)
      [ci,pr] = bootstrp('bimodalf',[1,2,0],iter,alpha,X,[],0,1);
    else
      [ci,pr,power] = bootstrp('bimodalf',[1,2,1],iter,alpha,X,[],0,1);
    end;
  end;

  return;

