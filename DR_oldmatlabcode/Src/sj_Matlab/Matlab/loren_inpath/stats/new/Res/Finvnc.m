% FINVNC: Returns the inverse of the cumulative distribution function for the
%         noncentral F distribution.
%
%     Syntax: f = finvnc(prob,v1,v2,lambda)
%
%         prob =   target cumulative probability.
%         v1,v2 =  degrees of freedom.
%         lambda = noncentrality parameter.
%         ---------------------------------------
%         f =      F, given prob.
%

% Note: this routine is very inefficient because it doesn't compute the
% inverse directly, but homes in on it iteratively by function minimization.
% Needs to be reprogrammed, but I can't come up with an algorithm.

% RE Strauss, 4/7/96

function f = finvnc(prob,v1,v2,lambda)
  tol = 1e-6;
  incr = 0.05;

  if (prob<0 | prob>1-tol)
    error('  FINVNC: Probability out of range');
  end;

%lambda
%prob

  fmean = v2.*(v1+lambda)./((v2-2).*v1);
  if (prob<0.5)
    maxf = fmean;
  else
    maxf = fmean*10;
  end;
%fmean
%maxf

%  x = fmin('finvncob',0,1-tol,[],prob,v1/2,v2/2,lambda);
  f = fmin('finvncob',0,maxf,[],prob,v1,v2,lambda);
%  f = x*v2/(v1-x*v1);

%prf = fdistnc(f,v1,v2,lambda)

  return;
