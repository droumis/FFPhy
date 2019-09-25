% MOSIMANN: Estimates the Mosimann size vector, which is independent of all 
%           possible log(ratios) of original data.
%
%   Syntax: [S,scrs] = mossiman(X)
%           
%           X =    [n x p] matrix of log-transformed data
%           =============================================
%           S =    [p x 1] vector of coefficients
%           scrs = [n x 1] vector of size scores
%

% Mosimann,JE & FC James. 1979. Evolution 33:444-459.
% Bookstein et al. 1985. Red book, pp. 27-32.

function [S,scrs] = mosimann(X)
  C = cov(X);
  [nvars,ncols] = size(C);

  invC = inv(C);
  const = 1/(sum(sum(invC)));

  S = (const .* (ones(1,nvars)*invC))';

  if (nargout>1)
    scrs = X*S;
  end;

  return;

