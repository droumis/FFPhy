% ALLOM: Converts column vectors of loadings to allometric coefficients by 
%        scaling to a mean of 1.
%
%     Usage: A = allom(L)
%
%        L = [p x c] matrix of size eigenvector coefficients.
%        -----------------------------------------------------------
%        A = [p x c] matrix of multivariate allometric coefficients.
%

% RE Strauss, 11/21/97

function A = allom(L)
  [P,C] = size(L);
  A = zeros(P,C);
  meanL = mean(L);
  for c = 1:C
    A(:,c) = L(:,c)./meanL(c);
  end;

  return;
