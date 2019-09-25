% COMPLEX: Converts a two-column real matrix to complex.
%
%     Usage: Z = complex(X)
%
%        Z = [n x 2] matrix of real values
%        X = [n x 1] vector of complex values
%

% RE Strauss, 6/17/93

function Z = complex(X)
   Z = X(:,1) + X(:,2)*sqrt(-1);
   return;
