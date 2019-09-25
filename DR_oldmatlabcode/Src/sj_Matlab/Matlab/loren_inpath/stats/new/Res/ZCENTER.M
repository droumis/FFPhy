% ZCENTER: Zero-center a data matrix by column.
%
%     Usage: Z = zcenter(X)
%
%           X = [n x p] data matrix.
%           ---------------------------------
%           Z = [n x p] centered data matrix.
%

% RE Strauss, 6/19/93

function Z = zcenter(X)
   [N,P] = size(X);
   Z = X - ones(N,1)*mean(X);

   return;
