% GRPCENTR: Zero-center the columns of a data matrix by group
%
%     Usage: Z = grpcentr(X,grps)
%
%         X =    [n x p] data matrix.
%         grps = group-identification vector of length n.
%         -----------------------------------------------
%         Z =    [n x p] group-centered data matrix.
%

% RE Strauss, 6/19/93
%   11/29/99 - changed calling sequence.

function Z = grpcentr(X,grps)
   G = design(grps);             % ANOVA-type design matrix
   mean_W = (G'*G)\G'*X;         % Within-group means
   Z = X - G*mean_W;             % Group-center all variables

   return;


