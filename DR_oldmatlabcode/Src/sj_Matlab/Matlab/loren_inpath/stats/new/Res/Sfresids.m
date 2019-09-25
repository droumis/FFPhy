% SFRESIDS: Size-invariant residuals.
%
%     Usage: [R,wload,aload,percvar] = sfresids(X,grps,{kind})
%
%        X =        [n x p] data matrix (obs x vars).
%        grps =     row or column vector of group identifiers.
%        kind =     kind of size vector: 'w' for within-group (default)
%                     or 'a' for among-group.
%        -----------------------------------------------------------------
%        R =        [n x p] matrix of size-invariant residuals.
%        wload =    [p x 1] vector of within-grp size-vector loadings (as 
%                     vector correlations).
%        aload =    [p x 1] vector of among-grp size-vector loadings.
%        percvar =  percent total variance for within-grp size vector.
%

% RE Strauss, 5/16/96
%   11/29/99 - changed calling sequence.

function [R,wload,aload,percvar] = sfresids(X,grps,kind)
   [nobs,nvars] = size(X);       % Numbers of observations & variables
   index = uniquef(grps);
   ngrps = length(index);        % Number of groups

   default_kind = 'w';           % Default = within

   within = 1;
   if (nargin < 4)
      kind = default_kind;
   elseif (kind ~= 'w')
      within = 0;
   end;

   % Size-invariant residuals from pooled within-group size vector

   if (within)
      Z = grpcentr(X,grps);      % Group-center all variables
   else
      Z = X;                     %   or not
   end;
   [wload,percvar,wsize] = pcacov(Z,1);  % Within-group size scores
   [aload,p,S] = pcacov(X,1);    % Among-group size scores
   B = linregr(wsize,Z);         % Regr slopes on within-grp size vector
   pred = S * B(2,:);            % Pred values based on among-grp slopes
   R = zcenter(X-pred);          % Centered size-invariant residuals

   return;

