% CORRPOOL: Calculates pooled within-group covariance matrix, and rescales it to 
%           a correlation matrix.  Missing values are replaced by within-group 
%           means.
%
%     Usage: Cp = corrpool(X,grps,use_ranks)
%
%         X =         [N x P] data matrix.
%         grps =      [N x 1] vector of group identifiers.
%         use_ranks = optional boolean flag indicating that rank correlations 
%                       are to be used.
%         -------------------------------------------------------------------
%         Cp =        pooled within-group correlation matrix.
%

% RE Strauss, 5/26/99
%   11/29/99 - changed calling sequence.

function Cp = corrpool(X,grps,use_ranks)
  if (nargin<3) use_ranks = []; end;

  if (isempty(use_ranks))
    use_ranks = 0;
  end;

  if (use_ranks)                    % Within-grp ranks for Spearman corrs
    X = ranks(X,grps);
  end;

  Cp = covpool(X,grps);             % Pooled within-group covar matrix
  Cp = covcorr(Cp);                 % Convert to correlation matrix

  return;
