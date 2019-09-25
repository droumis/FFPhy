% LAWLEY: Objective function for stepdisc().  Given the grouping vector and data 
%         matrix, returns measures of group separation.
%
%     Usage: [V,Vp] = lawley(X,grps)
%
%        X =    [n x p] data matrix (obs x vars).
%        grps = row or column vector of group identifiers.
%        ------------------------------------------------------------------------
%        V =    Rao's V = Lawley-Hotelling trace.  Proportional to the mean 
%                 Mahalanobis distance among groups.
%        Vp =   V divided by number of variables (expressed on per-variable basis).
%

% RE Strauss, 5/10/98
%   11/29/99 - changed calling sequence.

function [V,Vp] = lawley(X,grps)
  [nobs,p] = size(X);

  G = design(grps);                   % ANOVA-type design matrix
ngrps = size(G,2);

  mean_W = (G'*G)\G'*X;               % Within-group means
  devs_W = X - G*mean_W;              % Within-group deviations
  W = devs_W'*devs_W;                 % Within-group SSCP matrix

  devs_T = X - ones(nobs,1)*mean(X);  % Total deviations
  T = devs_T'*devs_T;                 % Total SSCP matrix
  B = T - W;                          % Among-sample SSCP matrix

  rcond_W = rcond(W);

  if (rcond_W < 1e-8)
    V = NaN;
    Vp = NaN;
  else
    invW = inv(W);
    V = trace(B*invW);                % Lawley-Hotelling trace = Rao's V
    Vp = V/p;
  end;

  return;
