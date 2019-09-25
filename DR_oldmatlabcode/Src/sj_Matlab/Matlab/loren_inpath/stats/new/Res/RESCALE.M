% RESCALE:  Given a group-identification vector and data matrix, rescales the
%           columns of the matrix by the pooled within-group stdev.
%
%     Usage: Y = rescale(X,grps)
%
%         X =    [N x P] data matrix.
%         grps = [N x 1] vector of group identifications.
%         -----------------------------------------------------------
%         Y =    [N x P] data matrix centered and rescaled by column.
%

% RE Strauss, 6/27/93
%   11/29/99 - reversed calling sequence

function Y = rescale(X,grps)
  G = design(grps);                 % Design matrix
  [N,ngrps] = size(G);
  [N,nvars] = size(X);

  mean_W = (G'*G)\G'*X;             % Within-group means
  Cp = zeros(nvars,nvars);          % Initialize pooled covar matrix

  for g = 1:ngrps
    index = find(G(:,g));             % Indices to nonzero elements
    n = length(index);                % Sample size of current group
    Y = X(index,:);                   % Get data for current group
    Ydev = Y - ones(n,1)*mean(Y);     % Deviations
    C = (1/(n-1))*Ydev'*Ydev;         % Within-group covariance matrix
    Cp = Cp + ((n-1) * C);            % Augment to pooled matrix
  end;
  Cp = Cp/(N-ngrps);                % Final pooled covariance matrix

  Y = zcenter(X) * inv(sqrt(diag(diag(Cp))));  % Recenter and rescale

  return;
