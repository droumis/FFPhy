% MANOVAF: Objective function for MANOVA, called by BOOTSTRP.
%          Returns the reciprocal of Wilks' lambda, because the null
%          hypothosis is rejected for small values of lambda.
%
%     Syntax: lambda = manovaf(X,grps,nu1,nu2)
%
%         x =       [n x p] data matrix.
%         grps =    [n x 1] vector of group memberships.
%         nu1,nu2 = arguments passed by bootstrp() but unused.
%         ----------------------------------------------------
%         lambda = reciprocal of Wilks' lambda.
%

% RE Strauss, 11/29/99

function lambda = anovaf(X,grps,nu1,nu2)
  [N,nvars] = size(X);

  G = design(grps);                   % Design matrix
  ngrps = size(G,2);                  % Number of groups

  totmean = ones(N,1)*mean(X);        % Matching matrix of grand means

  grpm = zeros(nvars,ngrps);          % Group means
  for k = 1:ngrps
    grpm(:,k) = mean(X(G(:,k)==1,:))';
  end;
  grpmean = G * grpm';                % Matching matrix of group means

  e = X - grpmean;                    % Within-group deviations
  g = totmean - grpmean;              % Among-group deviations

  sse =  e'*e;                        % Within-group sscp
  ssa =  g'*g;                        % Among-group sscp

  lambda = det(sse)/det(sse+ssa);     % Wilks' lambda
  lambda = 1/lambda;                  % Reciprocal of Wilks' lambda

  return;


