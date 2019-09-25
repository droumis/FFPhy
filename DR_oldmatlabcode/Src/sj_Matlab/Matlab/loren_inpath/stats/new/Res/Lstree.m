% LSTREE: Finds the least-squares additive tree for a distance matrix, using the 
%         algorithm of De Soete (1983).
%
%     Usage: addist = lstree(dist)
%
%         dist = [n x n] symmetric distance matrix among n taxa.
%         ---------------------------------------------------------
%         addist = best-fitting modified distance matrix.
%

% RE Strauss
%   12/28/99 -    changed fminu() to fmins().

function addist = lstree(dist)
  tol = 1e-6;
  if (sum(sum(dist-dist'))>tol | sum(diag(dist))>tol)
    error('LSTREE: input matrix not a distance matrix');
  end;

  delta = trilow(dist);                           % Extract vector of dists
  lend = length(delta);
  delta = zscore(delta)./sqrt(lend-1);            % Normalize SSD==1

  epsvar = var(delta)/3;                          % Variance of epsilon

  d0 = delta + randn(lend,1)*sqrt(epsvar);        % Initial D with noise
  [FDr,r] = lstreef(d0,0,delta);                  % Initial r

  d1 = fmins('lstreef',d0,[],[],r,delta);
  convcrit = sqrt(sum((d1-d0).^2));                 % Convergence criterion
convcrit

  while (convcrit > tol)                          % Iterate till converge
    d0 = d1;                                        % Save last distances
    r = 10*r;                                       % Update constraint weight
    d1 = fmins('lstreef',d0,[],[],r,delta);         % Minimize F(D,r)
    convcrit = sqrt(sum((d1-d0).^2));               % Update convergence criterion
convcrit
  end;

  addist = trisqmat(d1);                          % Return optimal solution

  return;
