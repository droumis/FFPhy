% MULTNORMPDF: Multivariate-normal probability density function.
%
%     Usage: fx = multnormpdf(X,mu,C)
%
%         X =  [n x p] matrix of points for which densities are desired.
%         mu = [m x p] matrix of means.
%         C =  [p x p] covariance matrix.
%         --------------------------------------------------------------
%         fx = [n x m] matrix of densities at the values of X, with one 
%                column for each mean.
%

% RE Strauss, 5/7/00

function fx = multnormpdf(X,mu,C)
  [n,p] = size(X);
  [m,mp] = size(mu);
  
  if (size(mu,2)~=p)
    error('  MULTNORMPDF: invalid matrix of means.');
  end;
  if (any(size(C)~=[p p])  | ~iscov(C))
    error('  MULTNORMPDF: invalid covariance matrix.');
  end;

  fx = zeros(n,m);
  den = (2*pi).^(1./(2*p)) * sqrt(det(C));
  invC = inv(C);

  for im = 1:m
    mmu = ones(n,1)*mu(im,:);                  
    d = X-mmu;

    for in = 1:n
      e = -0.5 * d(in,:)*invC*d(in,:)';
      fx(in,im) = (1./den) * exp(e);
    end;
  end;

  if (any(fx(:)>1))
    error('  MULTNORMPDF: densities > 1; check covariance matrix.');
  end;

  return;
