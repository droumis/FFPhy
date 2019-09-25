% RANDWALK:  Simulates a random walk (Brownian motion) in p-space with a 
%            specific degree of correlation among the random deviations.
%
%     Syntax:  [eigvals,A,V,nu,R,crds] = randwalk(N,P,devcorr)
%
%             N -       number of points along path (path length).
%             P -       dimension of space.
%             devcorr - correlation among random deviations [default = 0].
%             ------------------------------------------------------------
%             eigvals - [P-1] vector of sorted eigenvalues of coordinates, 
%                         scaled to min(eigvals)=1 (for P>1 only).
%             A -       asphericity (for P>1 only).
%             V -       ellipticity (for P>1 only.
%             nu -      slope of log(R(t)) vs t.
%             R -       [N x 1] vector of distances of path nodes from origin.
%             crds -    [N x P] matrix of path coordinates.
%

% RE Strauss, 12/10/96

function [eigvals,A,V,nu,R,crds] = randwalk(N,P,devcorr)
  if (nargin < 3) devcorr = []; end;

  if (isempty(devcorr))
    devcorr = 0;
  end;

  get_eig = 0;                            % Flags for output arguments
  get_R = 0;
  get_nu = 0;

  if (nargout >= 1 & P > 1)
    get_eig = 1;
  end;
  if (nargout >= 4)
    get_nu = 1;
  end;
  if (nargout >= 4 | get_nu)
    get_R = 1;
  end;

%  if (devcorr < eps)
    crds = cumsum(randn(N,P));
%  end;

  if (get_eig)
    eigvals = -sort(-eig(corrcoef(crds)));

    A = 0;
    for i = 1:(P-1)
      for j = (i+1):P
        A = A + (eigvals(i)-eigvals(j)).^2;
      end;
    end;
    A = A ./ ((P-1) * eigvals'*eigvals);

    V = prod(eigvals.^(1/P)) ./ mean(eigvals);

    eigvals = eigvals(1:P-1) ./ eigvals(P);
  end;

  if (get_R)
    R = eucl(crds,zeros(1,P));
  end;

  if (get_nu)
    t = log([1:N]);
    sst = t*t';
    ssRt = t*log(R);
    nu = ssRt ./ sst;
  end;

  return;
