% FACTORWF: Objective function for FACTORW.  Given the final model (original
%           covariance matrix, primary factor loadings, secondary factor
%           loadings, orthogonality constraints), calculates the sum of
%           squared deviations of observed covariances from predicted
%           (reproduced) covariances.
%

function sse = factorsf(fact,c,secnd,orthog)
  [nvar,nsec] = size(secnd);
  sf = zeros(nvar,nsec);

  pf = fact(1:nvar);                  % Unpack loadings vector into matrices
  low = nvar+1;
  for s = 1:nsec
    i = find(secnd(:,s));
    high = low+length(i)-1;
    sf(i,s) = fact(low:high);
    low = high+1;
  end;

  rc = pf*pf';                        % SSE contributed by reproduced covars
  r = vectcorr(pf,sf);
  r2 = r.^2;                            % Independent var of secondary factors
  for s = 1:nsec
    f = sf(:,s);                        % Secondary factor
    rc = rc + (f*f')*(1-r2(s));         % Reproduced covars
  end;

  diagc = diag(c);                      % Penalize loadings producing negative 
  diagrc = diag(rc);                    %   diagonal residual covariances
  neg_indx = find(diagrc > diagc);
  
  totneg = 0;
  if (length(neg_indx)>0)
    totneg = sum(diagrc(neg_indx) - diagc(neg_indx));
  end;

  if (totneg > 1e-5)
    sse = 1e6 * totneg;
  else
    skip = eye(nvar);                   % Ignore diagonals
    sse = sum((c(~skip)-rc(~skip)).^2); % SSE from predicted covars
%    sse = sse / sum(c(~skip).^2);       %   as proportion of original squared covars

    if (orthog>0)                       % SSE contributed by non-orthogonality
%      r1 = vectcorr(pf,sf);             % of secondary factors with primary factor
      sse = sse + mean(abs(r));            % Sum of absolute vector correlations
    end;

    if (orthog==2 & nsec>1)             % SSE contributed by non-orthogonality
      r2 = vectcorr(sf);                % of secondary factors with one another
      sse = sse + mean(trilow(abs(r2)));   % Sum of absolute vector correlations
    end;
  end;

  return;
