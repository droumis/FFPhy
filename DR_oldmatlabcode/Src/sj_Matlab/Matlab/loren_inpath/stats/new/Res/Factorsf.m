% FACTORSF: Objective function for FACTORS.  Given a primary factor and a
%           set of secondary factors, calculates the sum of squared deviations
%           of observed covariances from predicted (reproduced) covariances.
%

function sse = factorsf(f,c,pfactor,secnd,orthog)
  [nvar,nsec] = size(secnd);

  fact = zeros(nvar,nsec);            % Unpack loadings vector into matrix
  low = 1;
  for s = 1:nsec
    secvect = find(secnd(:,s));
    high = low+length(secvect)-1;
    fact(secvect,s) = f(low:high);
    low = high+1;
  end;

  skip = eye(nvar);                   % Ignore diagonals
  rc = zeros(nvar,nvar);
  r = vectcorr(pfactor,fact);
  r2 = r.^2;                          % Independent variances of secondary factors

  for s = 1:nsec                      % SSE contributed by reproduced covars
    ff = fact(:,s);                     % Secondary factor
%    rc = rc + (ff*ff')*(1-r2(s));       % Reproduced covars
    rc = rc + ff*ff';
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
    sse = sum((c(~skip)-rc(~skip)).^2); % SSE from predicted covars
%    sse = sse / sum(c(~skip).^2);       %   as proportion of original squared covars

    if (orthog>0)                       % SSE contributed by non-orthogonality
%      r1 = vectcorr(pfactor,fact);      % of secondary factors with primary factor
      sse = sse + mean(abs(r));          % Mean absolute vector correlations
    end;

    if (orthog==2 & nsec>1)             % SSE contributed by non-orthogonality
      r2 = vectcorr(fact);              % of secondary factors with one another
      sse = sse + mean(trilow(abs(r2)));  % Mean absolute vector correlations
      end;
    end;
  end;

  return;
