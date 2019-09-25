% << NOT COMPLETE >>
%
% CONTIN3: Three-way contigency-table analysis, based either on fixed
%          marginal totals (type 1, using conventional asymptotic chi-square
%          test) or floating marginal totals (type 2, randomized).
%
%     Syntax: [totprob,totx2,df,x2,prob,exp,pchisq] =
%                 contin3(dimen,obs,{fixed},{iter},{plotchis})
%
%          dimen =    vector specifying dimensions: [r c t].
%          obs =      [r x c x t] matrix of observed counts.  Matrix is
%                       stored as a horizontally concatenated series of
%                       layers of effective dimension [r x c*t], where the
%                       columns for layer t extend from cols [(t-1)*c+1]
%                       thru [t*c].
%          fixed =    optional flag indicating whether marginal (row & column)
%                       totals to be treated as fixed (=TRUE) or
%                       unconstrained (=FALSE).
%          iter =     optional number of randomization iterations [default = 5000].
%          plotchis = optional flag specifying plot of overall randomized
%                       chi-square distribution [default = FALSE].
%
%          totprob =  overall significance-level of observed table.
%          totx2 =    overall observed chi-square value.
%          df =       degrees of freedom for table.
%          x2 =       [r x c*t] matrix of observed chi-square values, by cell.
%          prob =     [r x c*t] matrix of chi-square probabilities, by cell.
%          exp =      [r x c*t] matrix of expected values.
%          pchisq =   asymptotic chi-square probability for type-1 analysis.
%

% Romesberg,HC & K Marshall. 1985. CHITEST: a Monte-Carlo computer program for
%   contingency table tests.  Computers & Geosciences 11:69-78.

function [totprob,totx2,df,x2,prob,exp,pchisq,sumpairx2,probpairx2]...
         = contin3(dimen,obs,fixed,iter,plotchis)
  TRUE = 1; FALSE = 0;

  if (nargin<3)
    fixed = TRUE;
  end;
  if (nargin<4)
    iter = 5000;
  end;
  if (nargin<5)
    plotchis = FALSE;
  end;

  r = dimen(1);
  c = dimen(2);
  t = dimen(3);
  ct = c*t;

  [nr,nct] = size(obs);
  if (nr==1 | nct==1)
    error('  Use 1-way analysis (CONTIN1) for vector of observations');
  end;
  if ((nct ~= ct) | (nr ~= r))
    error('  Incompatible dimensions');
  end;

  rowtot = zeros(r,t);                % Allocate marginal-total matrices
  coltot = zeros(t,c);
  tabtot = zeros(r,c);

  exp =  zeros(r,ct);                 % Allocate working matrices
  prob = zeros(r,ct);

  for it = 1:t                        % Marginal totals
    xtab = obs(:,((it-1)*c+1):(it*c));  % Extract current table
    coltot(it,:) = sum(xtab);           % Column totals
    rowtot(:,it) = sum(xtab')';         % Row totals
    tabtot = tabtot + xtab;             % Table totals
  end;
  tot = sum(sum(obs));                % Grand total

  for it = 1:t                        % Expected cell probabilities
    for ir = 1:r
      for ic = 1:c
        rprob = rowtot(ir,it);
        cprob = coltot(it,ic);
        tprob = tabtot(ir,ic);
        exp(ir,ic+(it-1)*c) = rprob*cprob*tprob;
      end;
    end;
  end;
  exp = exp * tot/sum(sum(exp));      % Convert to expected counts

  x2 = (obs-exp).^2./exp;             % 'Chi-square' statistic per cell
  totx2 = sum(sum(x2));               % Total 'chi-square' statistic
  df = (r-1)*(c-1)*(t-1);             % Degrees of freedom
x2
totx2
df

  if (fixed)                          % Asymptotic chi-square probability
    pchisq = 1-chi2cdf(totx2,df);     %   for fixed marginals only
  else
    pchisq = [];
  end;
pchisq

% ===================

  return;

ncell = r*c;
sumpairx2 = zeros(ncell,ncell);
probpairx2 = zeros(ncell,ncell);
%distpairx2 = zeros(ncell*ncell,iter);
for i = 1:r
  for j = 1:c
    nc1 = (j-1)*r+i;
    for ii = 1:r
      for jj = 1:c
        nc2 = (jj-1)*r+ii;
        sqsumx2 = sqrt(x2(i,j) + x2(ii,jj));
        sumpairx2(nc1,nc2) = sqsumx2;
      end;
    end;
  end;
end;

  if (fixed)                          % Asymptotic chi-square probability
    pchisq = 1-chi2cdf(totx2,df);     %   for fixed marginals only
  else
    pchisq = [];
  end;

  totprob = 0;
  if (plotchis)
    randx2 = zeros(iter,1);
  end;

  for it = 1:iter                     % Iterate random table construction
    obs = continrn(rowtot,coltot,fixed);  % Random table

    x2obs = (obs-exp).^2./exp;        % 'Chi-square' statistic per cell
    totx2obs = sum(sum(x2obs));       % Total 'chi-square' statistic
    if (plotchis)                     % Optionally save for plot of distribution
      randx2(it) = totx2obs;
    end;

for i = 1:r
  for j = 1:c
    nc1 = (j-1)*r+i;
    for ii = 1:r
      for jj = 1:c
        nc2 = (jj-1)*r+ii;
        sqsumx2 = sqrt(x2obs(i,j) + x2obs(ii,jj));
        if (sqsumx2 > sumpairx2(nc1,nc2))
          probpairx2(nc1,nc2) = probpairx2(nc1,nc2) + 1/iter;
        end;
%        distpairx2((nc2-1)*r*c+nc1,it) = sqsumx2;
      end;
    end;
  end;
end;


    [i,j] = find(x2obs > x2);         % Cell probabilities
    leni = length(i);
    if (leni>0)
      for ic = leni
        prob(i(ic),j(ic)) = prob(i(ic),j(ic))+1/iter;
      end;
    end;

    if (totx2obs > totx2)             % Total probability
      totprob = totprob + 1/iter;
    end;
  end;

  if (plotchis)                       % Plot randomized distribution
    histgram(randx2);
    xlabel('Randomized chi-square statistic');
  end;

  return;
