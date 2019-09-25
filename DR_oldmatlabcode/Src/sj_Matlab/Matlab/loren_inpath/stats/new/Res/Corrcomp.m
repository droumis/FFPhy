% CORRCOMP: Compare the correlation coefficients for two groups to determine 
%           whether they are significantly different.  The two groups must have 
%           data for the same set of variables, but can differ in sample size.
%           Uses Hotelling's (1953) small-sample modification of Fisher's 
%           z-transform.  Ignores missing data.
%
%     Usage: [pr,signif,R1,R2] = corrcomp(X1,X2,{multcomp},{alpha})
%
%         X1 =        [n x p] data matrix for group 1.
%         X2 =        [m x p] data matrix for group 2.
%         multcomp =  optional boolean value indicating, if true, that judgments 
%                       about statistical significance are to be based on the 
%                       sequential Bonferroni test [default = 0].
%         alpha =     optional alpha level for significance tests 
%                       [default = 0.05].
%         ----------------------------------------------------------------------
%         pr =        [p x p] matrix of probabilities of identity of correlation 
%                       coefficients.
%         signif =    corresponding boolean matrix of significance decisions.
%         R1 =        correlation matrix for group 1.
%         R2 =        correlation matrix for group 2.
%

% RE Strauss, 6/19/01
%   6/20/01 - ignore missing data by deleting observations.

function [pr,signif,R1,R2] = corrcomp(X1,X2,multcomp,alpha)
  if (nargin < 3) multcomp = []; end;
  if (nargin < 4) alpha = []; end;

  if (isempty(multcomp))
    multcomp = 0;
  end;
  if (isempty(alpha))
    alpha = 0.05;
  end;

  i = find(isfinite(rowsum(X1));
  X1 = X1(i,:);
  i = find(isfinite(rowsum(X2));
  X2 = X2(i,:);

  [n1,p] = size(X1);
  [n2,q] = size(X2);
  if (p~=q)
    error('  CORRCOMP: data matrices must have same number of variables.');
  end;

  R1 = corrcoef(X1);
  R2 = corrcoef(X2);

  [Z1,stderr1] = corrz(R1,n1,1);
  [Z2,stderr2] = corrz(R2,n2,1);
  v1 = stderr1*stderr1*n1;
  v2 = stderr2*stderr2*n2;
  pr = zeros(size(Z1));
  df = n1+n2-2;

  for i = 1:p-1
    for j = i+1:p
      m1 = Z1(i,j);
      m2 = Z2(i,j);
      t = abs(m2-m1)./sqrt((((n1-1).*v1+(n2-1).*v2)./(n1+n2-2)).*((n1+n2)./(n1.*n2)));
      pr(i,j) = 1-tcdf(t,df);
    end;
  end;

  pr = pr + pr';
  p = trilow(pr);

  if (multcomp & p>2)
    signif = trisqmat(seqbonf(p,alpha));
    signif = signif + signif';
  else
    signif = (pr <= alpha);
    signif = putdiag(signif,0);
  end;

  return;
