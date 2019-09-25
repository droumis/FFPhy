% CISignif: Determines the pairwise levels of statistical significance from data 
%           on means and confidence intervals only, based on two-tailed t-tests 
%           (if sample sizes are provided) or z-tests, assuming that the sample 
%           statistics being compared are normally distributed (or 
%           t-distributed, for finite sample sizes).
%
%     Usage: [pr,signif] = cisignif(statmean,statci,{sampsize},{ci_level},{alpha})
%
%         statmean =  [n x 1] vector of means
%         statci =    [n x 1] vector of half-intervals (+/- values) of confidence 
%                       intervals (if symmetric), or [n x 2] matrix of lower and 
%                       upper bounds (if asymmetric).
%         sampsize =  optional [n x 1] vector of corresponding sample sizes.  If 
%                       provided, pairwise t-tests are used; if not, pairwise 
%                       z-tests are used.
%         ci_level =  optional level of percent-confidence for intervals (assumed 
%                       identical for all distributions being tested 
%                       [default = 95].
%         alpha =     optional alpha-level for decisions of significance of 
%                       pairwise probabilities based on sequential Bonferroni 
%                       test.
%         -----------------------------------------------------------------------
%         pr =        [n x n] square symmetric matrix of pairwise probabilities 
%                       of identity.
%         signif =    corresponding square symmetric boolean matrix of 
%                       significance decisions at an overall level of alpha 
%                       [default = 0.05].
%

% RE Strauss, 6/4/01

function [pr,signif] = cisignif(statmean,statci,sampsize,ci_level,alpha)
  if (nargin < 3) sampsize = []; end;
  if (nargin < 4) ci_level = []; end;
  if (nargin < 5) alpha = []; end;

  get_signif = 0;
  if (nargout > 1)
    get_signif = 1;
  end;

  if (isempty(ci_level))
    ci_level = 0.95;
  end;
  if (isempty(alpha))
    alpha = 0.05;
  end;

  if (ci_level > 1)
    ci_level = ci_level / 100;
  end;

  ci_level = ci_level(1);

  statmean = statmean(:);
  n = length(statmean);

  if (size(statci,2) == 2)
    statci = rowmean(abs(statci-statmean*ones(1,2)));
  else
    statci = statci(:);
    if (length(statci) ~= n)
      error('  CISIGNIF: Incompatible mean and CI matrices.');
    end;
  end;

  if (~isempty(sampsize))
    sampsize = sampsize(:);
    if (length(sampsize) ~= n)
      error('  CISIGNIF: Incompatible mean and sample-size matrices.');
    end;
  else
    sampsize = 1e6*ones(size(statmean));      % N's large enough for z-test
  end;

  pr = ones(n,n);                             % Allocate output matrix
               
  s = statci ./ -tinv((1-ci_level)/2,sampsize); % Convert CI's to standard errors

  for i = 1:(n-1)
    for j = (i+1):n
      m1 = statmean(i);
      m2 = statmean(j);
      v1 = s(i)^2;
      v2 = s(j)^2;
      n1 = sampsize(i);
      n2 = sampsize(j);

      d = ((n1-1)*v1+(n2-1)*v2)/(n1+n2-2);
      t = (m2-m1)/sqrt(d);

      df = n1+n2-2;
      pval = tcdf(-abs(t),df);
      
      pr(i,j) = 2*pval;
      pr(j,i) = pr(i,j);
    end;
  end;
  if (get_signif)
    signif = seqbonf(pr,alpha);
  end;

  return;
