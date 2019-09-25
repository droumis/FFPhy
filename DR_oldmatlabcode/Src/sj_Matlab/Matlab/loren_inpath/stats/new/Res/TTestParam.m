% TTestParam:  Two-sample t-test, given only the parameters (means, standard deviations, 
%       and sample sizes) for the two groups.  For >2 groups, calculates all possible 
%       pairwise among k groups.
%
%     Usage: [pr,t,df,bonf] = ttestparam(m,s,n,{tail},{alpha})
%
%       m =     k-element vector of sample means.
%       s =     k-element vector of standard deviations.
%       n =     k-element vector of sample sizes.
%       tail =  optional indicator of direction of test:
%                 -1 = one-tailed test of mean1 <  mean2
%                  0 = two-tailed test of mean1 ~= mean2 [default]
%                 +1 = one-tailed test of mean1 >  mean2
%       alpha = optional overall alpha level, for k>2 groups [default=0.05],
%                 for the sequential Bonferroni test.
%       ---------------------------------------------------------------------
%       pr =    significance level.
%       t =     t-statistic value.
%       df =    degrees of freedom.
%       bonf =  sequential Bonferroni significance decisions at overall alpha 
%                 level, for k>2 groups.
%

% RE Strauss, 5/14/99
%   10/15/02 - renamed from tt2().
%   12/5/02 -  improved the documentation.

function [pr,t,df,bonf] = tt2(m,s,n,tail,alpha)
  if (nargin < 4) tail = []; end;
  if (nargin < 5) alpha = []; end;

  get_bonf = 0;
  if (nargout>3)
    get_bonf = 1;
  end;

  if (isempty(tail))
    tail = 0;
  end;
  if (isempty(alpha))
    alpha = 0.05;
  end;

  if (min(size(m))>1 | min(size(s))>1 | min(size(n))>1)
    error('  TTestParam: input must be vectors');
  end;

  k = length(m);
  if (length(s)~=k | length(n)~=k)
    error('  TTestParam: input must be vectors of identical length');
  end;

  if (alpha>1)
    alpha = alpha/100;
  end;
  bonf = [];

  if (k==2)                         % 2 groups
    df = sum(n)-2;
    s2p = ((n(1)-1)*s(1)^2+(n(2)-1)*s(2)^2)/df;
    nn = sum(n)/prod(n);
    t = (m(1)-m(2))/sqrt(s2p*nn);

    if (tail>0)
      pr = 1-tcdf(t,df);              % Right-tailed test
    elseif (tail<0)
      pr = tcdf(t,df);                % Left-tailed test
    else
      pr = 2*tcdf(-abs(t),df);        % Two-tailed test
    end;
  
  else                              % k groups
    df = zeros(k,k);
    t =  zeros(k,k);
    pr = zeros(k,k);
    bonf = zeros(k,k);

    for i = 1:(k-1)
      for j = (i+1):k
        nn = n([i,j]);
        df(i,j) = sum(nn)-2;
        s2p = ((n(i)-1)*s(i)^2+(n(j)-1)*s(j)^2)/df(i,j);
        nn = sum(nn)/prod(nn);
        t(i,j) = (m(i)-m(j))/sqrt(s2p*nn);

        if (tail>0)
          pr(i,j) = 1-tcdf(t(i,j),df(i,j));        % Right-tailed test
        elseif (tail<0)
          pr(i,j) = tcdf(t(i,j),df(i,j));          % Left-tailed test
        else
          pr(i,j) = 2*tcdf(-abs(t(i,j)),df(i,j));  % Two-tailed test
        end;

        df(j,i) = df(i,j);
        t(j,i) =  t(i,j);
        pr(j,i) = pr(i,j);
      end;
    end;

    if (get_bonf)
      p = trilow(pr);
      s = seqbonf(p,alpha);
      bonf = trisqmat(s);
    end;
  end;

  return;