% RatioProb: Tests whether an observed sex ratio (M/F) differs from unity.
%            Given the number of males and females, estimates the significance level
%            for the test of H0:M/F=1, and determines the minimum sample size needed
%            to detect a difference from unity.  The test may be left-tailed (for an
%            expected female bias), right-tailed (for a male bias), or 2-tailed (for
%            a difference).
%
%     Usage: [pr,signif,minN,ratio] = ratioprob(M,F,{tail},{alpha})
%
%         M =    vector (length k, for k tests) of male counts.
%         F =    corresponding vector of female counts.
%         tail = optional flag (scalar or vector) indicating directions of tests:
%                   -1 = left-tailed
%                    0 = 2-tailed [default]
%                   +1 = right-tailed
%         alpha = optional alpha level for minimum sample sizes and sequential
%                   Bonferroni adjustment (family-wide) [default = 0.05].
%         -----------------------------------------------------------------------
%         pr =   significance levels [right tail of the binomial B(N,0.5)].
%         signif = corresponding significance decisions based on sequential
%                   Bonferroni adjustment.
%         minN = minimum sample sizes (M+F) needed to significantly detect 
%                  differences of the observed magnitudes.
%         ratio = realized sex ratios (M/F).
%

% RE Strauss, 2/27/02

function [pr,signif,minN,ratio] = ratioprob(M,F,tail,alpha)
  if (nargin < 3) tail = []; end;
  if (nargin < 4) alpha = []; end;
  
  get_minN = 0;
  if (nargout > 2)
    get_minN = 1;
  end;

  if (isempty(tail))
    tail = 0;
  end;
  if (isempty(alpha))
    alpha = 0.05;
  end;

  M = M(:);
  F = F(:);
  tail = tail(:);
  
  [ok,M,F,tail] = samelength(M,F,tail);
  if (~ok)
    error('  RatioProb: input vectors must correspond in length.');
  end;
    k = length(M);

  ratio = M./F;
  X = max([M,F]')';
  N = M+F;
  
  for i = 1:k                           % Significance levels
    switch(tail(i))
      case -1
        pr = binocdf(X,N,0.5);
      case 0
        pr = binocdf(X,N,0.5);
        if (pr > 0.5)
          pr = 1-pr;
        end;
        pr = 2*pr;
      case 1
        pr = 1-binocdf(X,N,0.5);
    end;
  end;
  
  signif = seqbonf(pr,alpha,1);         % Sequential Bonferroni decisions
 
  if (get_minN)                         % Minimum sample sizes
    minN = 5*ones(k,1);
    for i = 1:k                             
      if ((M(i)>F(i)&tail(i)>0) | (M(i)<F(i)&tail(i)<0) | (M(i)~=F(i)&tail(i)==0))
        p = 1;
        while (p>alpha)
          minN(i) = minN(i)+1;
          f = prbcount([ratio(i),1]/(ratio(i)+1),minN(i),[],0,1);
          switch(tail(i))
            case -1
              p = binocdf(f(1),minN(i),0.5);
            case 0
              p = binocdf(f(1),minN(i),0.5);
              if (p > 0.5)
                p = 1-p;
              end;
              p = p*2;
            case 1
              p = 1-binocdf(f(1),minN(i),0.5);
          end;
        end;
      else
        minN(i) = Inf;
      end;
    end;
  end;
  
  return;
  