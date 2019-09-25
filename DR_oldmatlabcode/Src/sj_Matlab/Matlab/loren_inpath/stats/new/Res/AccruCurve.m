% AccruCurve: Given a vector of object identifiers for a set of samples, produces 
%             an average accrument curve by randomizing the sequence of samples.
%             Typical uses: samples are localities and objects are species; sample
%             are genotypes and objects are alleles.
%
%     Usage: asym = accrucurve(x,{iter},{noplot})
%
%         x =       vector of object identifiers.
%         iter =    optional number of permutation iterations [default = 1000].
%         noplot =  optional boolean flag indicating, if true, that plot is to be
%                     suppressed.
%         -----------------------------------------------------------------------
%         asym =    asymptote predicted by negative-exponential fit to the
%                     observed accrument curve.
%

% RE Strauss, 2/28/02

function asym = accrucurve(x,iter)
  if (nargin<2) iter = []; end;
  
  if (isempty(iter))
    iter = 1000;
  end;

  n = length(x);
  c = ones(iter,n);

  for it = 1:iter
    xp = x(randperm(n));
    for i = 2:n
      if (isin(xp(i),xp(1:(i-1))))
        c(it,i) = c(it,i-1);
      else
        c(it,i) = c(it,i-1)+1;
      end;
    end;
  end;
  
  cmean = mean(c);
  
  [P,pS,tp,r,curve] = vonbert([0:n],[0,cmean]);
  asym = P(7);
  
  plot([0:n],[0,cmean],'k',linspace(0,n)',curve,'k:');
  legend('Observed','Predicted',2);
  putxlab('Number of samples');
  putylab('Cumulative number of unique objects');

  return;
  