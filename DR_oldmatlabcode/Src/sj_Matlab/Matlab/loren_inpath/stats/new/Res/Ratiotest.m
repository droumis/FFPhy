% RATIOTEST: Tests whether a given ratio value is singificantly different from 
%            unity, as a function of sample size, based on either a one-tailed
%            binomial test or Fisher's exact test.  Varies sample size from 10 
%            to Nmax.  Note that the binomial probability treats the null ratio 
%            (1.0) as fixed, while the others treat it as having sampling 
%            variation.
%
%     Usage: [pval,N,realized] = ratiotest(ratio,{Nmax},{noplots})
%                 or
%            [pval,N,realized] = ratiotest(ratio,{[Nmin,Nmax]},{noplots})
%
%         ratio =       scalar value of ratio, interpreted as ratio of T to F 
%                         for a given condition.
%         Nmax =        optional maximum sample size [default = 100].
%         [Nmin,Nmax] = optional vector of minimum and maximum sample sizes 
%                         [default = 10,100].
%         noplots =     optional boolean value indicating, if true, that plots  
%                         are to be suppressed [default = 0].
%         ---------------------------------------------------------------------
%         pval =        3-column matrix of significance levels:
%                         col 1 - binomial tail test;
%                             2 - normal approximation;
%                             3 - Fisher's exact test.
%         N =           corresponding vector of sample sizes.
%         realized =    corresponding vector of ratio values, constrained by
%                         integer sample sizes.
%

% RE Strauss, 10/16/01

function [pval,N,realized] = ratiotest(ratio,Nmax,test,noplots)
  if (nargin < 2) Nmax = []; end;
  if (nargin < 3) test = []; end;
  if (nargin < 4) noplots = []; end;

  if (isempty(Nmax))
    Nmin = 10;
    Nmax = 100;
  end;
  if (isempty(test))
    test = 1;
  end;
  if (isempty(noplots))
    noplots = 0;
  end;

  if (length(Nmax)==1)
    Nmin = 10;
  else
    Nmin = min(Nmax);
    Nmax = max(Nmax);
  end;

  N = [Nmin:Nmax]';
  nvals = length(N);

  pFisher = zeros(nvals,1);
  pBinom =  zeros(nvals,1);
  pNormal = zeros(nvals,1);
  realized = zeros(nvals,1);

  for i = 1:nvals
    n = N(i);

    f = prbcount([ratio,1]/(ratio+1),n,[],0,1);
    realized(i) = f(1)/f(2);
    f = -sort(-f);

    e = prbcount([1,1]/2,n,[],0,1);
    e = -sort(-e);

    pFisher(i) = continex([f' e']);
    pBinom(i) =  1-binocdf(f(1),n,0.5);
    [d,s2d,ci,pNormal(i)] = propdiff(f(2),f(1),e(2),e(1));
  end;

  if (~noplots)
    figure;
    plot(N,realized,'k',[min(N),max(N)],[ratio,ratio],'k');
    putbnd(N,realized);
    putxlab('Sample size');
    putylab('Realized ratio');

    figure;
    plot(N,pBinom,'-',N,pNormal,'-',N,pFisher,'-',...
         [min(N),max(N)],[0.05 0.05],'--k');
    putbnd([N;N;N],[pBinom;pNormal;pFisher]);
    legend('Binomial','Normal','Fisher exact');
    putxlab('Sample size');
    putylab('Probability under H0');
  end;

  pval = [pBinom pNormal pFisher];

  return;
 