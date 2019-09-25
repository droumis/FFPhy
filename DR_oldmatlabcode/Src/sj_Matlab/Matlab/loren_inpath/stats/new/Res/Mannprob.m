% MANNPROB: Tail probabilities from the Mann-Whitney-Wilcoxon distribution.
%
%     Usage: pr = mannprob(U,N,{tail},{doplot})
%
%         U =       observed test-statistic values for the two groups; if only 
%                     a single value is provided, it is assumed to be the value 
%                     for the first group and that for the second group is 
%                     calculated from the sample sizes.
%         N =       2-element vector of sample sizes for the two groups.
%         tail =    optional flag indicating the "direction" of the test for 
%                     groups g1 & g2:
%                       -1: left-tailed  (H0: g1 <= g2)
%                        0: 2-tailed     (H0: g1 =  g2)  [default]
%                       +1: right-tailed (H0: g1 >= g2)
%         doplot =  optional boolean flag indicating, if true, that a plot of 
%                     the null distribution is to be produced.
%         ---------------------------------------------------------------------
%         pr =    significance level of the test.
%

% Based on Applied Statistics algorithm AS 62, Appl. Statist. 22(2), 1973.

% RE Strauss, 4/26/99
%   4/10/01 - miscellaneous improvements.

function pr = mannprob(U,N,tail,doplot)
  if (nargin < 3) tail = []; end;
  if (nargin < 4) doplot = []; end;

  if (isempty(tail))
    tail = 0;
  end;
  if (isempty(doplot))
    doplot = 0;
  end;

  if (length(U)==1)
    U(2) = prod(N) - U(1);
  end;

  U = U+1;                                % Adjust for U range = 0-prod(N)
  
  nmin = min(N);
  nmax = max(N);
  Umin = min(U);
  Umax = max(U);

  mn1 = prod(N)+1;
  n1 = nmax+1;
  freq = [ones(n1,1); zeros(mn1-n1,1)];

  lwrk = floor((mn1+1)/2 + nmin);
  work = zeros(lwrk,1);

  % Generate successively higher-order distributions

  in = nmax;
  for i = 2:nmin
    in = in+nmax;
    n1 = in+2;
    l = 1 + in/2;
    k = i;

    % Generate complete distribution from outside inwards

    for j = 1:l
      k = k+1;
      n1 = n1-1;
      summ = freq(j) + work(j);
      freq(j) = summ;
      work(k) = summ - freq(n1);
      freq(n1) = summ;
    end;
  end;

  freq = freq/sum(freq);                % Make distribution relative

  if (doplot)
    histgramb([0:prod(N)]',freq,[],[],1);
    putxlab('Sum of ranks');
    putylab('Frequency');
    puttitle(sprintf('Mann-Whitney-Wilcoxon distribution (N=%d)',sum(N)),14);
  end;

  % Cumulative frequency distribution

  cumfreq = cumsum(freq);

  if (tail>0)                         % Right-tailed probability
    if (isintegr(U(1)))
      pr = 1-cumfreq(U(1));
    else
      pr = 1-mean(cumfreq([floor(U(1)),ceil(U(1))]));
    end;
  elseif (tail<0)                     % Left-tailed probability
    if (isintegr(U(1)))
      pr = cumfreq(U(1));
    else
      pr = mean(cumfreq([floor(U(1)),ceil(U(1))]));
    end;
  else                                % Two-tailed probability
    if (isintegr(U(1)))
      pr = cumfreq(Umin) + 1-cumfreq(Umax);
    else
      pr = mean(cumfreq([floor(Umin),ceil(Umin)]));
      pr = pr + 1-mean(cumfreq([floor(Umax),ceil(Umax)]));
    end;
  end;

  return;
