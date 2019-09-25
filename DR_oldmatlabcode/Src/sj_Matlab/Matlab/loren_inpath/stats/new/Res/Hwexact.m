% HWEXACT: Given counts for 3 genotypes for 1 locus, 2 alleles, calculates
%          the exact significance probability using Fisher's exact test.
%          Optionally provides counts expected under Hardy-Weinberg
%          equilibrium, the observed HW-disequilibrium coefficient, and the
%          probabilities corresponding to all possible combinations of genotype
%          numbers.
%
%          Approximate upper limit on sample size is 125.
%
%          See: Weir, BS. 1996. Genetic Data Analysis II. Sinauer.  Pp 98-101.
%          Based on: Haldane, JBS. 1954. An exact test for randomness of mating.
%                      Journal of Genetics 52:631-635.
%
%     Usage: [signif,expect,diseq,distrib] = hwexact(obs)
%
%          obs -     vector of counts for phenotypes A1-A1, A1-A2, A2-A2.
%          ----------------------------------------------------------------------
%          signif -  significance level.
%          expect -  expected Hardy-Weinberg counts.
%          diseq -   disequilibrium coefficient.
%          distrib - [d x 6] distribution of genotype numbers and associated
%                      probabilities, cumulative probabilities, and
%                      HW-disequilibrium values; sorted into sequence by
%                      increasing probability.
%

% RE Strauss, 5/23/96
%   5/20/00 - replaced factorial products with sums of logs to avoid overflow.

function [signif,expect,diseq,distrib] = hwexact(obs)
  signif = [];
  diseq = [];
  distrib = [];

  get_distrib = 0;
  if (nargout>=4)
    get_distrib = 1;
  end;

  n11 = obs(1);                       % Obs number of phenotypes
  n12 = obs(2);
  n22 = obs(3);

  n = sum(obs);                       % Number of genotypes
  nale = 2*n;                         % Number of alleles
  n1 = 2*n11 + n12;                   % Obs number of #1 alleles
  n2 = nale - n1;                     % Obs number of #2 alleles

  p = n1/nale;                        % Allele frequencies
  q = 1-p;
  expect = [n*p*p, n*2*p*q, n*q*q];   % expectected HW counts

  nalef = sum(log((n+1):nale));           % Factorials
  n12f = sum(log(1:n12));
  n1f = sum(log((n11+1):n1));
  n2f = sum(log((n22+1):n2));
  tohet = n12*log(2);

  prhet = exp(n1f+n2f+tohet-n12f-nalef);                  

  if (rem(n1,2))                      % Number of #1 alleles is
    minh = 1;                           % Odd
  else
    minh = 0;                           % Even
  end;

  signif = 0;
  d = 0;

  for n12 = minh:2:n1                 % Vary number of heterozygotes
    n11 = (n1-n12)./2;
    n22 = n-n11-n12;

    if (n22 >= 0)
      n12f = sum(log(1:n12));                % Factorials
      n1f = sum(log((n11+1):n1));
      n2f = sum(log((n22+1):n2));
      tohet = n12*log(2);

      prh = exp(n1f+n2f+tohet-n12f-nalef);                    

      if (prh<=prhet)
        signif = signif + prh;
      end;

      if (get_distrib)
        d = d+1;
        distrib = [distrib; n11 n12 n22 prh];
        diseq =   [diseq; (n11-expect(1))./n];
      end;
    end;
  end;

  if (get_distrib)
    [y,i] = sort(distrib(:,4));         % Sort by probability
    distrib = [distrib(i,:) cumsum(y) diseq(i)];
  end;

  diseq = (obs(1)-expect(1))./n;

  return;
