% HWMULTINOM: Multinomial probability for the Hardy-Weinberg test for m alleles 
%             at 1 locus (Levene 1949).
%
%     Usage: pr = hwmultinom(alleles)
%
%             alleles = [n x 2] matrix of allele identifiers for n observations.
%             ------------------------------------------------------------------
%             pr =      probability of obtaining the sample, conditional on 
%                         observed allele counts and Hardy-Weinberg genotype 
%                         proportions.
%

% Levene, H. 1949. On a matching problem arising in genetics.  Annals of 
%   Mathematical Statistics 20:91-94.

% RE Strauss, 5/19/00

function pr = hwmultinom(alleles)
  [n,p] = size(alleles);
  if (p~=2)
    error('  HWMULTINOM: input matrix must have two columns.');
  end;

  [ua,fa] = uniquef(alleles);         % Allele identities and frequencies
  m = length(ua);                     % Number of alleles

  g = zeros(m,m);                     % Table of genotype frequencies
  for i = 1:n
    a1 = find(ua==alleles(i,1));
    a2 = find(ua==alleles(i,2));
    g(a1,a2) = g(a1,a2)+1;
    if (a1~=a2)
      g(a2,a1) = g(a2,a1)+1;
    end;
  end;
  ga = trilow(g);
  ng = length(ga);

  num = sum(log(2:n));
  den = num + sum(log([(n+1):(2*n)]));

  for i = 1:m
    num = num + sum(log(2:fa(i)));
  end;

  for i = 1:ng
    den = den + sum(log([2:ga(i)]));
  end;

  c = sum(ng) * log(2);

  pr = exp(c+num-den);

  return;
