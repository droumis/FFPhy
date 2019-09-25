% HWMULTLL: Log-likelihood statistic for Hardy-Weinberg proportions for 
%           m alleles at a single locus.
%
%     Usage: G2 = hwmultll(alleles)
%
%           alleles = [n x 2] matrix of allele identifiers for n observations.
%           ------------------------------------------------------------------
%           G2 =      likelihood-ratio statistic.
%

% RE Strauss, 5/20/00

function G2 = hwmultll(alleles)
  [n,p] = size(alleles);
  if (p~=2)
    error('  HWMULTLL: input matrix must have two columns.');
  end;

  [ua,fa] = uniquef(alleles);         % Allele identities and frequencies
  m = length(ua);                     % Number of alleles

  if (m==1)                           % Return if monomorphic
    G2 = 0;
    return;
  end;

  g = zeros(m,m);                     % Table of genotype frequencies
  for i = 1:n
    a1 = find(ua==alleles(i,1));
    a2 = find(ua==alleles(i,2));
    g(a1,a2) = g(a1,a2)+1;
    if (a1~=a2)
      g(a2,a1) = g(a2,a1)+1;
    end;
  end;
  homo = diag(g);
  het = trilow(g);

  i = find(homo==0);
  if (~isempty(i))
    homo(i) = [];
  end;
  i = find(het==0);
  if (~isempty(i))
    het(i) = [];
  end;

  lnL0 = sum(fa.*log(fa./(2*n))) + sum(het*log(2));
  lnL1 = sum(homo.*log(homo/n)) + sum(het.*log(het/n));

  G2 = -2*(lnL0 - lnL1);

  return;

