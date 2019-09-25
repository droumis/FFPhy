% HETCUM: Finds accumulation curves for mean heterozygosity for increasing 
%         numbers of populations.
%
%     Usage: [meanhet,sehet,seq] = hetcum(alleles,popl,{permute},{doplot})
%
%         alleles = [n x 2c] matrix of allele identifiers for n obs and c loci.
%         popl =    [n x 1] vector of group-membership identifiers for k groups.
%         permute = optional boolean flag indicating, if true, that the 
%                     sequence of populations is to be randomized before 
%                     accumulation.
%         doplot =  optional boolean value indicating, if true, that plot
%                     of accumulation curves is to be produced [default=0].
%         ----------------------------------------------------------------------
%         meanhet = vector of mean heterozygosity estimates for increasing 
%                     numbers of groups (1,...,k).
%         sehet =   corresponding vector of standard errors of mean 
%                     heterozygosity estimates.
%         seq =     sequence of group identifiers used in accumulation.
%

% RE Strauss, 8/25/97
%   5/25/00 - changed sequence of input arguments; added plot option; added 
%               random permutation of popls..

function [meanhet,sehet,seq] = hetcum(alleles,popl,permute,doplot)
  if (nargin < 3) permute = []; end;
  if (nargin < 4) doplot = []; end;

  if (isempty(permute))
    permute = 0;
  end;
  if (isempty(doplot))
    doplot = 0;
  end;
  
  popid = uniquef(popl);                % Unique population identifiers
  k = length(popid);                    % Number of populations

  if (permute)                          % Randomly permute identifiers
    popid = popid(randperm(k));
  end;
  seq = popid;

  meanhet = zeros(k,1);                 % Allocate output matrices
  sehet = zeros(k,1);

  i = find(popl==popid(1));             % Indices of first population
  for npop = 1:k                        % Take first 'npop' population
    i = [i; find(popl==popid(npop))];     % Add indices of next population
    p = popl(i);                          % Isolate first 'npop' popls
    a = alleles(i,:);
    [hetzyg,avghetzyg] = heterozyg(a,p);  % Find F-statistics
    r = size(avghetzyg,1);
    meanhet(npop) = avghetzyg(r,2);       % Stash results
    sehet(npop) = avghetzyg(r,3);
  end;

  if (doplot)
    figure;
    npop = 1:k;
    plot(npop,meanhet,'k');
    hold on;
    plot(npop,meanhet+(2*sehet),'k:');
    plot(npop,meanhet-(2*sehet),'k:');
    hold off;
    v = [0, 1.05*k, 0, 1.05*max(meanhet+(2*sehet))];
    axis(v);
    putxlab('Number of populations');
    putylab('Mean heterozygosity & 95% CI');
  end;

  return;
