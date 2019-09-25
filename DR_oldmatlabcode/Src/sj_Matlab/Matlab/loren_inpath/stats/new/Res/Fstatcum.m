% FSTATCUM: Finds accumulation curves for Wright's FST values, both "average" 
%           (across alleles) and maximum (of all alleles), for increasing 
%           numbers of populations.  A minimum of three populations is needed.
%
%     Usage: [avgfst,maxfst,seq] = fstatcum(alleles,popl,{permute},{doplot})
%
%           alleles = [n x 2c] matrix of allele identifiers for n obs and c loci.
%           popl =    [n x 1] vector of group-membership identifiers for k groups.
%           permute = optional boolean flag indicating, if true, that the 
%                       sequence of populations is to be randomized before 
%                       accumulation.
%           doplot =  optional boolean flag indicating, if true, that plot
%                       of accumulation curves is to be produced [default=0].
%           ------------------------------------------------------------------
%           avgfst =  vector (length k-1) of FST values, averaged across loci, 
%                       for the first 2,3,...,k populations.
%           maxfst =  vector (length k-1) of maximum FST (across loci) values 
%                       for the first 2,3,...,k populations.
%           seq =     sequence of group identifiers used in accumulation.
%

% RE Strauss, 8/22/97
%   5/25/00 - changed sequence of input arguments; changed handling of arguments 
%               for fstat(); added plot option; added random permutation of popls.

function [avgfst,maxfst,seq] = fstatcum(alleles,popl,permute,doplot)
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

  avgfst = zeros(k-1,1);                % Allocate output matrices
  maxfst = zeros(k-1,1);

  i = find(popl==popid(1));             % Indices of first population
  for npop = 2:k                        % Take first 'npop' population
    i = [i; find(popl==popid(npop))];     % Add indices of next population
    p = popl(i);                          % Isolate first 'npop' popls
    a = alleles(i,:);
    fval = fstat(a,p);                    % Find F-statistics
    nloci = size(fval,1)-1;
    fst = fval(nloci+1,3);
    avgfst(npop-1) = fst;                 % Stash results
    maxfst(npop-1) = max(fval(1:nloci,3));
  end;

  if (doplot)
    figure;
    npop = 2:k;
    plot(npop,avgfst,'k');
    hold on;
    plot(npop,maxfst,'k');
    hold off;
    putbnd(npop,maxfst);
    v = [1.8, 1.05*k, 0, 1.05*max(maxfst)];
    axis(v);
    putxlab('Number of populations');
    putylab('Overall and maximum Fst values');
  end;

  return;
