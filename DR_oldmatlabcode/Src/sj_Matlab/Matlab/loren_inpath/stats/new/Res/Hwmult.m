% HWMULT: Tests against Hardy-Weinberg proportions for m alleles at a single 
%         locus.  Defaults to the likelihood-ratio G-test; if iter>0, performs 
%         the permutation test of Pearson (1900; see Guo & Thompson 1992), but 
%         using the G2 statistic (Weir 1996:105-106) rather than the multinomial 
%         probability.  The chi-square approximation of the G2 statistic is not 
%         moderate to large numbers of alleles due to the underrepresentation of 
%         reliable for many heterozygote classes.
%
%     Usage: [pr,G2,df] = hwmult(alleles,{grps},{iter})
%
%         alleles = [n x 2c] matrix of allele identifiers, for n observations 
%                     and c loci.
%         grps =    optional grouping variable for k groups [default = null].
%         iter =    optional number of random permutations [default = 0].
%         -------------------------------------------------------------------
%         pr =      [k x c] matrix of probabilities for k groups and c loci.
%         G2 =      matching matrix of observed G-statistics.
%         df =      matching matrix of degrees-of-freedom for chi-square tests 
%                     (returned as null for permutation tests).
%         

% Pearson, K.  1900.  On the criterion that a given system of deviations from the
%   probable in the case of a correlated system of variables is such that it can be 
%   reasonably supposed to have arisen from random sampling.  Phil. Mag. 5th series,
%   50:157-175.
% Guo, SW & EA Thompson. 1992. Performing the exact test of Hardy-Weinberg 
%   proportion for multiple alleles.  Biometrics 48:361-372.

% RE Strauss, 5/21/00
%   6/18/02 - added Pearson 1900 reference to documentation.

function [pr,G2,df] = hwmult(alleles,grps,iter)
  if (nargin < 2) grps = []; end;
  if (nargin < 3) iter = []; end;

  [n,c2] = size(alleles);
  nloci = c2/2;
  if (~isintegr(nloci))
    error('  HWMULT: alleles matrix must have even number of columns.');
  end;

  if (isempty(grps))
    grps = ones(n,1);
  end;
  if (isempty(iter))
    iter = 0;
  end;

  if (length(grps) ~= n)
    error('  HWMULT: alleles matrix and grouping vector not compatible.');
  end;

  grpid = uniquef(grps);                % Groups
  ngrps = length(grpid);

  pr = NaN*ones(ngrps,nloci);
  G2 = zeros(ngrps,nloci);
  if (~iter)
    df = zeros(ngrps,nloci);
  else
    df = [];
  end;

  for ig = 1:ngrps
    g = find(grps==grpid(ig));            % Obs for current group

    for ic = 1:nloci
      a = alleles(g,[(2*ic-1):(2*ic)]);     % Isolate data for current locus
      
      i = find(~isfinite(sum(a'))); % Remove missing data
      if (~isempty(i))
        a(i,:) = [];
      end;

      G2_obs = hwmultll(a);                     % Observed G-statistic
      G2(ig,ic) = G2_obs;
G2

      if (~iter)
        ua = uniquef(a);
        m = length(ua);
        if (m>1)
          df_obs = m*(m-1)/2;
          pr(ig,ic) = 1-chi2cdf(G2_obs,df_obs);
          df(ig,ic) = df_obs;
        else
          pr(ig,ic) = 1;
          df(ig,ic) = 0;
        end;
      else
        if (G2_obs>0)
          p = 0;
          a = a(:);
          na = length(a);
          n = na/2;
          incr = 1./iter;

          for it = 1:iter                           % Iterate
            a = a(randperm(na));                      % Randomly permute
            prg = hwmultll([a(1:n) a(n+1:na)]);       % Multinomial prob of random genotypes
            if (prg >= G2_obs)
              p = p + incr;
            end;
          end;

          pr(ig,ic) = p;
        else
          pr(ig,ic) = 1;
        end;
      end;  % if iter
      save hwmult pr G2 df;
    end;  % for loci
  end;  % for groups

  return;
