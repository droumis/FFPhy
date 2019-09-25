% MICSATGD: Genetic distances based on microsatellite loci, alleles of which
%           vary in number of repeats of simple DNA sequence.  Bootstraps
%           the distance matrix and UPGMA cluster analysis.
%
%     Usage: [dist,CI,clst_supt,tot_supt,topo] = 
%                 micsatgd(grps,alleles,{iter},{supt_cut},{CI_level})
%
%           grps =      [n x 1] vector of group-membership identifiers.
%           alleles =   [n x 2c] matrix of allele frequencies for n obs and c loci.
%           iter =      number of bootstrap iterations (default=0).
%           supt_cut =  threshhold for cluster frequency support [0-1] (default=0.5).
%           CI_level =  level of confidence intervals (default=95).
%           --------------------------------------------------------------------------
%           dist =      [k x k] symmetric matrix of intergroup genetic distances.
%                         If bootstrapping is done, the lower triangular matrix 
%                         contains bootstrap standard errors.
%           CI =        [k x k] matrix of lower (lower triangular matrix) and
%                         upper (upper triangular matrix) confidence limits of
%                         genetic distances for k groups.
%           clst_supt = [(k-2) x k] matrix, with one row for all but the base
%                         node, specifying group membership and frequency support at 
%                         each node; col 1 is the freq support, cols 2:k are clusters.
%           tot_supt =  [(k-2) x k] matrix specifying group membership and frequency support 
%                         for any bootstrapped clusters having frequency support 
%                         greater than or equal to the specified threshhold; 
%                         col 1 is the freq support, cols 2:k are clusters.
%           topo =      [(k-1) x 4] matrix summarizing dendrogram topology:
%                         col 1 = 1st OTU/cluster being grouped at current step
%                         col 2 = 2nd OTU/cluster
%                         col 3 = ID of cluster being produced
%                         col 4 = distance at node
%

% Goldstein, D.B., A.R. Linares, L.L. Cavalli-Sforza, and M.W. Feldman.
%   1995.  An evaluation of genetic distance for use with microsatellite
%   loci.  Genetics 139:463-471.

% RE Strauss, 3/11/97
%   9/20/99 - update handling of null input arguments.
%   4/27/00 - convert NaN's to zeros.
%   5/14/00 - changed from micsatfr() to allelefreq().

function [dist,CI,clst_supt,tot_supt,topo] = ...
            micsatgd(grps,alleles,iter,supt_cut,CI_level)
  if (nargin < 3) iter = []; end;
  if (nargin < 4) supt_cut = []; end;
  if (nargin < 5) CI_level = []; end;

  if (isempty(iter))                  % Default input arguments
    iter = 0;                           % No bootstrapping
  end;
  if (isempty(supt_cut))
    supt_cut = 0.5;                     % Frequency-support threshhold = 0.5
  end;
  if (isempty(CI_level))
    CI_level = 0.95;                    % 95% CI
  elseif (CI_level>1)
    CI_level = CI_level / 100;
  end;

  CI_level = CI_level + (1-CI_level)/2;  % Center the CI

  grpid = uniquef(grps);              % Group identifiers
  ngrps = length(grpid);              % Number of groups
  [nobs,i] = size(alleles);           % Number of observations and loci
  nloci = floor(i/2);

  [r,c] = find(~isfinite(alleles));   % Convert NaN to zeros
  if (~isempty(r))
    for i = 1:length(r)
      alleles(r(i),c(i)) = 0;
    end;
  end;

  freq = allelefreq(grps,alleles);    % Allele frequencies
  dist = micsatds(grps,nloci,freq);   % Genetic distance matrix
  [topo,clst_supt] = upgma(dist);     % Plot dendrogram

  if (iter == 0)                      % If no bootstrap, quit
    clst_supt = [];
    tot_supt = [];
    CI = [];
    return;
  end;

  % Bootstrap genetic distances and cluster analyses

  clst_supt = [ones(size(clst_supt,1),1) clst_supt];  % Attach col 1,
  tot_supt = clst_supt;                               %   = freqs of cluster support
  [sr,sc] = size(tot_supt);           % Sizes of support matrices 
  cr = sr;

  iter = iter+1;                      % Include initial solution
  dist_distrib = zeros(iter,ngrps*ngrps); % Allocate matrix for booted dists
  dist_distrib(1,:) = dist(:)';       % Stash observed distances

  for b = 2:iter                 % Bootstrap iterations
    A = bootsamp(alleles,grps);
    freq = micsatfr(grps,A);            % Allele frequencies
    d = micsatds(grps,nloci,freq);      % Genetic distance matrix
    dist_distrib(b,:) = d(:)';          % Stash entire dist matrix as row
    [t,s] = upgma(d,[],1);              % Get UPGMA support matrix

    for i = 1:cr                      % Update cluster-support information
      for j = 1:cr                        % Update cluster if in clst_supt
        if (all(s(i,:)==clst_supt(j,2:sc)))
          clst_supt(j,1) = clst_supt(j,1) + 1;
        end;
      end;

      found = 0;
      for j = 1:sr                        % Check if cluster is in support matrix
        if (all(s(i,:)==tot_supt(j,2:sc)))
          tot_supt(j,1) = tot_supt(j,1) + 1;
          found =  1;
        end;
      end;
      if (~found)                         % If not, append
        tot_supt = [tot_supt; 1 s(i,:)];
        sr = sr + 1;
      end;
    end;
  end;  % Bootstrap

  clst_supt(:,1) = clst_supt(:,1) / iter; % Convert counts to proportions
  tot_supt(:,1) = tot_supt(:,1) / iter;

  [Y,i] = sort(-clst_supt(:,1));          % Sort cluster-support freqs in
  clst_supt = clst_supt(i,:);             %   descending sequence

  indx = find(tot_supt(:,1) >= supt_cut); % Keep support-freqs >= threshhold
  if (isempty(indx))
    tot_supt = [];
  else
    tot_supt = tot_supt(indx,:);
    [Y,i] = sort(-tot_supt(:,1));         % Sort support freqs in descending seq
    tot_supt = tot_supt(i,:);
  end;

  dist_distrib = sort(dist_distrib);      % Sort columns of dists separately
  dist_stderr = reshape(std(dist_distrib),ngrps,ngrps); % Standard errors
  for i=1:(ngrps-1)                       % Plug stderrs into lower triangular
    for j=(i+1):ngrps                     %   distance matrix
      dist(j,i) = dist_stderr(j,i);
    end;
  end;

  high_limit = ceil((iter*CI_level));     % CI indices
  low_limit = iter - high_limit +1;

  high_dist = reshape(dist_distrib(high_limit,:),ngrps,ngrps);  % Pull out
  low_dist = reshape(dist_distrib(low_limit,:),ngrps,ngrps);    %   confid limits

  CI = zeros(ngrps,ngrps);                % Allocate CI matrix
  for i=1:(ngrps-1)
    for j=(i+1):ngrps
      CI(i,j) = high_dist(i,j);
      CI(j,i) = low_dist(j,i);
    end;
  end;

  return;
