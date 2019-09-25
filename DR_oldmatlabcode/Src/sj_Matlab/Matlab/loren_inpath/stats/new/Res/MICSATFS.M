% MICSATFS: Bootstrapped measures of population subdivision based on microsatellite
%           loci, analogous to F-statistics but based on nested ANOVA of allele sizes.
%
%     Usage: [R,R_std,R_low,R_high] = micsatfs(alleles,grps,{iter},{CI_level})
%
%           alleles =   [n x 2c] matrix of allele values for n obs and c loci.
%           grps =      [n x 1] vector of group-membership identifiers.
%           iter =      number of bootstrap iterations [default=0].
%           CI_level =  level of confidence intervals [default=95].
%           --------------------------------------------------------------------------
%           R =         [c+1 x 7] matrix of modified "F-statistics".  Rows 1 to c
%                         represent loci, row c+1 represents average across loci.
%                         Col 1 = W (genetic variance among alleles within individuals)
%                             2 = I (genetic variance among individuals within groups)
%                             3 = S (genetic variance among groups)
%                             4 = T (total genetic variance)
%                             5 = R(IS) = W/T
%                             6 = R(IT) = I/T
%                             7 = R(ST) = S/T
%           R_std =     corresponding bootstrapped standard errors
%           R_low =     corresponding low confidence limits
%           R_high =    corresponding high confidence limits
%

% Slatkin, M.  1995.  A measure of population subdivision based on microsatellite
%   allele frequencies.  Genetics 139:457-462.

% RE Strauss, 3/9/97
%   11/29/99 - changed calling sequence.
%   4/27/00 -  convert NaN's to zeros.

function [R,R_std,R_low,R_high] = micsatfs(alleles,grps,iter,CI_level)
  if (nargin < 3) iter = []; end;
  if (nargin < 4) CI_level = []; end;

  if (isempty(iter))                  % Default input arguments
    iter = 0;                         %   No bootstrapping
  end;
  if (isempty(CI_level))
    CI_level = 0.95;                  %   95% CI
  elseif (CI_level>1)
    CI_level = CI_level / 100;
  end;

  CI_level = CI_level + (1-CI_level)/2;  % Center the CI

  [nobs,i] = size(alleles);           % Number of observations and loci
  nloci = floor(i/2);

  ncol = 7;
  nrow = nloci+1;
  R = zeros(nrow,ncol);

  [r,c] = find(~isfinite(alleles));   % Convert NaN to zeros
  if (~isempty(r))
    for i = 1:length(r)
      alleles(r(i),c(i)) = 0;
    end;
  end;

  for c = 1:nloci                     % Get genetic variance components by locus
    X = alleles(:,(2*c-1:2*c));         % Isolate current locus

    [W,I,S,T] = micsatvs(X,grps);

    R(c,1:4) = [W,I,S,T];               % Stash variance components in matrix
    R(c,5:7) = [W,I,S]/T;               % R(IS), R(IT), R(ST)
  end;

  R(nrow,1:4) = mean(R(1:nloci,1:4));
  R(nrow,5:7) = R(nrow,1:3)/R(nrow,4);

  if (iter == 0)                      % If no bootstrap, quit
    R_std = [];
    R_low = [];
    R_high = [];
    return;
  end;

  % Bootstrap genetic-variance statistics

  iter = iter+1;                       % Include initial solution
  R_distrib = zeros(iter,nrow*ncol);   % Allocate matrix for booted dists
  R_distrib(1,:) = R(:)';              % Stash observed distances

  for b = 2:iter                       % Bootstrap iterations
    A = bootsamp(alleles,grps);          % Get bootstrapped sample

    for c = 1:nloci                      % Get genetic variances by locus
      X = A(:,(2*c-1:2*c));              % Isolate current locus
      [W,I,S,T] = micsatvs(X,grps);

      R(c,1:4) = [W,I,S,T];               % Stash variance components in matrix
      R(c,5:7) = [W,I,S]/T;               % R(IS), R(IT), R(ST)
    end;

    R(nrow,1:4) = mean(R(1:nloci,1:4));
    R(nrow,5:7) = R(nrow,1:3)/R(nrow,4);

    R_distrib(b,:) = R(:)';             % Stash solution
  end;  % Bootstrap

  R = reshape(R_distrib(1,:),nrow,ncol);     % Restore original solution
  R_std = reshape(std(R_distrib),nrow,ncol); % Standard errors
  R_distrib = sort(R_distrib);               % Sort columns of distribution separately

  ci = bootci(R_distrib,R(:)',[],CI_level);  % Confidence limits

  R_low =  reshape(ci(1,:),nrow,ncol);
  R_high = reshape(ci(2,:),nrow,ncol);  

  return;   
