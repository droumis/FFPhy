% FSTAT: Estimation of Wright's F-statistics (FST [theta], FIS [f], FIT [F]), 
%        as described by Weir & Cockerham (1984).  Estimates F-statistics 
%        separately for each locus and averaged across loci.  Optionally 
%        bootstraps across individuals within populations.
%
%     Usage: [fval,CI_fval] = fstat(alleles,popl,{iter},{CI_level})
%
%         alleles = [n x 2c] matrix of allele identifiers for n obs and c loci.
%         popl =    [n x 1] vector of group-membership identifiers for k groups.
%         iter =      optional number of bootstrap iterations [default=0].
%         CI_level =  confidence level for bootstrapped confidence intervals 
%                       [default = 95].
%         ----------------------------------------------------------------------
%         fval =      [c+1 x 3] matrix of F estimates by locus:
%                       col 1 = Fis estimates
%                           2 = Fit estimates
%                           3 = Fst estimates
%                       The last row are the averages across loci.
%         CI_fval =   [c+1 x 6] matrix of lower and upper confidence limits:
%                       col 1-2 = bounds for Fis
%                           3-4 = bounds for Fit
%                           5-6 = bounds for Fst
%                       The last row are the averages across loci.
%

% Weir, BS and CC Cockerham. 1984. Estimating F-statistics for the analysis of 
%   population structure.  Evolution 38:1358-1370.

% RE Strauss, 8/22/00

function [fval,CI_fval] = fstat(alleles,popl,iter,CI_level)
  if (nargin < 3) iter = []; end;
  if (nargin < 4) CI_level = []; end;

  if (isempty(iter))
    iter = 0;
  end;
  if (isempty(CI_level))
    CI_level = 0.95;
  end;
  if (CI_level > 1)
    CI_level = CI_level/100;
  end;
  alpha = 1-CI_level;

  nloci = size(alleles,2)/2;
  fvect = fstatf(alleles,popl)';
  fval = reshape(fvect,nloci+1,3);
  
  CI_fval = [];
  if (iter)
    ci = bootstrp('fstatf',[1 0 0 0],iter,alpha,alleles,popl)';
    ci = ci(:);
    CI_fval = reshape(ci,nloci+1,6);
    CI_fval = CI_fval(:,[1 4 2 5 3 6]);
  end;

  return;
