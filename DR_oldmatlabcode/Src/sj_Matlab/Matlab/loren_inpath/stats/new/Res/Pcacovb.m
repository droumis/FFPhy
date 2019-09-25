% PCACOVB: Objective function for bootstrapping pcacov().  Returns a single 
%           row vector containing loadings and percvar values.
%
%     Usage: retstr = pcacovb(X,not_used1,not_used2,not_used3,npc,loadtype,origload)
%
%         X =           [n x p] data matrix (obs x vars).
%         grps =        row or column vector of group identifiers.
%         npc =         number of leading discriminant functions for
%                         which scores are desired (default = groups-1).
%         loadtype =    optional boolean flag indicating the scaling for the 
%                         loadings: 
%                           0: vector correlations [default];
%                           1: regression coefficients;
%                           2: squared loadings sum to unity.
%         origload =    [p x ndf] matrix of loadings from original analysis.
%         --------------------------------------------------------------------------
%         retstr =      row vector containing loadings and percvar results.
%

% RE Strauss, 11/21/99, modified from discrimb.m
%   5/2/00 -    isolated scores & loadings in LOADSCRS.

function retstr = pcacovb(X,nu1,nu2,nu3,npc,loadtype,origload)
  covmat = cov(X);                      % Covariance matrix
  [evects,evals] = eigen(covmat);
  percvar = 100 * evals / sum(evals);   % Percents variance
  percvar = percvar(1:npc);             % Retain subset

  loadings = loadscrs(X,evects,npc,loadtype);

  for d = 1:npc                         % Check if direction consistent with original
    if (corr(loadings(:,d),origload(:,d),2)<0)
      loadings(:,d) = -loadings(:,d);
    end;
  end;

  retstr = [loadings(:)' percvar'];

  return;

