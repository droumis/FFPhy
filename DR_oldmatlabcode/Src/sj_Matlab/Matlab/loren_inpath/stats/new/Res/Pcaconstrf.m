% PCACONSTRF: Objective function for PCACONSTR. 
%
%     Usage: [soln,scores] = pcaconstrf(X,nu1,nu2,nu3,iv,npc,loadtype)
%
%         X =           [n x q+p] matrix of independent+dependent variables.
%         nu1,nu2,nu3 = unused.
%         iv =          vector of indices of independent variables.
%         npc =         number of leading principal components to be 
%                         returned.
%         loadtype =    boolean flag indicating the scaling for the loadings: 
%                           0: vector correlations;
%                           1: regression coefficients;
%                           2: squared loadings sum to unity.
%         --------------------------------------------------------------------
%         soln = row vector containing [loadings(:); percvar(:)]'
%           loadings =    [p x npc] matrix of principal components (columns).
%           percvar =     [p x 1] vector of percents variance-explained
%                         for principal components.
%         scores =      [n x npc] matrix of PCA scores (columns).
%

% Little, RJA & DB Rubin. 1987. Statistical Analysis with Missing Data. Wiley.
%   Section 6.5, pp. 112-115.

% RE Strauss, 6/1/00
%   11/28/00 - convert to multiple regression rather than matrix sweep.
%    2/ 9/00 - pass 'npc' to pcacov();
%              return regression parameters and statistics.

function [soln,scores,b,stats] = pcaconstrf(X,nu1,nu2,nu3,iv,npc,loadtype)
  Y = X;
  Y(:,iv) = [];                           % Matrix of dependent variables

  [b,stats,pred,resid] = linregr(X(:,iv),Y);  % Multiple regression
  [loadings,percvar,scores] = pcacov(resid,npc);

%  m = mean(X);                            % Means
%  C = cov(X);                             % Covariances

%  [D,E,residC] = sweepreg(m,C,iv);        % Sweep independent vars from matrix

%  [evects,evals] = eigen(residC);         % PCA of residual cov matrix
%  [loadings,scores] = loadscrs(Y,evects,npc,loadtype);

%  percvar = 100 * evals / sum(evals);     % Percents variance
%  percvar = percvar(1:npc);               % Retain subset

  soln = [loadings(:); percvar(:)]';

  return;

