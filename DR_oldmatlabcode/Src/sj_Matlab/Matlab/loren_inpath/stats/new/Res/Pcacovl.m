% PCACOVL: Calculate loadings and scores for a principal component analysis.
%
%     Usage: [loadings,scores] = pcacovl(X,evects,npc,loadtype)
%
%         X =         [n x p] data matrix.
%         evects =    [p x npc] matrix of eigenvectors.
%         npc =       optional number of leading principal components to be 
%                       returned [default = number of significant eigenvalues 
%                       based on broken-stick null model].
%         loadtype =  optional boolean flag indicating the scaling for the 
%                       loadings: 
%                         0: vector correlations [default];
%                         1: regression)coefficients;
%                         2: squared loadings sum to unity.
%         -------------------------------------------------------------------
%         loadings =  [p x npc] matrix of loadings.
%         scores =    [n x npc] matrix of scores.
%

% RE Strauss, 5/2/00

function [loadings,scores] = pcacovl(X,evects,npc,loadtype)
  [n,p] = size(X);
  scores = score(X,evects,npc);             % Scores for subset of PCs

  switch(loadtype)                          % Loadings
    case 0,                                   % Vector correlations
      loadings = corr(X,scores);            

    case 1,                                   % Regression coefficients
      loadings = zeros(p,npc);
      for ipc = 1:npc
        S = [ones(n,1) scores(:,ipc)];
        for ip = 1:p
          y = X(:,ip);
          b = inv(S'*S)*S'*y;
          loadings(ip,ipc) = b(2);
        end;
      end;

    case 2,
      loadings = sumsqscale(evects(:,1:npc)); % Squared coeffs sum to unity
  end;

  return;
