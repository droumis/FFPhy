% RANDTREE: Provides the ancestor function for a random tree based on Rohlf's M.
%
%     Usage: anc = randtree(T)
%
%           T = number of terminal taxa.
%           ----------------------------
%           anc = ancestor function.
%

% RE Strauss, 12/26/99

function anc = randtree(T)
  M = floor(rand*treenum(T));
  anc = treemanc(T,M);

  return;
