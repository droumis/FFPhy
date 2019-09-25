% CONTMAPF: Objective function for CONTMAPI(), which iteratively maps a 
%           character onto a tree.
%
%     Usage: treelen = contmapf(xnode,x,brlen,k,links)
%
%           xnode =  [m x 1] vector of node values for m internal nodes.
%           x =      [n x 1] vector of terminal values for n taxa.
%           w =      [2n-2 x 1] vector of branch weights, in same order as links.
%           k =      Minkowski k
%           links =  [2n-2 x 1] list of branches of tree:
%                      col 1: ancestor nodes.
%                          2: descendant nodes.
%           ---------------------------------------------------------------------
%           treelen = sum of squared branch lengths of tree.
%

% RE Strauss, 6/23/98
%   12/17/99 - added Minkowski k.
%    5/22/01 - use reciprocal weights^k.

function treelen = contmapf(xnode,x,w,k,links)
  x = [x; xnode];                     % Concatenate node & terminal taxon values
  a = links(:,1);                     % Isolate ancestors & descendants of branches
  d = links(:,2);

  brlen = abs(x(a)-x(d));
  treelen = ((1./w(:)').^k * (brlen.^k))';
%  treelen = ((1./w(:)') * (brlen.^k))';

  return;

