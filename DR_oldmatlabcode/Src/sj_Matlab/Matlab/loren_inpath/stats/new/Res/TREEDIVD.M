% TREEDIVD: Recursive algorithm to find links (internodes) and levels on a
%           tree, given the ancestor function.
%
%     Usage: [links,anc,curr_node,level]
%               = treedivd(links,anc,curr_node,level)
%
%           links =     3-column matrix of ancestors, descendants, and levels;
%                         initially pass to function as a null vector, [].
%           anc =       ancestor function, with ancestor of root node 
%                         specified as 0.
%           curr_node = current node; initially pass as root node.
%           level =     current level; initially pass as zero.
%

% RE Strauss, 8/17/98

function [links,anc,curr_node,level] = treedivd(links,anc,curr_node,level)
  desc = find(anc==curr_node);    % Get descendants of node
  if (isempty(desc))              % If none, return
    return;
  end;
  lev = level + 1;                % Next level in tree

  links = [links; curr_node desc(1) lev];   % Add links to list
  links = [links; curr_node desc(2) lev];

  [links,anc,curr_node,level] = treedivd(links,anc,desc(1),lev);
  [links,anc,curr_node,level] = treedivd(links,anc,desc(2),lev);

  return;
