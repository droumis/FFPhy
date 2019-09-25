% LNKTOANC: Recursive algorithm to convert set of links (internodes) to an ancestor 
%           function and corresponding set of branch lengths.
%
%     Usage: [anc,brlen,links,curr_node] = lnktoanc(anc,brlen,links,curr_node)
%
%           anc =       ancestor function, with ancestor of root node 
%                         specified as 0; initially pass as null vector.
%           brlen =     corresponding vector of branch lengths; initially pass as null 
%                         vector.
%           links =     3-column matrix of ancestors & descendants (not necessarily in 
%                         that order) and branch lengths. 
%           curr_node = current node; initially pass as root node.
%

function [anc,brlen,links,curr_node] = lnktoanc(anc,brlen,links,curr_node)
  if (isempty(links))
    return;
  end;
  nlinks = size(links,1);           % Current number of links

  if (isempty(anc))                 % If first call, allocate 'anc' & 'brlen' vectors
    anc = zeros(1,nlinks+1);
    brlen = anc;
  end;

  lnk = [links; links(:,[2 1 3])];  % Double links matrix, switching cols 1,2 in second copy

  i = find(lnk(:,1)==curr_node);    % Find descendants of current node
  if (length(i)>2)
    error(sprintf('  LNKTOANC: tree must be dichotomous at node %d',curr_node));
  end;
  if (isempty(i))
    return;
  end;
  desc = lnk(i,2);

  anc(desc) = [curr_node; curr_node];       % Stash values of ancestor and branch length
  brlen(desc) = lnk(i,3);                   %   for each of descendants

  j = (i > nlinks);                         % Remove links from input matrix
  i(j) = i(j) - nlinks;
  links(i,:) = [];

  [anc,brlen,links,curr_node] = lnktoanc(anc,brlen,links,desc(1));
  [anc,brlen,links,curr_node] = lnktoanc(anc,brlen,links,desc(2));

  return;
