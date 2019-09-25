% TREESPACE: Plots tree on a 2D character space, given an ancestor function.
%
%     Usage: treespace(crds,anc,labels,fontsize)
%
%           crds =     coordinates for terminal taxa and nodes, including root.
%           anc =      vector specifying ancestor function, with ancestor of 
%                        root specified as 0.
%           labels =   optional matrix of labels (rows) for terminal taxa
%                        [default = numeric labels].
%           fontsize = optional font size for labels [default = 10].
%                        

% RE Strauss, 4/24/00

function treespace(crds,anc,labels,fontsize)
  if (nargin < 3) labels = []; end;
  if (nargin < 4) fontsize = []; end;

  if (isempty(fontsize))
    fontsize = 10;
  end;

  n_anc = length(anc);
  n_taxa = ceil(n_anc/2);
  n_nodes = n_anc - n_taxa;
  n_edges = n_anc - 1;

  [n,p] = size(crds);
  if (p ~= 2)
    error('  TREESPACE: coordinates must be two-dimensional');
  end;
  if (n ~= n_anc)
    error('  TREESPACE: coordinate matrix and ancestor vector not compatible');
  end;

  if (~isempty(labels))
    [n,p] = size(labels);
    if (n ~= n_taxa)
      error('  TREESPACE: labels valid only for terminal taxa');
    end;
  end;

  % Get list of ancestors, descendants, and node levels by recursion
  root = find(anc==0);
  if (isempty(i))
    error('  TREESPACE: ancestor of root node must be specified as 0');
  end;

  [links,anc,curr_node,level] = treedivd([],anc,root,0);
  anc = links(:,1);          % Ancestral ends of links
  desc = links(:,2);         % Descendant ends of links
  n_anc = n_anc-1;

  scatter(crds(1:n_taxa,:));

  deltax = 0.02*range(crds(:,1));
  for t = 1:n_taxa
    if (isempty(labels))
      h = text(crds(t,1)+deltax,crds(t,2),num2str(t));
    else
      h = text(crds(t,1)+deltax,crds(t,2),num2str(labels(t,:)));
    end;
    set(h,'fontsize',fontsize);
  end;

  hold on;
  plot(crds(root,1),crds(root,2),'ok');
  for i = 1:n_anc
    b = [anc(i) desc(i)];
    plot(crds(b,1),crds(b,2),'k');
  end;
  hold off;

  return;
