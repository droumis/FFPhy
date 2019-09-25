% CLADPLOT: Plots a horizontal squared tree given an ancestor function 
%           and corresponding branch lengths.  The ancestor and branch lengths 
%           of the root should be specified as 0.
%
%     Usage: cladplot(anc,brlen)
%

% RE Strauss, 5/28/97

function cladplot(anc,brlen)
  n_anc = length(anc);
  n_taxa = ceil(n_anc/2);
  n_nodes = n_anc - n_taxa;
  n_edges = n_anc - 1;

  % Get list of ancestors, descendants, and node levels by recursion
  i = find(anc==0);
  if (isempty(i))
    error('  Ancestor of root node must be specified as 0');
  end;

  [links,anc,curr_node,level] = treedivd([],anc,i,0);
  anc = links(:,1);          % Ancestral ends of links
  desc = links(:,2);         % Descendant ends of links
  level = links(:,3);        % Levels in tree
  level = max(level)-level;  % Put level zero at tips
  max_level = max(level);

  t = find(desc <= n_taxa);  % Get indices to taxa
  taxon = desc(t);           % Get taxa

  x = zeros(n_anc,1);        % Allocate coords of taxa and nodes
  y = zeros(n_anc,1);

  y(taxon) = [1:n_taxa];     % Get y-coords of terminal taxa

  for lev = 0:max_level      % Get y-coords of nodes, from tips to root
    for e = 1:2:n_edges-1
      if (level([e:e+1]) == [lev lev]')  % Find pair of edges at level
        a = anc(e);
        d1 = desc(e);
        d2 = desc(e+1);
        y(a) = mean([y(d1) y(d2)]);
      end;
    end;
  end;

  for e = 1:2:n_edges-1     % Get x-coords of nodes, from root to tips
    a = anc(e);
    d1 = desc(e);
    d2 = desc(e+1);
    x(d1) = x(a) + brlen(d1);
    x(d2) = x(a) + brlen(d2);
  end;

  figure;
  hold on;                    % Plot tree
    for e = 1:n_anc-1
      a = anc(e);
      d = desc(e);
      plot([x(a) x(a)],[y(a) y(d)]);
      plot([x(a) x(d)],[y(d) y(d)]);
    end;

    deltax = 0.015 * max(x);
    deltay = 0;
    for t = 1:n_taxa          % Add taxon identifiers
      text(x(t)+deltax,y(t)+deltay,num2str(t));
    end;
    axis('off');
  hold off;

  return;
