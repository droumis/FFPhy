% TREERF: Objective function for treeroot().  Given ancestor and branch-length vectors 
%         representing a tree and information about the root, reroots the tree and 
%         finds the variance of patristic distances of terminal taxa to the root.
%
%     Usage: varpd = treerf(p,outgroup,anc,brlen)
%
%           p =         value indicating the proportion of the distance from the 
%                         outgroup taxon to its ancestor at which the root is to be 
%                         positioned.
%           outgroup =  index of outgroup taxon (position of outgroup in 'anc').
%           anc =       vector specifying ancestor function, with ancestor of 
%                         root specified as 0.
%           brlen =     corresponding vector of branch lengths.
%           ------------------------------------------------------------------------
%           varpd =     variance of patristic distances from terminal taxa to root.
%

% RE Strauss, 6/17/98
%   10/27/02 - change call from fmin() to fminbnd().

function varpd = treerf(p,outgroup,anc,brlen)
  ntaxa = (length(anc)+1)/2;                      % Number of terminal taxa
  [ranc,rbrlen] = reroot(anc,brlen,outgroup,p);   % Reroot tree
    
  v = zeros(ntaxa,1);
  for t = 1:ntaxa                       % Get variance of patristic dists to terminal taxa
    c = t;
    patdist = 0;
    while (anc(c) > 0)
      patdist = patdist + rbrlen(c);
      c = anc(c);
    end;
    v(t) = patdist;
  end;

  varpd = var(v);

  return;

