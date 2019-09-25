% TREEROOT: Given an ancestor function and corresponding vector of branch lengths, 
%           re-roots the tree so as to minimize the variance among the patristic 
%           distances from each taxon to the root.  (Solution might not be unique for 
%           discrete branch lengths.)  This variance provides a lower bound for the 
%           actual heterogeneities in evolutionary rate among lineages.
%           See Farris (1972: 658).
%
%     Usage: [new_anc,new_brlen,minvar] = treeroot(anc,brlen)
%
%           anc =       vector specifying ancestor function, with ancestor of root 
%                         specified as 0.
%           brlen =     corresponding vector of branch lengths.
%           -----------------------------------------------------------------------
%           new_anc =   ancestor function of rerooted tree.
%           new_brlen = corresponding vector of branch lengths.
%           minvar =    variance of patristic distances from terminal taxa to root.
%

% Farris, JS. 1972. Estimating phylogenetic trees from distance matrices.  Am. Nat.
%   106:645-668.

% RE Strauss, 6/17/98
%   10/27/02 - change call from fmin() to fminbnd().

function [new_anc,new_brlen,minvar] = treeroot(anc,brlen)
  ntaxa = (length(anc)+1)/2;
  root = find(anc==0);                    % Find current root
  [links,anc,curr_node,level] = treedivd([],anc,root,0);

  ancestors = links(:,1);
  descendants = links(:,2);

  pvals = [0, 0.5, 1-1e-6]';
  minv = zeros(length(descendants)*length(pvals),3);
  vh = 0;
  v = zeros(ntaxa,1);

  for i = 1:length(descendants)           % Loop thru all descendants, rerooting tree at each
    d = descendants(i);
    for j = 1:length(pvals)                 % For each internode, reroot at both ends
      p = pvals(j);
      varpd = treerf(p,d,anc,brlen);        % Reroot tree and get variance of patristic dists
      vh = vh+1;
      minv(vh,:) = [d p varpd];             % Stash in results matrix
    end;
  end;

  [minvar,i] = min(minv(:,3));               % Outgroup node having least variance
  newroot = minv(i,1);
  popt = fminbnd('treerf',0,1,[],newroot,anc,brlen);

  minvar = treerf(popt,newroot,anc,brlen);
  [new_anc,new_brlen] = reroot(anc,brlen,newroot,popt);

  return;
