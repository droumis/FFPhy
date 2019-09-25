% <<< Incomplete >>>

% UPGMAFLX: Optimized flexible unweighted pair-group hierarchical cluster analysis 
%           of a distance matrix.  Optimizes the flexible parameter beta so as to
%           minimize the sum of squared deviations of patristic distances from
%           original distances.
%           Ref. Belbin, Faith & Milligan (1992).
%
%     Usage: [topology,support] = upgmaflx(dist,{labels},{suppress},{fontsize})
%
%         dist =     [n x n] symmetric distance matrix.
%         labels =   optional [n x q] matrix of group labels for dendrogram.
%         suppress = optional flag; if present and true (=1), suppresses
%                       graphical dendrogram output.
%         fontsize = optional font size for labels [default = 10].
%         -----------------------------------------------------------------------
%         topology = [(n-1) x 4] matrix summarizing dendrogram topology:
%                       col 1 = 1st OTU/cluster being grouped at current step
%                       col 2 = 2nd OTU/cluster
%                       col 3 = ID of cluster being produced
%                       col 4 = distance at node
%         support =  [(n-2) x (n-1)] matrix, with one row for all but the base 
%                       node, specifying group membership (support) at each node.
%

% Belbin, L, DP Faith & GW Milligan. 1992. A comparison of two approaches to beta-
%   flexible clustering. Multivariate Behavioral Research 27:417-433.

% RE Strauss, 1/3/99

function [topology,support] = upgmaflx(dist,labels,suppress,fontsize)

  beta = constr('upgmaf',beta,[],-0.9,+0.9,[],

X=CONSTR('FUN',X,OPTIONS,VLB,VUB,'GRADFUN')

  return;

