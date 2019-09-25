% TREEWALK: Produce a random tree by a random walk through a p-dimensional space 
%           via random step-lengths from a Minkowski distribution (for which k=2 
%           corresponds to a Brownian-motion model).  Random walk may be along a 
%           given tree (specified by an ancestor function and, optionally, 
%           corresponding branch lengths that are proportional to times between 
%           cladogenetic events) or produced according to expected times between 
%           cladogenetic events. 
%
%     Usage: [w_anc,w_brlen] = trewalk(P,{stepstd},{tmax},{k},{anc},{brlen})
%
%           P =       dimension of space.
%           stepstd = standard deviation of step-lengths along branches 
%                       [default = 1].
%           tmax =    number of total time-steps for tree [default = 1000].
%           k =       Minkowski parameter [default = 2].
%           anc =     ancestor function for specified tree.
%           brlen =   if 'anc' is given: corresponding branch lengths (time 
%                       durations) for specified tree;
%                     if 'anc' is not given: 2-element vector of mean and 
%                       standard deviation of expected branch lengths (in time 
%                       steps).
%           ------------------------------------------------------------------
%           w_anc =   ancestor function for resulting tree (same as anc, if 
%                       specified).
%           w_brlen = branch lengths of resulting tree.
%

% RE Strauss, 12/26/99

function [w_anc,w_brlen] = trewalk(P,stepstd,tmax,k,anc,brlen)

  
  return;
