% RANDUPG:  Produces random UPGMA trees from Euclidean distances of uniformly distributed 
%           character states.  Returns ancestor functions.
%
%

function [anc,i] = randupg(T,ntrees)
  
  



      Usage: [anc,brlen,links,curr_node] = lnktoanc(anc,brlen,links,curr_node)
 
            anc =       ancestor function, with ancestor of root node 
                          specified as 0; initially pass as null vector.
            brlen =     corresponding vector of branch lengths; initially pass as null 
                          vector.
            links =     3-column matrix of ancestors & descendants (not necessarily in 
                          that order) and branch lengths. 
            curr_node = current node; initially pass as root node.

  return;
