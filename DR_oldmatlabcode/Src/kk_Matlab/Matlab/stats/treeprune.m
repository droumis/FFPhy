function Tree = treeprune(Tree,varargin)
%TREEPRUNE Produce a sequence of subtrees by pruning.
%   T2 = TREEPRUNE(T1,'level',LEVEL) takes a decision tree T1 as created
%   by the TREEFIT function, and a pruning level LEVEL, and returns the
%   decision tree T2 pruned to that level.  The value LEVEL=0 means
%   no pruning.  Trees are pruned based on an optimal pruning scheme
%   that first prunes branches giving less improvement in error cost.
%
%   T2 = TREEPRUNE(T1,'nodes',NODES) prunes the nodes listed in the NODES
%   vector from the tree.  Any T1 branch nodes listed in NODES become
%   leaf nodes in T2, unless their parent nodes are also pruned.  The
%   TREEDISP function can display the node numbers for any node you select.
%
%   T2 = TREEPRUNE(T1) returns the decision tree T2 that is the full,
%   unpruned T1, but with optimal pruning information added.  This is
%   useful only if you created T1 by pruning another tree, or by using
%   the TREEFIT function with pruning set 'off'.  If you plan to prune
%   a tree multiple times along the optimal pruning sequence, it is more
%   efficient to create the optimal pruning sequence first.
%
%   Pruning is the process of reducing a tree by turning some branch
%   nodes into leaf nodes, and removing the leaf nodes under the
%   original branch.
%
%   Example:  Display full tree for Fisher's iris data, as well as
%   the next largest tree from the optimal pruning sequence.
%      load fisheriris;
%      t = treefit(meas,species,'splitmin',5);
%      view(t);
%      t1 = treeprune(t,'level',1);
%      view(t1);
%
%   See also TREEFIT, TREETEST, TREEDISP, TREEVAL.

%   Copyright 1993-2006 The MathWorks, Inc. 
%   $Revision: 1.3.4.4 $  $Date: 2006/12/15 19:30:21 $

if isa(Tree,'struct')
    Tree = classregtree(Tree);
end

Tree = prune(Tree,varargin{:});