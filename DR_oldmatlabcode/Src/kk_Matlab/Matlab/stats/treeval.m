function [id,nodes,idnames]=treeval(Tree,varargin)
%TREEVAL Compute fitted value for decision tree applied to data.
%   YFIT = TREEVAL(TREE,X) takes a classification or regression tree TREE
%   as produced by the TREEFIT function, and a matrix X of predictor
%   values, and produces a vector YFIT of predicted response values.
%   For a regression tree, YFIT(j) is the fitted response value for a
%   point having the predictor values X(j,:).  For a classification tree,
%   YFIT(j) is the class number into which the tree would assign the point
%   with data X(j,:).  To convert the number into a class name, use the
%   third output argument (see below).
%
%   YFIT = TREEVAL(TREE,X,SUBTREES) takes an additional vector SUBTREES of
%   pruning levels, with 0 representing the full, unpruned tree.  TREE must
%   include a pruning sequence as created by the TREEFIT or TREEPRUNE function.
%   If SUBTREES has K elements and X has N rows, then the output YFIT is an
%   N-by-K matrix, with the Jth column containing the fitted values produced by
%   the SUBTREES(J) subtree.  SUBTREES must be sorted in ascending order.
%   (To compute fitted values for a tree that is not part of the optimal
%   pruning sequence, first use TREEPRUNE to prune the tree.)
%
%   [YFIT,NODE] = TREEVAL(...) also returns an array NODE of the same size
%   as YFIT containing the node number assigned to each row of X.  The
%   TREEDISP function can display the node numbers for any node you select.
%
%   [YFIT,NODE,CNAME] = TREEVAL(...) is valid only for classification trees.
%   It retuns a cell array CNAME containing the predicted class names.
%
%   NaN values in the X matrix are treated as missing.  If the TREEVAL
%   function encounters a missing value when it attempts to evaluate the
%   split rule at a branch node, it cannot determine whether to proceed to
%   the left or right child node.  Instead, it sets the corresponding fitted
%   value equal to the fitted value assigned to the branch node.
%
%   Example: Find predicted classifications for Fisher's iris data.
%      load fisheriris;
%      t = treefit(meas, species);  % create decision tree
%      sfit = treeval(t,meas);      % find assigned class numbers
%      sfit = t.classname(sfit);    % get class names
%      mean(strcmp(sfit,species))   % compute proportion correctly classified
%
%   See also TREEFIT, TREEPRUNE, TREEDISP, TREETEST.

%   Copyright 1993-2006 The MathWorks, Inc. 
%   $Revision: 1.2.4.3 $  $Date: 2006/11/11 22:55:57 $
if isa(Tree,'struct')
    Tree = classregtree(Tree);
end

if isequal(Tree.method,'regression')
    [id,nodes] = eval(Tree,varargin{:});
else
    [idnames,nodes,id] = eval(Tree,varargin{:});
end