function outfig = treedisp(Tree,varargin)
%TREEDISP Show classification or regression tree graphically.
%   TREEDISP(T) takes as input a decision tree T as computed by the
%   TREEFIT function, and displays it in a figure window.  Each branch
%   in the tree is labeled with its decision rule, and each terminal node
%   is labeled with the predicted value for that node.  You can click on
%   any node to get more information about it, as specified by the
%   pop-up menu at the top of the figure.
%
%   TREEDISP(T,'PARAM1',val1,'PARAM2',val2,...) specifies optional
%   parameter name/value pairs:
%
%      'names'       A cell array of names for the predictor variables,
%                    in the order in which they appear in the X matrix
%                    from which the tree was created (see TREEFIT)
%      'prunelevel'  Initial pruning level to display
%
%   For each branch node, the left child node corresponds to the points
%   that satisfy the condition, and the right child node corresponds to
%   the points that do not satisfy the condition.
%
%   Example:  Create and graph classification tree for Fisher's iris data.
%             The names are abbreviations for the column contents (sepal
%             length, sepal width, petal length, and petal width).
%      load fisheriris;
%      t = treefit(meas, species);
%      treedisp(t,'names',{'SL' 'SW' 'PL' 'PW'});
%
%   See also TREEFIT, TREETEST, TREEPRUNE, TREEVAL.

%   Copyright 1993-2006 The MathWorks, Inc.
%   $Revision: 1.4.4.11 $  $Date: 2006/11/11 22:55:53 $

if isa(Tree,'struct')
    Tree = classregtree(Tree);
end

view(Tree,varargin{:});