function Tree=treefit(X,y,varargin)
%TREEFIT Fit a tree-based model for classification or regression.
%   T = TREEFIT(X,Y) creates a decision tree T for predicting response Y
%   as a function of predictors X.  X is an N-by-M matrix of predictor
%   values.  Y is either a vector of N response values (for regression),
%   or a character array or cell array of strings containing N class
%   names (for classification).  Either way, T is binary tree where each
%   non-terminal node is split based on the values of a column of X.  NaN
%   values in X or Y are taken to be missing values, and observations with
%   any missing values are not used in the fit.
%
%   T = TREEFIT(X,Y,'PARAM1',val1,'PARAM2',val2,...) specifies optional
%   parameter name/value pairs:
%
%   For all trees:
%      'categorical' Vector of indices of the columns of X that are to be
%                   treated as unordered categorical variables
%      'method'     Either 'classification' (default if Y is text) or
%                   'regression' (default if Y is numeric)
%      'names'      A cell array of names for the predictor variables,
%                   in the order in which they appear in the X matrix
%                   from which the tree was created (see TREEFIT)
%      'splitmin'   A number N such that impure nodes must have N or more
%                   observations to be split (default 10)
%      'prune'      'on' (default) to compute the full tree and the optimal
%                   sequence of pruned subtrees, or 'off' for the full tree
%                   without pruning
%
%   For classification trees only:
%      'cost'       Square matrix C, C(i,j) is the cost of classifying
%                   a point into class j if its true class is i (default
%                   has C(i,j)=1 if i~=j, and C(i,j)=0 if i=j).  Alternatively
%                   this value can be a structure S having two fields:  S.group
%                   containing the group names as a character array or cell
%                   array of strings, and S.cost containing the cost matrix C.
%      'splitcriterion'  Criterion for choosing a split, either 'gdi' (default)
%                   for Gini's diversity index, 'twoing' for the twoing rule,
%                   or 'deviance' for maximum deviance reduction
%      'priorprob'  Prior probabilities for each class, specified as a
%                   vector (one value for each distinct group name) or as a
%                   structure S with two fields:  S.group containing the group
%                   names as a character array or cell array of strings, and
%                   S.prob containing a a vector of corresponding probabilities
%
%   Example:  Create classification tree for Fisher's iris data.
%      load fisheriris;
%      t = treefit(meas, species, 'names',{'SL' 'SW' 'PL' 'PW'})
%      view(t);
%
%   See also TREEDISP, TREEPRUNE, TREETEST, TREEVAL.

%   Reference:  Breiman et al. (1993), "Classification and Regression
%   Trees," Chapman and Hall, Boca Raton.

%   Copyright 1993-2006 The MathWorks, Inc. 
%   $Revision: 1.6.4.11 $  $Date: 2006/12/15 19:30:20 $

Tree = classregtree(X,y,varargin{:});