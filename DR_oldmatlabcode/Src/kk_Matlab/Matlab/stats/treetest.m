function [varargout] = treetest(Tree,varargin)
%TREETEST Compute error rate for tree.
%   COST = TREETEST(T,'resubstitution') computes the cost of the tree T
%   using a resubstitution method.  T is a decision tree as created by
%   the TREEFIT function.  The cost of the tree is the sum over all
%   terminal nodes of the estimated probability of that node times the
%   node's cost.  If T is a classification tree, the cost of a node is
%   the sum of the misclassification costs of the observations in
%   that node.  If T is a regression tree, the cost of a node is the
%   average squared error over the observations in that node.  COST is
%   a vector of cost values for each subtree in the optimal pruning
%   sequence for T.  The resubstitution cost is based on the same
%   sample that was used to create the original tree, so it under-
%   estimates the likely cost of applying the tree to new data.
%
%   COST = TREETEST(T,'test',X,Y) uses the predictor matrix X and
%   response Y as a test sample, applies the decision tree T to that
%   sample, and returns a vector COST of cost values computed for the
%   test sample.  X and Y should not be the same as the learning sample,
%   which is the sample that was used to fit the tree T.
%
%   COST = TREETEST(T,'crossvalidate',X,Y) uses 10-fold cross-validation to
%   compute the cost vector.  X and Y should be the learning sample, which
%   is the sample that was  used to fit the tree T.  The function
%   partitions the sample into 10 subsamples, chosen randomly but with
%   roughly equal size.  For classification trees the subsamples also have
%   roughly the same class proportions.  For each subsample, TREETEST fits
%   a tree to the remaining data and uses it to predict the subsample.  It
%   pools the information from all subsamples to compute the cost for the
%   whole sample.
%
%   [COST,SECOST,NTNODES,BESTLEVEL] = TREETEST(...) also returns the vector
%   SECOST containing the standard error of each COST value, the vector
%   NTNODES containing number of terminal nodes for each subtree, and the
%   scalar BESTLEVEL containing the estimated best level of pruning.
%   BESTLEVEL=0 means no pruning (i.e. the full unpruned tree).  The best
%   level is the one that produces the smallest tree that is within one
%   standard error of the minimum-cost subtree.
%
%   [...] = TREETEST(...,'PARAM1',val1,'PARAM2',val2,...) specifies
%   optional parameter name/value pairs chosen from the following:
%
%      'nsamples'   The number of cross-validation samples (default 10)
%      'treesize'   Either 'se' (the default) to choose the smallest
%                   tree whose cost is within one standard error of the
%                   minimum cost, or 'min' to choose the minimal cost tree
%                   (not meaningful for resubstitution error calculations)
%
%   Example:  Find best tree for Fisher's iris data using cross-validation.
%             The solid line shows the estimated cost for each tree size,
%             the dashed line marks 1 standard error above the minimum,
%             and the square marks the smallest tree under the dashed line.
%      % Start with a large tree
%      load fisheriris;
%      t = treefit(meas,species','splitmin',5);
%
%      % Find the minimum-cost tree
%      [c,s,n,best] = treetest(t,'cross',meas,species);
%      tmin = treeprune(t,'level',best);
%
%      % Plot smallest tree within 1 std. error of minimum cost tree
%      [mincost,minloc] = min(c);
%      plot(n,c,'b-o', n(best+1),c(best+1),'bs',...
%           n,(mincost+s(minloc))*ones(size(n)),'k--');
%      xlabel('Tree size (number of terminal nodes)')
%      ylabel('Cost')
%
%   See also TREEFIT, TREEDISP, TREEPRUNE, TREEVAL.

%   Copyright 1993-2006 The MathWorks, Inc. 
%   $Revision: 1.5.4.5 $  $Date: 2006/11/11 22:55:56 $

if nargin==0
    error('stats:treetest:TooFewInputs',...
          'At least one input required.')
end

if ~isa(Tree,'classregtree')
    if isa(Tree,'struct')
        Tree = classregtree(Tree);
    else
        error('stats:treetest:BadTree',...
              'First argument must be a decision tree.');
    end
end

[varargout{1:max(1,nargout)}] = test(Tree,varargin{:});