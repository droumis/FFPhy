function [idname,nodes,id]=eval(Tree,X,subtrees)
%EVAL Compute fitted value for decision tree applied to data.
%   YFIT = EVAL(T,X) takes a classification or regression tree T and a
%   matrix X of predictor values, and produces a vector YFIT of predicted
%   response values. For a regression tree, YFIT(J) is the fitted response
%   value for a point having the predictor values X(J,:).  For a
%   classification tree, YFIT(J) is the class into which the tree would
%   assign the point with data X(J,:).
%
%   YFIT = EVAL(T,X,SUBTREES) takes an additional vector SUBTREES of
%   pruning levels, with 0 representing the full, unpruned tree.  T must
%   include a pruning sequence as created by the CLASSREGTREE constructor or
%   the PRUNE method. If SUBTREES has K elements and X has N rows, then the
%   output YFIT is an N-by-K matrix, with the Ith column containing the
%   fitted values produced by the SUBTREES(I) subtree.  SUBTREES must be
%   sorted in ascending order. (To compute fitted values for a tree that is
%   not part of the optimal pruning sequence, first use PRUNE to prune the
%   tree.)
%
%   [YFIT,NODE] = EVAL(...) also returns an array NODE of the same size
%   as YFIT containing the node number assigned to each row of X.  The
%   VIEW method can display the node numbers for any node you select.
%
%   [YFIT,NODE,CNUM] = EVAL(...) is valid only for classification trees.
%   It returns a vector CNUM containing the predicted class numbers.
%
%   NaN values in the X matrix are treated as missing.  If the EVAL method
%   encounters a missing value when it attempts to evaluate the split rule
%   at a branch node, it cannot determine whether to proceed to the left or
%   right child node.  Instead, it sets the corresponding fitted value
%   equal to the fitted value assigned to the branch node.
%
%   For a tree T, the syntax [...]=T(X) or [...]=T(X,SUBTREES) also invokes
%   the EVAL method.
%
%   Example: Find predicted classifications for Fisher's iris data.
%      load fisheriris;
%      t = classregtree(meas, species);  % create decision tree
%      sfit = eval(t,meas);              % find assigned class names
%      mean(strcmp(sfit,species))        % proportion correctly classified
%
%   See also CLASSREGTREE, CLASSREGTREE/PRUNE, CLASSREGTREE/VIEW, CLASSREGTREE/TEST.

%   Copyright 2006-2007 The MathWorks, Inc. 
%   $Revision: 1.1.6.1.2.1 $  $Date: 2007/01/30 02:22:33 $

[nr,nc] = size(X);
if nc~=Tree.npred
   error('stats:treeval:BadInput',...
         'The X matrix must have %d columns.',Tree.npred);
end

if nargin<3
   subtrees = 0;
elseif prod(size(subtrees))>length(subtrees)
   error('stats:treeval:BadSubtrees','SUBTREES must be a vector.');
elseif any(diff(subtrees)<0)
   error('stats:treeval:BadSubtrees','SUBTREES must be sorted.');
end

if ~isempty(Tree.prunelist)
   prunelist = Tree.prunelist;
elseif ~isequal(subtrees,0)
   error('stats:treeval:NoPruningInfo',...
         'The decision tree TREE does not include a pruning sequence.')
else
   prunelist = repmat(Inf,size(Tree.node));
end

ntrees = length(subtrees);
nodes = doapply(Tree,X,1:nr,1,zeros(nr,ntrees),subtrees,prunelist,ntrees);
id = Tree.class(nodes);

if isequal(Tree.method,'classification')
   idname = Tree.classname(id);
else
   idname = id;
end


%------------------------------------------------
function nodes = doapply(Tree,X,rows,thisnode,nodes,subtrees,prunelist,endcol)
%DOAPPLY Apply classification rule to specified rows starting at a node.
%   This is a recursive function.  Starts at top node, then recurses over
%   child nodes.  THISNODE is the current node at each step.
%
%   NODES has one row per observation and one column per subtree.
%
%   X, NODES, PRUNELIST, and SUBTREES are the same in each recursive call
%   as they were in the top-level call.  ROWS describes the subset of X and
%   NODES to consider.  1:ENDCOL are colums of NODES and the elements of
%   SUBTREES to consider.

splitvar      = Tree.var(thisnode);
cutoff        = Tree.cut(thisnode);
assignedclass = Tree.class(thisnode);
kids          = Tree.children(thisnode,:);
catsplit      = Tree.catsplit;
prunelevel    = prunelist(thisnode);
ntrees        = size(nodes,2);   % number of trees in sequence

% For how many of the remaining trees is this a terminal node?
if splitvar==0      % all, if it's terminal on the unpruned tree
   ncols = endcol;
else                % some, if it's terminal only after pruning
   ncols = sum(subtrees(1:endcol) >= prunelevel);
end
if ncols>0          % for those trees, assign the node level now
   nodes(rows,(endcol-ncols+1:endcol)) = thisnode;
   endcol = endcol - ncols;
end

% Now deal with non-terminal nodes
if endcol > 0
   % Determine if this point goes left, goes right, or stays here
   x = X(rows,abs(splitvar));   
   if splitvar>0                % continuous variable
      isleft = (x < cutoff);
      isright = ~isleft;
      ismissing = isnan(x);
   else                         % categorical variable
      isleft = ismember(x,catsplit{cutoff,1});
      isright = ismember(x,catsplit{cutoff,2});
      ismissing = ~(isleft | isright);
   end

   subrows = rows(isleft & ~ismissing);  % left child node
   if ~isempty(subrows)
      nodes = doapply(Tree,X,subrows,kids(1),nodes,subtrees,prunelist,endcol);
   end
   
   subrows = rows(isright & ~ismissing); % right child node
   if ~isempty(subrows)
      nodes = doapply(Tree,X,subrows,kids(2),nodes,subtrees,prunelist,endcol);
   end

   subrows = rows(ismissing);            % missing, treat as leaf
   if ~isempty(subrows)
      nodes(subrows,1:endcol) = thisnode;
   end
end

   
