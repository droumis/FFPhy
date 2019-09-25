% ClusterBinary: Produces and randomizes a UPGMA or neighbor-joining tree, given the 
%                original binary data matrix (obs x vars) and the name of a function that 
%                will return a symmetric distance matrix.  Identifies the clusters that
%                are statistically significant.
%
%     Syntax: [dist,topo,support,tot_support] = 
%                 clusterbinary(X,'distfunc',{iter},{treetype},{outgroup},{bootobs}, ...
%                         {p1},{p2},{p3},{p4},{p5})
%
%          X =       [p x n] data matrix, for which the n taxa (rows) are to be 
%                      clustered.
%          'distfunc' =  name of function (in single quotes), to be called as:
%
%                      dist = func(X)
%
%                      The function should return an [n x n] symmetric 
%                      distance matrix representing pairwise 
%                      'distances' among rows (observations).
%          iter =     optional number of bootstrap iterations [default=0].
%          treetype = flag indicating the kind of cluster analysis to be done:
%                       0 = upgma [default];
%                       1 = neighbor-joining.
%          outgroup = optional outgroup for neighbor-joining tree 
%                       [default = null].
%          bootobs =  optional boolean flag indicating, if true, that observations 
%                       are to be bootstrapped (within groups) rather than 
%                       variables; requires that a grouping vector be passed as 
%                       argument 'p1' [default = 0].
%          p1-p5 =    up to five additional parameters to be passed to 
%                       'distfunc'; assumes that X is the first input argument.
%          -------------------------------------------------------------------------
%          dist =    [k x k] symmetric matrix of intergroup distances.
%          topo =    [(k-1) x 4] matrix summarizing dendrogram topology:
%                      col 1 = 1st OTU/cluster being grouped at current step
%                      col 2 = 2nd OTU/cluster
%                      col 3 = ID of cluster being produced
%                      col 4 = distance at node
%          support = [(k-2) x k] matrix, with one row for all but the base
%                       node, specifying group membership and frequency support at 
%                       each node; col 1 is the freq support, cols 2:k are clusters.
%                     If iter=0, a [(k-2) x (k-1)] matrix is returned, where 
%                       cols 1:(k-1) are clusters.
%          tot_support = matrix similar to 'support' but listing all clusters 
%                          encountered during bootstrapping.
%

% RE Strauss, 8/17/98
%   11/14/00 -  correct bootstrap sampling to sample across variables (cols) 
%                 rather than rows (obs).
%   1/21/00 -   fixed problem with bootstrapping neighbor-joining tree; 
%                 added optional outgroup.
%   5/7/01 -    corrected problem with passing additional parameters.

function [distorig,topo,support,tot_supt] = ...
                     clusterbinary(X,func,iter,treetype,outgroup,bootobs,p1,p2,p3,p4,p5)

  if (nargin < 3) iter = []; end;
  if (nargin < 4) treetype = []; end;
  if (nargin < 5) outgroup = []; end;
  if (nargin < 6) bootobs = []; end;

  p1_arg = 7;                         % Position of 'p1' in argument list

  get_tot_supt = 0;                   % Set flag for total-support matrix
  if (nargout > 3)
    get_tot_supt = 1;
  end;

  if (isempty(iter))                  % Default input parameter values
    iter = 0;
  end;
  if (isempty(treetype))
    treetype = 0;
  end;
  if (isempty(bootobs))
    bootobs = 0;
  end;

  if (bootobs)
    if (~exist('p1'))
      error('  ClusterBinary: grouping vector required to bootstrap observations.');
    end;
  end;

  [p,n] = size(X);

  evalstrX = [func, '(X'];
  evalstrA = [func, '(A'];
  for i=1:(nargin-p1_arg+1)
    evalstrX = [evalstrX,',p',int2str(i)];
    evalstrA = [evalstrA,',p',int2str(i)];
  end;
  evalstrX = [evalstrX,')'];
  evalstrA = [evalstrA,')'];

  distorig = eval(evalstrX);            % Get distance matrix
  switch (treetype)                     % Get original dendrogram
    case 0,
      [topo,support] = upgma(distorig,[],1);            % UPGMA
    case 1,
      [anc,brlen,support] = addtree(distorig,outgroup,[],1);   % Neighbor-joining
      topo = [];
    otherwise
      error('  ClusterBinary: invalid tree type.');
  end;

  if (~iter)                            % If no bootstrap, quit
    tot_supt = [];
    return;
  end;

  % Bootstrap characters and cluster analyses

  support = [ones(size(support,1),1) support];  % Attach col 1,
  tot_supt = support;                           %   = freqs of cluster support
  [sr,sc] = size(tot_supt);           % Sizes of support matrices 
  cr = sr;
  iter = iter+1;                      % Include initial solution

  for b = 2:iter                      % Bootstrap iterations
    if (bootobs)
      A = bootsamp(X,p1);               % Bootstrap observations
    else
      A = bootsamp(X')';                % Bootstrap variables
    end;
    dist = eval(evalstrA);            % Get distance matrix

    if (treetype==0)                  % Get support matrix
      [t,s] = upgma(dist,[],1);         
    elseif (treetype==1)
      [a,b,s] = addtree(dist,1,[],1);     % Root additive tree at taxon 1
    end;

    for i = 1:cr                      % Update cluster-support information
      for j = 1:cr                        % Update cluster if in support
        if (all(s(i,:)==support(j,2:sc)))
          support(j,1) = support(j,1) + 1;
        end;
      end;

      if (get_tot_supt)
        found = 0;
        for j = 1:sr                        % Check if cluster is in total-support matrix
          if (all(s(i,:)==tot_supt(j,2:sc)))
            tot_supt(j,1) = tot_supt(j,1) + 1;
            found =  1;
          end;
        end;
        if (~found)                         % If not, append
          tot_supt = [tot_supt; 1 s(i,:)];
          sr = sr + 1;
        end;
      end;
    end;
  end;  % Bootstrap

  support(:,1) = support(:,1) ./ iter; % Convert counts to proportions

  [Y,i] = sort(-support(:,1));          % Sort cluster-support freqs in
  support = support(i,:);               %   descending sequence

  if (get_tot_supt)
    tot_supt(:,1) = tot_supt(:,1) ./ iter;
    [Y,i] = sort(-tot_supt(:,1));
    tot_supt = tot_supt(i,:);
  end;

  return;




