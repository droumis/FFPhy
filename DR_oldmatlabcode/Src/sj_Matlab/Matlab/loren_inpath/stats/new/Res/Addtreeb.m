% ADDTREEB: Produces and bootstraps an additive tree, given the original 
%           data matrix and the name of a function that will return a 
%           symmetric distance matrix.  
%           Bootstraps the characters across OTUs.
%
%     Syntax: [anc,brlen,dist,support] = addtreeb(X,outgrp,'func',iter)
%
%          X =        [p x n] data matrix, for which the n taxa (rows) are to be 
%                       clustered.
%          outgrp =   optional index of outgroup taxon, used only to root the tree.    
%                       If null, Farris' (1972: 658) minimum rate-heterogeneity 
%                       criterion is used to estimate the root.  See addtree().
%         'func' =    name of function (in single quotes), to be called as:
%
%                       dist = func(X)
%
%                       The function should return an [n x n] symmetric distance 
%                       matrix representing pairwise 'distances' among taxa (rows).
%          iter =     number of bootstrap iterations.
%          --------------------------------------------------------------------------
%          anc =      [1 x 2n-1] vector specifying ancestor function, with ancestor 
%                       of root specified as 0.
%          brlen =    optional corresponding vector of branch lengths.
%          dist =     [k x k] symmetric matrix of intergroup distances.
%          support =  [(k-1) x k] matrix specifying group membership and frequency 
%                       support at each node; col 1 is the freq support, cols 2:k are 
%                       clusters.
%

function [anc_orig,brlen,dist_orig,support] = addtreeb(X,outgrp,func,iter)
  [ntaxa,nchars] = size(X);
  nnodes = ntaxa - 1;

  dist_orig = feval(func,X);                      % Get distance matrix
  [anc_orig,brlen] = addtree(dist_orig,outgrp);   % Get original tree

  combtaxa = zeros(nnodes,ntaxa);                 % Allocate support-matrix components
  support = zeros(nnodes,1);
  cladesize = zeros(nnodes,1);

  for node = 1:nnodes                             % Determine clade membership for each node
    tips = treetips(anc_orig,node+ntaxa);           % Terminal taxa within clade
    ntips = length(tips);
    cladesize(node) = ntips;                        % Number of taxa within clade
    combtaxa(node,1:ntips) = tips';                 % Stash list of terminal taxa as row
  end;

  [cladesize,i] = sort(cladesize);                % Sort taxon combinations by cladesize
  combtaxa = combtaxa(i,:);
cladesize
combtaxa

  for it = 1:iter                                 % Bootstrap characters
    Xb = bootsamp(X')';                             % Bootstrapped data matrix
    dist = feval(func,Xb);                          % Distance matrix from bootstrapped data
    anc = addtree(dist,outgrp,[],1);                % Get ancestor function of bootstrapped tree

    for node = 1:nnodes                             % Cycle thru nodes of bootstrapped tree
      tips = treetips(anc,node+ntaxa)';               % Terminal taxa within clade
      ntips = length(tips);
      ni = find(cladesize==ntips);                    % Find original clade with same taxa
      for i = 1:length(ni)
        if (combtaxa(i,1:ntips)==tips)                % If exists,
          support(i) = support(i) + 1/iter;           %   increment freqency of support
        end;
      end;
    end;
  end;

  support = [support combtaxa];                   % Concatenate matrices for output

  return;

