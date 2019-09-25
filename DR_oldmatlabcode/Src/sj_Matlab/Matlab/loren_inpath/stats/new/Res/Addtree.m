% ADDTREE: Finds an additive tree networks from a distance matrix using the Saitou & 
%         Nei (1987) neighbor-joining method, simplified as in Studier & Kepler (1988).  
%         Uses the Kuhner & Felsenstein (1994) method to adjust negative branch lengths.  
%         The network is rooted at a specified outgroup to form an additive tree.  
%         If no outgroup is specified, Farris' (1972: 658) minimum rate-heterogeneity 
%         criterion is used to estimate the root.  Returns an ancestor function and 
%         vector of corresponding branch lengths.  
%
%     Usage: [anc,brlen,support] = addtree(D,{outgroup},{labels},{noplot | fontsize})
%
%         D =         [n x n] symmetric distance matrix for n taxa.
%         outgroup =  optional index of outgroup taxon, used only to root the tree.    
%                       If omitted or null, Farris' (1972: 658) minimum 
%                       rate-heterogeneity criterion is used to estimate the root.
%         labels =    optional matrix of labels (rows of string matrix) for terminal 
%                       taxa [default = sequence number in data matrix].
%         noplot =    optional boolean flag (0,1) indicating that tree is not to be 
%                       drawn [default = 0].
%         fontsize =  optional font size for printed labels on tree [default = 10].
%         --------------------------------------------------------------------------------
%         anc =       [1 x 2n-1] vector specifying ancestor function, with ancestor 
%                       of root specified as 0.
%         brlen =     optional corresponding vector of branch lengths.
%         support =   [(n-2) x (n-1)] matrix, with one row for all but the base 
%                       node, specifying group membership (support) at each node.
%

% If a branch length is negative, it is set to zero and the difference transferred to 
% the adjacent branch.  Moves the position of the node, but doesn't affect the total distance 
% between an adjacent pair of nodes.  But if the branch length of the adjacent branch then 
% becomes negative, it is set to zero, lengthing the distance between adjacent pair of nodes 
% and thus the total tree length.

% Farris, JS. 1972. Estimating phylogenetic trees from distance matrices.  Am. Nat.
%   106:645-668.
% Kuhner, MK & J Felsenstein. 1994. A simulation comparison of phylogeny algorithms 
%   under equal and unequal evolutionary rates.  Mol. Biol. Evol. 11:459-468.
% Saitou, N & M Nei. 1987. The neighbor-joining method: a new method for reconstructing 
%   phylogenetic trees.  Mol. Biol. Evol. 4:406-425.
% Studier, JA and KJ Kepler. 1988. A note on the neighbor-joining algorithm of Saitou 
%   and Nei.  Mol. Biol. Evol. 5:729-731.
% Swofford, DL, GJ Olsen, PJ Waddell, & DM Hillis. 1996. Phylogenetic inference.  
%   Pp 407-514 (ch 11), in DM Hillis, C Moritz & BK Mable, Molecular Systematics, 2nd ed.
%   Sinauer.
%

% RE Strauss, 8/17/98

function [anc,brlen,support] = addtree(D,outgroup,labels,fontsize)
  [r,c] = size(D);                      % Check distance matrix

  if (r~=c | sum(diag(D))>0)
    error('Addtree: distance matrix must be square symmetric, zeros on diagonal');
  else
    N = r;
  end;

  if (nargin < 2)                       % Check for other input arguments
    outgroup = [];
  end;
  if (nargin < 3)
    labels = [];
  end;
  if (nargin < 4)
    fontsize = [];
  end;

  get_support = 0;
  if (nargout > 2)
    get_support = 1;
    support = zeros(N-2,N-1);
  end;

  if (isempty(outgroup))                % Default input arguments
    outgroup = 0;
  end;

  if (isempty(labels))
    labels = tostr(1:N);
  end;

  if (isempty(fontsize))
    fontsize = 10;
    noplot = 0;
  else
    if (fontsize > 1)
      noplot = 0;
    else
      noplot = fontsize;
      fontsize = 10;
    end;
  end;

  taxon = 1:N;                            % Initial terminal-taxon labels
  
  anc = zeros(1,2*N-1);                   % Ancestor function
  brlen = anc;                            % Branch lengths

  for step = 1:(N-1)                      % Find tree
    next_node = N+step;
    n = N-step+1;
    r = sum(D);                             % Net adj divergences from all other taxa

    if (step == N-1)
      i = 1;
      j = 2;
    else
      M = zeros(n,n);                       % Find min rate-corrected distance
      for i = 1:(n-1)
        for j = (i+1):n
          M(j,i) = D(i,j) - (r(i)+r(j))/(n-2);
        end;
      end;

      [c,j,i] = trilow(M);      
      M = min(c);                           % Find taxon pair having min rate-corr dist
      k = find(c-M < eps);

      if (length(k)>1)                        % >1 min
        d = trilow(D);                          % Find min original dist
        [dmin,di] = min(d(k));
        i = i(k(di));
        j = j(k(di));
      else                                    % Unique min value
        i = i(k);
        j = j(k);
      end;
    end;

    anc(taxon(i)) = next_node;                   % Attach taxa to common node
    anc(taxon(j)) = next_node;

    d = D(i,j);                           % Branch lengths of taxa to node
    if (step == N-1)
      brlen(taxon(i)) = d/2;
      brlen(taxon(j)) = d/2;
    else
      brlen(taxon(i)) = d/2 + (r(i)-r(j))/(2*n-4);
      brlen(taxon(j)) = d - brlen(i);

      if (brlen(taxon(i))<0)                % If negative branch length,
        brlen(taxon(j)) = brlen(taxon(j))+brlen(taxon(i));    % Transfer to adjacent branch
        if (brlen(taxon(j))<0)                % If adjacent branch negative,
          brlen(taxon(j)) = 0;                %   set equal to zero (lengthing the tree)
        end;
        brlen(taxon(i)) = 0;                   % Set original negative brlen equal to zero
      elseif (brlen(taxon(j))<0)
        brlen(taxon(i)) = brlen(taxon(i))+brlen(taxon(j));   
        if (brlen(taxon(i))<0)                
          brlen(taxon(i)) = 0;                
        end;
        brlen(taxon(j)) = 0;                  
      end;

      dnew = zeros(n-2,1);
      kk = 0;
      for k = 1:n                           % Dists from node to other taxa
        if (k~=i & k~=j)
          kk = kk+1;
          dnew(kk) = 0.5*(D(i,k)+D(j,k)-d);
        end;
      end;

      D([i,j],:) = [];                      % Remove taxa from distance matrix
      D(:,[i,j]) = [];
      D = [D dnew; dnew' 0];

      taxon([i,j]) = [];
      taxon = [taxon next_node];
    end;
  end;

  if (outgroup)                                 % If outgroup taxon is specified,
    [anc,brlen] = reroot(anc,brlen,outgroup);   %   root there,
  else
    [anc,brlen] = treeroot(anc,brlen);          %   else root by min var(patristic dists)
  end;

  if (get_support)
    r = 0;
    for node = (N+1):(2*N-1)
      if (anc(node))
        tips = treetips(anc,node)';
        r = r+1;
        support(r,1:length(tips)) = tips;
      end;
    end;
  end;

  if (~noplot)
    treeplot(anc,brlen,labels,0,fontsize);
  end;

  return;

