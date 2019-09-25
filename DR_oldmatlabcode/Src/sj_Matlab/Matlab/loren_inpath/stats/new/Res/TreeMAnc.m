% TREEMANC: Given the number of terminal taxa T and a value of Rohlf's topology 
%           number n integer M, recovers the tree topology as an ancestor function.
%
%     Usage: anc = treemanc(T,M)
%
%           T =   number of terminal taxa.
%           M =   Rohlf's topology number, from 0 to treenum(T)-1.
%           ------------------------------------------------------
%           anc = ancestor function of length 2T-1 (row vector).
%

% Rohlf, FJ. 1983. Numbering binary trees with labeled terminal vertices.  Bull. Math. 
%   Biol. 45:433-40.

function anc = treemanc(T,M)
  N = treevect(T,M);          % Get Rohlf's N-tuple describing topology

  n_anc = 2*T-1;              % Length of ancestor function
  anc = zeros(1,n_anc);       % Initialize ancestor function

  n_edge = 2*T-2;             % Number of edges
  edge = zeros(n_edge,2);     % Initialize edge connections

  h = T+1;                    % Current node
  edge(1,:) = [1,h];          % Initial edges
  edge(2,:) = [2,h];

  % Find end-points of edges of tree
  for t = 3:T
    e = N(t)+1;               % Insertion edge for next taxon
    e_max = 2*(t-2)+1;        % Max edge
    h = h+1;                  % Next node

    if (e < e_max)            % Insert taxon into tree
      for i = (n_edge-2):-1:(e+1) % Shift subsequent edges downward
        edge(i+2,:) = edge(i,:);
      end;
      h_save = edge(e,2);
      edge(e,2) = h;
      edge(e+1,:) = [t,h];
      edge(e+2,:) = [h,h_save];
    else % (append)           % Append taxon onto tree
      h_save = edge(e-1,2);
      edge(e,:) = [h_save,h];
      edge(e+1,:) = [t,h];
    end;
  end;

  % Get ancestor function from edges
  for e = 1:n_edge
    anc(edge(e,1)) = edge(e,2);
  end;
  return;
