% TOPOANC: Converts a dendrogram-type topology matrix to a cladogram-type
%           ancestor-function vector and associated branch-length vector
%
%     Usage: [anc,brlen] = topoanc(topology)
%
%         topology = [(n-1) x 4] matrix summarizing dendrogram topology:
%                       col 1 = 1st OTU/cluster being grouped at current step
%                       col 2 = 2nd OTU/cluster
%                       col 3 = ID of cluster being produced
%                       col 4 = distance at node
%         anc =      [2n-1] ancestor function
%         brlen =    corresponding lengths of branches (edges) connecting
%                      descendants with ancestors
%

% RE Strauss, 7/14/95

function [anc,brlen] = topoanc(topology)
  n = size(topology,1)+1;             % Number of terminal OTUs
  anc = zeros(1,2*n-1);               % Allocate return matrices
  brlen = anc;

  for i = 1:(n-1)
    d1 = topology(i,1);
    d2 = topology(i,2);
    a =  topology(i,3);
    nodedist = topology(i,4);

    anc(d1) = a;
    anc(d2) = a;

    prevdist = 0;
    if (d1>n)                         % If current descendant was a previous
      k = find(topology(:,3)==d1);    %   ancestor, get it's previous node distance
      prevdist = topology(k,4);
    end;
    brlen(d1) = nodedist - prevdist;

    prevdist = 0;
    if (d2>n)
      k = find(topology(:,3)==d2);
      prevdist = topology(k,4);
    end;
    brlen(d2) = nodedist - prevdist;
  end;

  return;
