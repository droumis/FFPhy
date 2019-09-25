% VISITADJ: Recursive function to visit (identify) the connected vertices of a 
%           graph via an adjacency matrix, beginning with a specified vertex k.
%
%     syntax: [visit,adjacency] = visitadj(k,visit,adjacency)
%
%           k -         initial vertex.
%           visit -     [n x 1] boolean vector of visited vertices
%           adjacency - boolean adjacency matrix specifying edges
%

function [visit,adjacency] = visitadj(k,visit,adjacency)
  con = find(adjacency(k,:));         % Find edges for current vertex
  lencon = length(con);

  if (lencon>0)
    visit(con) = ones(lencon,1);        % Identify connections as visited
    visit(k) = 1;
    adjacency(con,k) = zeros(lencon,1); % Remove edges from adjacency matrix
    adjacency(k,con) = zeros(1,lencon);

    for i = 1:lencon                    % Recursively follow connections
      visit = visitadj(con(i),visit,adjacency);
    end;
  end;

  return;
