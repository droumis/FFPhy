% FRIEDRAF: Applies the Friedman-Rafesky (1979) distribution-free runs test
%           for the inequality of two or more multivariate distributions.
%           Randomized significance level is evaluated by permuting group labels.
%
%     Syntax: [prob,runs] = friedraf(X,grps,iter)
%
%           X -    [N x p] data matrix.
%           grps - [N x 1] vector of group labels.
%           iter - number of permutation iterations [default = 5000].
%           -----------------------------------------------------------
%           prob - significance level.
%           runs - observed number of multivariate runs.
%

% Good, PG. 1994. Permutation Tests, p. 71-73. Springer Series in Statistics,
%   Springer-Verlag.

% RE Strauss, 7/4/96
%   11/29/99 - changed calling sequence.

function [prob,runs] = friedraf(X,grps,iter)
  if (nargin < 3) iter = []; end;

  if (isempty(iter))
    iter = 5000;
  end;

  N = length(grps);

  edges = mstree(X,1);                % MST for combined sample

  grp1 = grps(edges(:,1));            % Find grp labelss for endpoints
  grp2 = grps(edges(:,2));
  g = (grp1==grp2);                   % Find edges connecting points of same grp
  e = edges(g,:);                     % Reduce edge list
  lene = size(e,1);

  if (lene == 0)                      % If no edges,
    runs = N;                           % Runs = number of points
  
  else                                % If have edges,
    pt = zeros(N,1);                    % Boolean list of points in edge list
    pt(e(:)) = ones(2*lene,1);

    e = [e; fliplr(e)];                 % Create sparse adjacency matrix
    adj = sparse(e(:,1),e(:,2),ones(2*lene,1));

    m = 1;
    runs = sum(pt==0);                  % Include unclustered points as runs

    while (m > 0)                       % Determine number of 'runs' of connected pts
      [m,j] = max(pt);                    % Find next unexamined pt
      if (m > 0)
        pt(j) = 0;                        % Remove pt from list
        if (sum(adj(j,:)) > 0)            % If connected to at least one other pt,
          runs = runs + 1;                  % Increment number of runs
          visit = zeros(N,1);               % Boolean list of pts to be visited
          visit(j) = 1;
          visit = visitadj(j,visit,adj);    % Visit points in same run
          i = find(visit);                  % Identify points visited
          pt(i) = zeros(length(i),1);       % Remove points from list
        end;
      end;
    end;
  end;

  prob = 0;
  incr = 1/iter;

  for it = 1:iter                     % Randomized significance level
    grp = grps(randperm(N));            % Randomly permute labels

    grp1 = grp(edges(:,1));            % Find grp labelss for endpoints
    grp2 = grp(edges(:,2));
    g = (grp1==grp2);                   % Find edges connecting points of same grp
    e = edges(g,:);                     % Reduce edge list
    lene = size(e,1);

    if (lene == 0)                      % If no edges,
      r = N;                              % Runs = number of points
    else                                % If have edges,
      pt = zeros(N,1);                    % Boolean list of points in edge list
      pt(e(:)) = ones(2*lene,1);

      e = [e; fliplr(e)];                 % Create sparse adjacency matrix
      adj = sparse(e(:,1),e(:,2),ones(size(e,1),1));

      m = 1;
      r = sum(pt==0);                     % Include unclustered points as runs

      while (m > 0)                       % Determine number of 'runs' of connected pts
        [m,j] = max(pt);                    % Find next unexamined pt
        if (m > 0)
          pt(j) = 0;                        % Remove pt from list
          if (sum(adj(j,:)) > 0)            % If connected to at least one other pt,
            r = r + 1;                      % Increment number of runs
            visit = zeros(N,1);               % Boolean list of pts to be visited
            visit(j) = 1;
            visit = visitadj(j,visit,adj);    % Visit points in same run
            i = find(visit);                  % Identify points visited
            pt(i) = zeros(length(i),1);       % Remove points from list
          end;
        end;
      end;
    end;

    if (r >= runs)                      % Increment tail probability
      prob = prob + incr;
    end;
  end;

  return;

