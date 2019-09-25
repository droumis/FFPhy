% ROUTING:  Given a set of m points, finds the shortest serial (branchless) connection 
%           among n nodes of a net, which need not be complete.  The start node must be 
%           specified, and the path can terminate in a specified end node or may terminate 
%           with any node.  If the start and end nodes are identical, the path describes 
%           a Hamiltonian circuit.
%
%           Uses the heuristic method of Fencl (1972), ACM Algorithm 456.
%
%     Usage: [p,l] = routing(p,d,sn,{en},{r})
%
%         p = vector of nodes to be connected.
%         d = [m x m] distance matrix, with zeros on the diagonal.  The matrix need not 
%               be symmetric, and non-existent links in the network can be specified by
%               large positive values (>n*max(d)).
%         sn = starting node number.
%         en = ending node number [0 or null if undetermined].
%         r = number of repeated trials [max = 2n, default = n].
%         -------------------------------------------------------------------------------
%         p = optimal connection of nodes.
%         l = total length of connection.
%

% Fencl, Z. 1972. Algorithm 456, Routing Problem [H]. Collected Algorithms from CACM.

function [p,l] = routing(p,d,sn,en,r)
  n = length(p);                      % Number of nodes to be connected
  [m,c] = size(d);

  if (nargin < 4)
    en = [];
  end;
  if (nargin < 5)
    r = [];
  end;

  if (isempty(en))
    en = 0;
  end;
  if (isempty(r))
    r = n;
  end;

  large = n * max(max(d));            % Large number

  if (en == 0)                        % Define non-existing arcs by assigning their
    id = d(:,sn);                     %   distances large negative values
    d(:,sn) = -large * ones(m,1);
    d(sn,sn) = 0;
  end;

  if (sn~=en & en~=0)
    id = zeros(m,1);
    id(1) = d(en,sn);
    d(en,sn) = -large;
  end;

  l = large;
  for irs = 1:r                       % Run trials
irs

    for js = 2:n                        % Build tour by consecutive insertion of nodes
      mininc = large;
      je = js-1;                          % Last node entered

      for j = js:n                        % Find unconnected node having minimum increment
        jp = p(j);                          % Potential new node
        for i = 1:je                        % Cycle thru connected nodes
          ip1 = p(i);
          if (i < je)
            ip2 = p(i+1);
          else
            ip2 = p(1);
          end;
          inc = d(ip1,jp) + d(jp,ip2) - d(ip1,ip2);
          if (inc < mininc)
            jsave = j;
            isave = i;
            mininc = inc;
          end;
        end;
      end;  % for j = js:n
      
      i = isave;
      j = jsave-1;

      while (j ~= i)                      % Stretch tour by inserting chosen node
        ip = p(j);
        p(j) = p(j+1);
        p(j+1) = ip;
        j = j-1;
      end;
    end;  % for js = 2:n
% 90

    n1 = n-1;
    if (n > 2)                            % Adjust tour by the 3-opt method,
      for k = 1:n1                        %   varying consecutive chain length k
        icount = 0;
    
        while (icount < n)                  % Shift consecutive chain throughout sequence of n nodes
          icor = 0;
          
          for j = 1:n                         % Calc chain length in forward & backward directions
            l1 = 0;
            lr = 0;

            if (k > 1)
              i = j;
              k1 = 1;

              while (k1 < k)
                if (i > n)
                  i = i - n;
                end;
                ip = p(i);
                ip1 = i+1;
                if (ip1 > n)
                  ip1 = 1;
                end;
                ip1 = p(ip1);
                l1 = l1 + d(ip,ip1);
                lr = lr + d(ip1,ip);
                i = i+1;
                k1 = k1 + 1;
              end;
            end;

            mininc = large;                 % For each positioned chain (as is & inverted),
            j1 = j + k - 1;                 %   check all arcs if insertion improves tour
            if (j1 > n)
              j1 = j1 - n;
            end;

            for i = 1:n
              case1 = (j<=j1 & (i>=j & i<=j1));
              case2 = (j>j1 & (i<=j1 | i>=j));
              if (~case1 | ~case2)
                ip = p(i);
                jp = p(j);
                jp1 = p(j1);
                ip1 = i+1;
                if (ip1 > n)
                  ip1 = 1;
                end;
                je = ip1;
                if (ip1 == j)
                  ip1 = j1+1;
                end;
                if (ip1 > n)
                  ip1 = 1;
                end;
                ip1 = p(ip1);
                ln = l1;
                ir = 0;

% 130
                doloop = 1;
                while (doloop)
                  inc = d(ip,jp) + ln + d(jp1,ip1) - d(ip,ip1);
                  case3 = (inc>mininc | (inc==mininc & (je~=j | (je==j & ir==1))));
                  if (~case3)
                    i1 = i;
                    ir1 = ir;
                    mininc = inc;
                  end;
% 140
                  if (ir ~= 1)
                    ir = 1;
                    ln = lr;
                    js = jp;
                    jp = jp1;
                    jp1 = js;
                  else
                    doloop = 0;
                  end;
                end;
              end;
            end;
% 150
            i = i1+1;
            if (i > n)
              i = 1;
            end;

            if (i~=j | ir1~=0)
              icor = 1;                     % Reinsert chain of length k starting in j
              if (ir1 == 0)                 %   between nodes p(i1) and p(i1+1)
                js = j;
                je = 0;
              else
                js = j1;
                je = -1;
              end;
              k1 = 0;
% 160              
              doloop = 2;
              k1 = k1 + 1;

              while (doloop)
                if (k1 > k)
                  doloop = 0;
                else  
                  if (doloop==2 | ip-i1==0)
                    i = js;
                    js = js + je;
                    if (js < 1)
                      js = n;
                    end;
                    doloop = 1;
                  end;
                
                  ip = i+1;
                  if (ip > n)
                    ip = 1;
                  end;
                  jp = p(i);
                  p(i) = p(ip);
                  p(ip) = jp;
                  i = i+1;
                  if (i > n)
                    i = 1;
                  end;
  
                  if (ip-i1 == 0)
                    k1 = k1+1;
                  end;
                end;
              end;
            end;
          end;
% 190

          if (icor ~= 0)
            icount = icount+1;
          end;
        end;
      end;  % for k = 1:n1
% 200

    end;  % if (n > 2)
% 210

    for i = 1:n                             % Orient tour with sn in p(1)
      if (p(1) ~= sn)
        js = p(1);
        p(1:(n-1)) = p(2:n);
        p(n) = js;
      end;
    end;
% 230

    l1 = 0;                                 % Calculate tour length
    for i = 1:(n-1)
      ip = p(i);
      ip1 = p(i+1);
      l1 = l1 + d(ip,ip1);
    end;

    ip = p(1);
    if (sn == en)
      l1 = l1 + d(ip1,ip);
    end;
l1

    if (l1 < l)                             % Save solution, if better, and save
      l = l1;                               %   new initiate node
      q = p;
l
q
    end;

    j = irs+1;
    if (j > n)
      j = j - n;
    end;
    js = p(1);
    p(1) = p(j);
    p(j) = js;

  end;  % for irs = 1:r
% 280

  p = q;
    
  end;

  return;


