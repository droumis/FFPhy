% PATHMIN:  Given the distances among a set of n points, finds the shortest 
%           serial (branchless) connection among n points of a net, which need 
%           not be complete.  The start point must be specified, and the path 
%           can terminate in a specified end point or may terminate with any point.  
%           If the start and end points are identical, the path describes a 
%           Hamiltonian circuit.
%
%           Uses the heuristic method of Fencl (1972), ACM Algorithm 456.  If 
%           only the start point is given, uses the nearest-neighbor path as the 
%           initial estimate, otherwise constructs the initial path by the 
%           consecutive insertion of points (local optimum).
%
%     Usage: [minpath,minlen] = pathmin(dist,sn,{en},{p},{ntrials})
%
%         dist =    [n x n] distance matrix, with zeros on the diagonal.  The 
%                     matrix need not be symmetric, and non-existent links in 
%                     the network can be specified by large positive values 
%                     (>n*max(d)).
%         sn =      starting point number.
%         en =      ending point number [0 or null if undetermined].
%         p =       vector of points to be connected [default = 1:n].
%         ntrials = number of repeated trials [max = 2n, default = n].
%         --------------------------------------------------------------------
%         minpath = optimal connection of points.
%         minlen =  total length of connection.
%

% Fencl, Z. 1972. Algorithm 456, Routing Problem.
%   Collected Algorithms from CACM.

% RE Strauss, 7/7/98
%   9/27/00 - check against nearest-neighbor path.

function [minpath,minlen] = pathmin(dist,sn,en,p,ntrials)
  if (nargin < 3) en = []; end;
  if (nargin < 4) p = []; end;
  if (nargin < 5) ntrials = []; end;

  [m,c] = size(dist);
  if (m~=c | sum(diag(dist)))
    error('  PATHMIN: distance matrix must be square with zero diagonals');
  end;

  if (isempty(en))
    en = 0;
  end;
  if (isempty(p))
    p = 1:m;
  end;

  n = length(p);                      % Number of points
  if (max(p) > m)
    error('  PATHMIN: point labels exceed order of distance matrix');
  end;

    if (isempty(ntrials))
    ntrials = n;
  end;

  large = n * max(max(dist));         % Large value

  if (en == 0)                        % Connect starting pt to potential ending pts
    dist(:,sn) = -large * ones(m,1);  %   by large negative distances
    dist(sn,sn) = 0;
  elseif (sn~=en)
    dist(en,sn) = -large;
  end;

  minlen = large;
  for trial = 1:ntrials               % Run trials

    if (en==0 | en~=sn)                 % If no endpoint & no circuit,
      pp = pathnn(dist(p,p));             % Find nearest-neighbor path
      p = p(pp);
    else                                % Else
      for js = 2:n                        % Build path by consecutive insertion of points
        mininc = large;
        je = js-1;                          % Last point entered

        for j = js:n                        % Find unconnected point having minimum increment
          jp = p(j);                          % Potential new point
          for i = 1:je                        % Cycle thru connected points
            ip1 = p(i);
            if (i < je)
              ip2 = p(i+1);
            else
              ip2 = p(1);
            end;
            inc = dist(ip1,jp) + dist(jp,ip2) - dist(ip1,ip2);
            if (inc < mininc)
              jsave = j;
              isave = i;
              mininc = inc;
            end;
          end;
        end;  % for j = js:n
        
        i = isave;
        j = jsave-1;
  
        while (j ~= i)                      % Stretch path by inserting chosen point
          ip1 = p(j);
          p(j) = p(j+1);
          p(j+1) = ip1;
          j = j-1;
        end;
      end;  % for js = 2:n
    end;

    if (n > 2)                            % Adjust path by the 3-opt method,
      for k = 1:n-1                       %   varying consecutive chain length k
        icount = 0;
    
        while (icount < n)                  % Shift consecutive chain throughout sequence of n points
          icor = 0;
          
          for j = 1:n                         % Calc chain length in forward & backward directions
            lenf = 0;
            lenr = 0;

            if (k > 1)
              i = j;
              k1 = 1;

              while (k1 < k)
                if (i > n)
                  i = i - n;
                end;
                ip1 = p(i);
                ip2 = i+1;
                if (ip2 > n)
                  ip2 = 1;
                end;
                ip2 = p(ip2);
                lenf = lenf + dist(ip1,ip2);
                lenr = lenr + dist(ip2,ip1);
                i = i+1;
                k1 = k1 + 1;
              end;
            end;

            mininc = large;                 % For each positioned chain (as is & inverted),
            j1 = j + k - 1;                 %   check all arcs if insertion improves path
            if (j1 > n)
              j1 = j1 - n;
            end;

            for i = 1:n
              case1 = (j<=j1 & (i>=j & i<=j1));
              case2 = (j>j1 & (i<=j1 | i>=j));

              if (~case1 & ~case2)
                ip1 = p(i);
                jp = p(j);
                jp1 = p(j1);
                ip2 = i+1;
                if (ip2 > n)
                  ip2 = 1;
                end;
                je = ip2;
                if (ip2 == j)
                  ip2 = j1+1;
                end;
                if (ip2 > n)
                  ip2 = 1;
                end;
                ip2 = p(ip2);
                ln = lenf;
                ir = 0;

                doloop = 1;
                while (doloop)
                  inc = dist(ip1,jp) + ln + dist(jp1,ip2) - dist(ip1,ip2);
                  case3 = (inc>mininc | (inc==mininc & (je~=j | (je==j & ir==1))));

                  if (~case3)
                    i1 = i;
                    ir1 = ir;
                    mininc = inc;
                  end;

                  if (ir ~= 1)
                    ir = 1;
                    ln = lenr;
                    js = jp;
                    jp = jp1;
                    jp1 = js;
                  else
                    doloop = 0;
                  end;
                end;
              end;
            end;  % for i = 1:n

            i = i1+1;
            if (i > n)
              i = 1;
            end;

            if (i~=j | ir1~=0)
              icor = 1;                     % Reinsert chain of length k starting in j
              if (ir1 == 0)                 %   between points p(i1) and p(i1+1)
                js = j;
                je = 0;
              else
                js = j1;
                je = -1;
              end;
              k1 = 0;

              doloop = 2;
              k1 = k1 + 1;

              while (doloop)
                if (k1 > k)
                  doloop = 0;
                else  
                  if (doloop==2 | ip1-i1==0)
                    i = js;
                    js = js + je;
                    if (js < 1)
                      js = n;
                    end;
                    doloop = 1;
                  end;
                
                  ip1 = i+1;
                  if (ip1 > n)
                    ip1 = 1;
                  end;
                  jp = p(i);
                  p(i) = p(ip1);
                  p(ip1) = jp;
                  i = i+1;
                  if (i > n)
                    i = 1;
                  end;
  
                  if (ip1-i1 == 0)
                    k1 = k1+1;
                  end;
                end;
              end;
            end;
          end;

          if (icor == 0)
            icount = n;
          else
            icount = icount+1;
          end;

        end;  % while count < n
      end;  % for k = 1:n-1

    end;  % if (n > 2)

    for i = 1:n                             % Orient path with sn in p(1)
      if (p(1) ~= sn)
        js = p(1);
        p(1:(n-1)) = p(2:n);
        p(n) = js;
      end;
    end;

    len = 0;                                 % Calculate path length
    for i = 1:(n-1)
      ip1 = p(i);
      ip2 = p(i+1);
      len = len + dist(ip1,ip2);
    end;

    ip1 = p(1);
    if (sn == en)
      len = len + dist(ip2,ip1);
    end;

    if (len < minlen)                       % Save solution, if better, and save
      minlen = len;                         %   new initiate point
      q = p;
    end;

    t = trial+1;                            % Put new point at the lead
    if (t > n)
      t = t - n;
    end;
    s = p(1);
    p(1) = p(t);
    p(t) = s;

  end;  % for trial = 1:ntrials

  minpath = q;
  if (sn == en)
    if (size(minpath,1)==1)
      minpath = [minpath minpath(1)];
    else
      minpath = [minpath; minpath(1)];
    end;
  end;

  return;


