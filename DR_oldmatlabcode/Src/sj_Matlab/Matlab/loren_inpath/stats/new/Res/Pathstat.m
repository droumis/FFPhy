% PATHSTAT: Calculates some basic statistics of a search path, given
%           the point coordinates and path sequence.
%           If a step length on the path is less than threshold, the actual 
%           step length is returned in the 'steplen' vector and it is added to 
%           the total path length, but it is excluded from other statistics.
%
%     Usage: [stats,npts,steplen] = pathstat(crds,thresh)
%
%           crds =   [n x 2] matrix of point coordinates.
%           thresh = threshold length of step lengths, as a proportion of mean 
%                      step length [default = mean step length / 3].
%           ------------------------------------------------------------------
%           stats =  column vector of descriptive statistics:
%                        1) total path length
%                        2) minimum step length
%                        3) maximum step length
%                        4) mean step length
%                        5) std of step length
%                        6) autocorrelation of step length
%                        7) number of path crossings
%                        8) std of path distances between crossings
%                        9) angular vector autocorrelation
%                       10) minimum absolute angular deviation
%                       11) maximum absolute angular deviation
%                       12) mean of the absolute angular deviations
%                       13) std of the absolute angular deviations
%                       14) mean of the signed angular deviations
%                       15) std of the signed angular deviations
%           npts =    number of effective points in path (more distant than 
%                       threshold).
%           steplen = vector of step lengths.
%        

% RE Strauss, 3/25/98   
%   12/16/99 - produce NaN's if insufficient points in path.
%   12/16/99 - increased number of output statistics from 10 to 14.
%    5/26/00 - exit immediately if path too short; corrected some criteria for 
%                insufficient points in path; added std of AAD.
%    6/14/00 - corrected path-crossing statistics.

function [stats,npts,step] = pathstat(crds,thresh)
  if (nargin < 2) thresh = []; end;

  [n,p] = size(crds);

  if (p~=2)
    error('  PATHSTAT: two-dimensional coordinates only.');
  end;

  stats = NaN*ones(15,1);
  path = (1:n)';
  step = [];

  if (n < 2)
    return;
  end;

  % Step lengths

  step = steplength(crds);              % Calc step lengths
  totlen = sum(step);

  if (isempty(thresh))                  % Threshold
    thresh = totlen./(3*n);
  else
    thresh = thresh * totlen./n;
  end;

  i = find(step > thresh);              % Omit short step lengths
  steplen = step(i);
  crds = crds([1;i+1],:);
  n = length(steplen)+1;
  npts = n;

  autocorr = NaN;
  if (n > 5)
    s1 = steplen(1:(n-2));              % Autocorrelation of step lengths
    s2 = steplen(2:(n-1));
    autocorr = corrcoef([s1,s2]);
    autocorr = autocorr(1,2);
  end;

  stats(1) = totlen;
  stats(2) = min(steplen);
  stats(3) = max(steplen);
  stats(4) = mean(steplen);
  stats(5) = std(steplen);
  stats(6) = autocorr;

  % Path crossings

  if (n > 3)
    crosslist = [];
    for i = 1:(n-3)                     % Cycle forward thru path segments
      segi = [crds(i,:), crds(i+1,:)];
      segs = [crds(i+2:n-1,:), crds(i+3:n,:)];
      [intr,x,y] = intrsect(segi,segs);   % Find intersections with remaining segments

      if (any(intr))                      % If any intersections,  
        seg = find(intr)+(i+1);
        for k = 1:length(seg)
          j = seg(k);                       % Get segment having intersection
          xc = x(seg(k)-(i+1));             % Get intersection coordinates
          yc = y(seg(k)-(i+1));
                                          % Find distance along segment to intersection
          deltaxi = abs(crds(i+1,1)-crds(i,1));
          deltayi = abs(crds(i+1,2)-crds(i,2));
          deltaxj = abs(crds(j+1,1)-crds(j,1));
          deltayj = abs(crds(j+1,2)-crds(j,2));

          if (deltaxi > deltayi)
            crosslist = [crosslist; ...
                         i+(abs(xc-crds(i,1))./deltaxi)];
          else
            crosslist = [crosslist; ...
                         i+(abs(yc-crds(i,2))./deltayi)];
          end;

          if (deltaxj > deltayj)
            crosslist = [crosslist; ...
                         j+(abs(xc-crds(j,1))./deltaxj)];
          else
            crosslist = [crosslist; ...
                         j+(abs(yc-crds(j,2))./deltayj)];
          end;
        end;
      end;
    end;

    lenlist = length(crosslist);        % Check for crosspoints
    intercross = [];
    if (lenlist > 0)                    % Calc inter-crosspoint distances
      crosslist = sort(crosslist);        % Sort list of intersections
      crosslist = [1; crosslist; n-(1e-6)];      % Add first and last points

      intercross = zeros(lenlist+1,1);
      for i = 1:(lenlist+1)
        p1 = floor(crosslist(i));         % Crosspoints
        p2 = floor(crosslist(i+1));
        r1 = crosslist(i) - p1;           % Proportional residuals
        r2 = crosslist(i+1) - p2;
  
        if (p1 == p2)                     % If crosspoints on same step
          d = steplen(p1) .* (r2-r1);
        else                              % If crosspoints on different steps
          d = steplen(p1) .* (1-r1);        % Distance along first step
          d = d + steplen(p2) .* r2;        % Distance along last step
          if (p2-1 > p1)                    % Intermediate steps, if any
            d = d + sum(steplen((p1+1):(p2-1)));  % Intermediate steps
          end;
        end;
        intercross(i) = d;
      end;
    end;

    stats(7) = lenlist/2;
    if (~isempty(intercross))
      stats(8) = std(intercross);
    else
      stats(8) = 0;
    end;
  else
    stats(7:8) = [0 0];
  end;

  % Angular deviations between successive steps

  autocorr = NaN;
  angdev = NaN;

  if (n > 2)
    angdev = angledev(crds(1:(n-2),:),crds(2:(n-1),:),crds(3:n,:)); % Angular deviations
    angdev = abs(angdev);
    b = (angdev < pi/2);
    angdev(b) = pi/2 - angdev(b);
    angdev(~b) = -(angdev(~b) - pi/2);

    if (n > 5)                          % Autocorrelation of angular deviations
      a1 = angdev(1:(n-3));               
      a2 = angdev(2:(n-2));
      autocorr = corrcoef([a1,a2]);
      autocorr = autocorr(1,2);
    end;
  end;

  stats(9) = autocorr;

  if (isfinite(angdev))
    stats(10) = min(abs(angdev));
    stats(11) = max(abs(angdev));
    stats(12) = mean(abs(angdev));
    stats(13) = std(abs(angdev));
    stats(14) = mean(angdev);
    stats(15) = std(angdev);
  end;

  return;

