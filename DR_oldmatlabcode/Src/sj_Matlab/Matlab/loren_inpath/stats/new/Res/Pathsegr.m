% PATHSEGR: Given a set of coordinates describing a path in p dimensions and a 
%           line-segment length, finds the residual left after stepping off the 
%           line segment along the path.  Steps to the end of the path but does not 
%           overstep the end.
%
%     Usage: [resid,nsteps,stepcrds] = pathsegr(path,steplen)
%
%           path =     [n x p] coordinates of the path, for n points in p dimensions.
%           steplen =  length of line segment.
%           --------------------------------------------------------------------------
%           resid =    residual distance left over at the end of the path.
%           nsteps =   number of steps of the line segment along the path (before 
%                        overrunning the path).
%           stepcrds = coordinates which the line segments intersect the path.
%

% RE Strauss, 4/11/98

function [resid,nsteps,stepcrds] = pathsegr(path,steplen)
  [n,p] = size(path);
  tol = 1e-5;

%workpath = path;
%close all;

  stepcrds = path(1,:);
  endpath = 0;
  step = 0;

  while (~endpath)
    step = step + 1;
    
    pathdist = path - ones(n,1)*stepcrds(step,:);     % Dists of path crds from current step
    pathdist = sqrt(sum((pathdist.^2)')') - steplen;  % Resid dists from current step

    if (all(pathdist<0))                              % If all resid dists are negative,
      endpath = 1;                                    %   we're at end of path

    else
      i = min(find(pathdist>0));                      % Closest path crd with positive resid

      A = stepcrds(step,:);
      B = path(i-1,:);
      C = path(i,:);
      c = sqrt(sum((B-C).*(B-C)));
      pl = 0;
      pu = c;

      if (A==B)                                         % If next intsect on next path segment,
        Pl = A;
        Pu = C;
        P = (Pl+Pu)./2;
        p = (pl+pu)./2;
        a = sqrt(sum((P-A).*(P-A))) - steplen;

        while (abs(a) > tol)                              % Binary search for intersection point
          if (a>0)
            Pu = P;
            pu = p;
          else
            Pl = P;
            pl = p;
          end;
          P = (Pl+Pu)./2;
          p = (pl+pu)./2;
          a = sqrt(sum((P-A).*(P-A))) - steplen;
        end;
  
      else                                              % If next intsect not on next path segment,
        cos_alpha = cos(pi - abs(angledev(A,B,C)));
        b = pathdist(i-1) + steplen;
  
        Pl = B;
        Pu = C;
        P = (Pl+Pu)./2;
        p = (pl+pu)./2;
      
        a = sqrt(b*b + p*p - 2*b*p*cos_alpha) - steplen;  % Use law of cosines to find dist to path
        while (abs(a) > tol)                              % Binary search for intersection point
          if (a>0)
            Pu = P;
            pu = p;
          else
            Pl = P;
            pl = p;
          end;
          P = (Pl+Pu)./2;
          p = (pl+pu)./2;
          a = sqrt(b*b + p*p - 2*b*p*cos_alpha) - steplen;
        end;
      end;
  
      stepcrds(step+1,:) = P;
      path(1:i-1,:) = ones(i-1,1)*stepcrds(step+1,:);
    end;
  end;

  nsteps = size(stepcrds,1)-1;                        % Final number of steps
  d = stepcrds(nsteps+1,:)-path(n,:);                 % Residual path length
  resid = sqrt(sum(d.*d));   

%figure;
%plot(workpath(:,1),workpath(:,2),'k');
%hold on;
%plot(stepcrds(:,1),stepcrds(:,2),'ko');
%hold off;
%axis('equal');

  return;



