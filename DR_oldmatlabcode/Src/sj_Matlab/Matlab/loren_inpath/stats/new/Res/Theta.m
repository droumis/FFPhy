% THETA: Returns value (0-360) having the same order properties as the angle
%         made by p0 & pt with respect to the horizontal; does not return the
%         actual angle.  Useful for angular sorting of points.  Returns a flag 
%         if point coordinates of p0 & pt are identical.
%
%     Syntax: [value,ident] = theta(p0,pt)
%
%         p0 =    [1 x 2] vector specifying reference-point coordinates
%         pt =    [n x 2] matrix specifying point coordinates
%         ----------------------------------------------------------------
%         value = [n x 1] vector of corresponding angular values
%         ident = [n x 1] vector of identity flags (TRUE [=1] if the point is 
%                   identical to p0, FALSE [=0] if not)
%

% RE Strauss, 1/29/96

% Sedgewick, R. 1988  Algorithms, 2nd ed. Addison Wesley (p. 353).

function [value,ident] = theta(p0,pt)
  [n,p] = size(pt);
  value = zeros(n,1);
  ident = zeros(n,1);

  for i = 1:n
    if (p0(1,:)==pt(i,:))
      value(i) = 0;
      ident(i) = 1;
    else
      dx = pt(i,1)-p0(1,1);
      ax = abs(dx);
      dy = pt(i,2)-p0(1,2);
      ay = abs(dy);

      if ((ax<eps) & (ay<eps))
        t = 0;
      else
        t = dy / (ax+ay);
      end;

      if (dx<0)
        t = 2-t;
      elseif (dy<0)
        t = 4+t;
      end;

      value(i) = t*90;
    end;
  end;

  return;

