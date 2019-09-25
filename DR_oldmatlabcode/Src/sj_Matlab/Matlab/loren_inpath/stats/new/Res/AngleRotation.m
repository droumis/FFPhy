% AngleRotation: Given two line segments PQ and PR, with P being the common vertex,
%     finds the counterclockwise angle of rotation from PQ to PR.  
%     Returns NaN if P is congruent with Q or R.
%
%     Usage: theta = anglerotation(P,Q,R,{do_neg})
%
%         P,Q,R =   [n x 2] matrices specifying the cartesian coordinates of the
%                     points P, Q, and R.  If any matrix is size [1 x 2], 
%                     representing a single point, it is expanded to the same size 
%                     as the others.
%         do_neg =  optional boolean flag indicating, if true, that angles less than
%                     pi radians are to be returned as negative angles [default = 0].
%         ---------------------------------------------------------------------------
%         theta =   [n x 1] vector of counterclockwise angle of rotation, in radians.
%

% RE Strauss, 12/5/00
%   5/5/03 - misc. small changes; 
%            allow for negative angles;  
%            allow for P,Q,R to be matrices of N points rather than single points.

function theta = anglerotation(P,Q,R,do_neg)
  if (~nargin) help anglerotation; return; end;
  
  if (nargin < 4) do_neg = []; end;
  
  if (isempty(do_neg)) do_neg = 0; end;

  threshold = 1e-6;
  
  if (isvector(P))
    P = P(:)';
  end;
  if (isvector(Q))
    Q = Q(:)';
  end;
  if (isvector(R))
    R = R(:)';
  end;
  
  if (size(P,2)~=2 | size(P,2)~=2 | size(R,2)~=2)
    error('  AngleRotation: points must be 2-dimensional.');
  end;
  
  [ok,P,Q,R] = samelength(P,Q,R);
  if (~ok)
    error('  AngleRotation: incompatible point matrices.');
  end;
  
  [N,c] = size(P);
  theta = zeros(N,1);                   % Allocate output matrix
  
  for i = 1:N
    p1 = P(i,1); p2 = P(i,2);
    q1 = Q(i,1); q2 = Q(i,2);
    r1 = R(i,1); r2 = R(i,2);

    Dpq = eucl(P(i,:),Q(i,:));
    Dpr = eucl(P(i,:),R(i,:));
  
    if (Dpq<threshold | Dpr<threshold)
      theta(i) = NaN;
    else
      thetaQx = acos((q1-p1)/Dpq);
      thetaQy = asin((q2-p2)/Dpq);
      if (thetaQy >= 0)
        thetaQ = thetaQx;
      else
        thetaQ = 2*pi - thetaQx;
      end;

      thetaRx = acos((r1-p1)/Dpr);
      thetaRy = asin((r2-p2)/Dpr);
      if (thetaRy >= 0)
        thetaR = thetaRx;
      else
        thetaR = 2*pi - thetaRx;
      end;
  
      theta(i) = thetaR - thetaQ;
  
      while (theta(i)<0)
        theta(i) = theta(i) + 2*pi;
      end;
  
      if (do_neg & theta(i) > pi)
        theta(i) = theta(i) - 2*pi;
      end;
    end;
  end;

  return;