% AngDeviation: Given two line segments PQ and PR, with P being the common vertex,
%     finds the counterclockwise angle from PQ to PR.  If two or more points are
%     congruent, returns NaN.  Clockwise angles are returned as negative values.
%
%     Usage: theta = angdeviation(P,Q,R)
%

function theta = angdeviation(P,Q,R)
  threshold = 1e-6;
  
  p1 = P(1);
  p2 = P(2);
  q1 = Q(1);
  q2 = Q(2);
  r1 = R(1);
  r2 = R(2);

  Dpq = eucl(P,Q);                      % Step 1
  Dpr = eucl(P,R);
  
  if (Dpq<threshold | Dpr<threshold)    % Step 2
    theta = NaN;
    return;
  end;
  
  thetaQx = acos((q1-p1)/Dpq);          % Step 3
  thetaQy = asin((q2-p2)/Dpq);
  if (thetaQy >= 0)
    thetaQ = thetaQx;
  else
    thetaQ = 2*pi - thetaQ;
%     thetaQ = -(2*pi - thetaQx);
  end;

  thetaRx = acos((r1-p1)/Dpr);          % Step 4
  thetaRy = asin((r2-p2)/Dpr);
  if (thetaRy >= 0)
    thetaR = thetaRx;
  else
    thetaR = 2*pi - thetaRx;
%     thetaR = -(2*pi - thetaRx);
  end;
  
  theta = thetaQ - thetaR;              % Step 5
  
  if (theta<0)
    theta = theta + 2*pi;
  end;

  return;