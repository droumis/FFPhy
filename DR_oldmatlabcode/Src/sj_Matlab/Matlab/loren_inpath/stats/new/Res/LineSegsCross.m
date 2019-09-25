% LineSegsCross: Given two line segments, P-Q and R-S, returns true (1) if the
%       line segments cross and false (0) if they don't.  Returns NaN if the
%       line segments are parallel or congruent.
%
%     Usage: bx = linesegscross(Q,R,S,T)
%
%         Q,R - endpoints of the first line segment; each is a vector of two
%                 coordinates.
%         S,T - endpoints of the second line segment; each is a vector of two
%                 coordinates.
%         -------------------------------------------------------------------
%         bx -  true if the segments cross, false if not.
%

function bx = linesegscross(Q,R,S,T)
  threshhold = 1e-6;

  q1 = Q(1);
  q2 = Q(2);
  r1 = R(1);
  r2 = R(2);
  s1 = S(1);
  s2 = S(2);
  t1 = T(1);
  t2 = T(2);

  denom = (q2-r2)*(t1-s1) - (r1-q1)*(s2-t2);
  if (abs(denom) < threshhold)
    bx = NaN;
    return;
  end;
  
  p1 = ((r1-q1)*(s1*t2-s2*t1)-(t1-s1)*(q1*r2-q2*r1))/denom;
  p2 = ((s2-t2)*(q1*r2-q2*r1)-(q2-r2)*(s1*t2-s2*t1))/denom;
  P = [p1 p2];
  
  PQ = eucl(P,Q);
  PR = eucl(P,R);
  QR = eucl(Q,R);
  
  PS = eucl(P,S);
  PT = eucl(P,T);
  ST = eucl(S,T);

  if (abs(PQ+PR-QR)<threshhold & abs(PS+PT-ST)<threshhold)
    bx = 1;
  else
    bx = 0;
  end;
  
  return;
  