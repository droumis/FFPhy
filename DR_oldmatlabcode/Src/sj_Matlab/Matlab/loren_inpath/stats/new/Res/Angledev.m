% AngleDev: Given two points P & Q defining a ray, and a third 'test' point R,
%           returns the signed angular deviation (in radians, -pi <= theta <= pi) 
%           of R from the ray (ie, the change in angle taken when moving along 
%           the path from P to Q to R).  Positive deviations are counterclockwise, 
%           negative deviations are clockwise.
%           Optionally returns the counterclockwise angular deviation 
%           (0 <= theta <= 2*pi).
%
%     Syntax: theta = angledev(P,Q,R,{no_neg})
%
%           P -       [n x 2] or [1 x 2] x,y coordinates of the P point.  
%                       If a single point is passed, it is used for all Q and 
%                       test points.
%           Q -       [n x 2] or [1 x 2] x,y coordinates of the Q point.  
%                       If a single point is passed, it is used for all 
%                       P and test points.
%           R -       [n x 2] matrix of point coordinates.
%           no_neg -  optional boolean no_neg indicating that the angles returned 
%                       are to be counterclockwise angular deviations [default=0].
%           ----------------------------------------------------------------------
%           theta - column vector of angular deviations, in radians.
%

% RE Strauss, 6/26/96
%   5/4/03 - rewrite, calling anglerotation().

function theta = angledev(P,Q,R,no_neg)
  if (~nargin) help angledev; return; end;
  
  if (nargin < 4) no_neg = []; end;
  
  if (isempty(no_neg)) no_neg = 0; end;

  if (isvector(P))
    P = P(:)';
  end;
  if (isvector(Q))
    Q = Q(:)';
  end;
  if (isvector(R))
    R = R(:)';
  end;
  
  [ok,P,Q,R] = samelength(P,Q,R);
  
  [N,vc] = size(P);
  [N,rc] = size(Q);
  [N,cc] = size(R);
  
  if (vc~=2 | rc~=2 | cc~=2)
    error('  ANGLEDEV: at least one input matrix is wrong size');
  end;
  
%   S = Q+(Q-P);
%   figure;
%   subplot(1,2,1);
%   plot([P(1),Q(1),R(1),S(1)],[P(2),Q(2),R(2),S(2)],'ko');
%   hold on;
%   plot([P(1),Q(1),R(1)],[P(2),Q(2),R(2)],'k');
%   plot([Q(1),S(1)],[Q(2),S(2)],'k:');
%   hold off;
%   axis equal;
  
  % Move P and Q to origin, carrying along line segments.
  
  for i = 1:N
    R(i,:) = R(i,:) - Q(i,:);
    Q(i,:) = Q(i,:) - P(i,:);
    P(i,:) = [0,0];
  end;

%   subplot(1,2,2);
%   plot([P(1),Q(1),R(1)],[P(2),Q(2),R(2)],'ko');
%   hold on;
%   plot([P(1),Q(1)],[P(2),Q(2)],'k');
%   plot([P(1),R(1)],[P(2),R(2)],'k');
%   hold off;
%   axis equal;

  % Find angles of rotation

  theta = anglerotation(P,Q,R,~no_neg);
  
  return;
