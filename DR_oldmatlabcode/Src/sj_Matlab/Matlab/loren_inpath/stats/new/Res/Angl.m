% ANGL: Given three points Q,P,R in two dimensions, finds the smaller absolute
%       angle (in radians) at P between the line segments PQ and PR (range 0-pi).  
%       Or, optionally returns the counterclockwise angle from PQ to PR 
%       (range 0-2*pi).
%
%     Usage: theta = angl(Q,P,R,counterclock)
%
%       Q,P,R = matching [n x 2] matrices of point coordinates.  If P is a 
%                 single point, it is expanded to a size compatible with Q and R.
%       counterclock = optional boolean flag indicating, if true, that the 
%                 counterclockwise angle is to be returned.  If false, the 
%                 smaller absolute angle is returned [default].
%       -------------------------------------------------------------------------
%       theta = [n x 1] vector of angles.
%

% RE Strauss, 10/1/94
%   12/2/98 -  return theta=NaN for congruent points.
%   3/19/00 -  added option for counterclockwise angles.
%   12/10/01 - fixed bug in returning counterclockwise angles with congruent points.

% Bookstein et al. 1985, Appendix A.1.6.

function theta = angl(Q,P,R,counterclock)
  if (nargin < 2) P = []; end;
  if (nargin < 3) R = []; end;
  if (nargin < 4) counterclock = []; end;

  if (isempty(counterclock))
    counterclock = 0;
  end;

  [nQ,pQ] = size(Q);
  [nP,pP] = size(P);
  [nR,pR] = size(R);

  if (pQ~=2 & nQ==2)          % Transpose input matrices if necessary
    Q = Q';
    [nQ,pQ] = size(Q);
  end;
  if (pP~=2 & nP==2)          
    P = P';
    [nP,pP] = size(P);
  end;
  if (pR~=2 & nR==2)          
    R = R';
    [nR,pR] = size(R);
  end;

  if (nP==1)                  % If P is a single point, expand it
    P = ones(nQ,1)*P;
    nP = nQ;
  end;

  if (nQ~=nR | nP~=nR)
    error('  ANGL: input matrices not compatible');
  end;

  Q = Q-P;                                % Translate P to origin
  R = R-P;

  Dpq = sqrt(Q(:,1).^2 + Q(:,2).^2);      % Euclidean distances
  Dpr = sqrt(R(:,1).^2 + R(:,2).^2);

  theta = NaN*ones(nP,1);
  i = find((Dpq>eps) & (Dpr>eps));

  if (~isempty(i))
    if (counterclock)
      [rQ,thetaQ] = polarcrd(Q(i,1),Q(i,2));
      [rR,thetaR] = polarcrd(R(i,1),R(i,2));
      theta(i) = thetaR - thetaQ;
    else
      theta(i) = acos((Q(i,1).*R(i,1) + Q(i,2).*R(i,2))./(Dpq.*Dpr));
    end;

    i = find(theta < 0);
    while (~isempty(i))
      theta(i) = theta(i) + 2*pi;
      i = find(theta < 0);
    end;

    b = find(theta < 1e-6);
    if (sum(b)>0)
      theta(b) = zeros(length(b),1);
    end;
  else
    theta = NaN*ones(nP,1);
  end;

  return;
