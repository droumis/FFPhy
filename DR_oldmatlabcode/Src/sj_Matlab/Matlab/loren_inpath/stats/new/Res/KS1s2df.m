% KS1S2Df: Finds the proportions of points in the four quadrants centered on a given point.
%          Called by KS1S2D().
%

% RE Strauss, 5/28/03

function quadprop = ks1s2df(pts)
  npts = size(pts,1);
  q1 = zeros(npts,1);
  q2 = zeros(npts,1);
  q3 = zeros(npts,1);
  q4 = zeros(npts,1);
  
  for ip = 1:npts                     % For each point,
    cp = pts(ip,:);                     % Current point
    cpx = cp(1);
    cpy = cp(2);
    p = pts - ones(npts,1)*cp;          % Center all on current point
    q1 = sum(p(:,1)>cpx & p(:,2)>cpy);
    q2 = sum(p(:,1)<cpx & p(:,2)>cpy);
    q3 = sum(p(:,1)<cpx & p(:,2)<cpy);
    q4 = sum(p(:,1)>cpx & p(:,2)<cpy);
  end;
  quadprop = [q1 q2 q3 q4]/(npts-1);

  return;
  