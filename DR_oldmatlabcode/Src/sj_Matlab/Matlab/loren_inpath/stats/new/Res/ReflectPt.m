% ReflectPt: Reflects a point about a line, given two points defining a line segment.
%
%     Usage: rpt = reflectpt(pt,line)
%
%         pt =   2-element vector of point coordinates.
%         line = 4-element vector of two sets of point coordinates: [p1 p2, q1 q2].
%         -------------------------------------------------------------------------
%         rpt =  2-element row vector of reflected point coordinates.
%

% Bookstein et al. 1985. Morphometrics in Evolutionary Biology, appendix A.1.11.
%   Note: couldn't get A.1.11 to work.

% RE Strauss, 4/22/03
%   5/14/03 - added error message.

function rpt = reflectpt(pt,line)
  pt = pt(:);
  line = line(:);
  
  if (length(pt)~=2 | length(line)~=4)
    error('  ReflectPt: point or line vector is wrong length.');
  end;

  pts = [line(1:2)'; line(3:4)'; pt'];
  [newpts,theta] = regrot(1,2,pts);
  newpts(3,2) = -newpts(3,2);
  pts = rotate(newpts,-theta,newpts(1,:));
  pts = pts + ones(3,1)*line(1:2)';
  rpt = pts(3,:);

  return;
  