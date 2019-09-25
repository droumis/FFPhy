%   newloc = PROJECTPOINT(point, linepoints)
%           Projects the point on to the line segments defined by linepoints.
%           Finds the closest end point of a line segment, projects the point
%           to the one or two segments connected to that endpoint, and chooses
%           the location of the point of the intersection of the perpendicular
%           to the line segment drawn througn the point. If no
%           such point can be found, the end point of the closest segment
%           is used. 
%           Returns newloc, a 1x5 vector with x, y, distance, onseg,
%           and segnum elements
function [newloc] = JY_projectpoint(p, lpnts)

% this function really should be rewritten so that it first computes the pixel
% coordinates of all of the points on the lines and then finds the minimum
% distance from point the the set of pixels along linepoints. 

% find the endpoint closest to p
ptmp = zeros(size(lpnts));
ptmp(:,1) = p(1);
ptmp(:,2) = p(2);

nseg = size(lpnts,1) - 1;
startp = lpnts(1:end-1,:);
endp = lpnts(2:end,:);

d = zeros(size(startp,1),1);
m = zeros(size(startp,1),1);
a = zeros(size(startp,1),1);
c = zeros(size(startp,1),1);
onseg = zeros(size(startp,1),1);
perpm = zeros(size(startp,1),1);
e = zeros(size(startp,1),1);
intp = zeros(size(startp,1),3);


% determine which lines are vertical
vertical = (endp(:,1) == startp(:,1));
horizontal = (endp(:,2) == startp(:,2));
vert = find(vertical);
nvert = find(~vertical);
horiz = find(horizontal);
nvnh = find(~vertical & ~horizontal);

% if the line is vertical, the distance is the difference in x coordinates
d(vert) = abs(startp(vert,1) - p(1));

% if the line is horizontal, the distance is the difference in y coordinates
d(horiz) = abs(startp(horiz,2) - p(2));

% for the rest, get the equation for the lines ax + by = c, setting b to 1
if (~isempty(nvnh))
    m(nvnh) = (endp(nvnh,2) - startp(nvnh,2))./(endp(nvnh,1) - startp(nvnh,1));
    a(nvnh) = -m(nvnh);
    c(nvnh) = startp(nvnh,2) - m(nvnh) .* startp(nvnh,1);
    % compute the distance
    d(nvnh) = abs(a(nvnh) * p(1) + p(2) - c(nvnh)) ./ sqrt(a(nvnh).^2 + 1);
end




%d= abs(det([endp-startp;p-startp]))/norm(endp-startp);

intp(:,3) = d;

% determine, for each point, whether the perpendicular projection is within
% the segment

% there are three cases
% Case 1: vertical line
%    check if the y coordinates are between startp and endp
if (~isempty(vert))
    onseg(vert) = ((p(2) > min([startp(vert,2) endp(vert,2)]')') & ... 
                   (p(2) < max([startp(vert,2) endp(vert,2)]')'));
    intp(vert,1) = startp(vert,1);
    intp(vert,2) = p(2);
end

% Case 2: horizontal line
%    check if the x coordinates are between startp and endp
if (~isempty(horiz))
    onseg(horiz) = ((p(1) > min([startp(horiz,1) endp(horiz,1)]')') & ...
                   (p(1) < max([startp(horiz,1) endp(horiz,1)]')'));
    intp(horiz,1) = p(1);
    intp(horiz,2) = startp(horiz,2);
end

% Case 3: neither horizontal nor vertical
%    compute the equation of the perpendicular line (y = -1/m*x + e) and solve
%    for the intersection point
if (~isempty(nvnh))
    perpm(nvnh) =  1 ./ a(nvnh);
    e(nvnh) = p(2) - perpm(nvnh) * p(1); 
    intp(nvnh,1) = (e(nvnh) - c(nvnh)) ./ (m(nvnh) - perpm(nvnh)); 
    intp(nvnh,2) = m(nvnh) .* intp(nvnh,1) + c(nvnh);
    onseg(nvnh) = ((intp(nvnh,1) > min([startp(nvnh,1) endp(nvnh,1)]')') & ...
                   (intp(nvnh,1) < max([startp(nvnh,1) endp(nvnh,1)]')') & ...
                   (intp(nvnh,2) > min([startp(nvnh,2) endp(nvnh,2)]')') & ...
                   (intp(nvnh,2) < max([startp(nvnh,2) endp(nvnh,2)]')'));
end

intp(:,4) = onseg;
intp(:,5) = (1:size(intp(:,1),1))';

% for all of the points that are not on a segment, change the distance to be
% the distance to the closest endpoint
tmp = find(~onseg);

% calculate the distance from the point to each of the endpoints of segments
tmppnt = zeros(size(lpnts,1), 2);
tmppnt(:,1) = p(1);
tmppnt(:,2) = p(2);
d = dist(tmppnt, lpnts);

% set the distance to each segment where there is no on segment perpendicular
% projection to be the distance to the closest endpoint and the point to be
% that endpoint
intp(tmp,3) = d(tmp);
intp(tmp,1:2) = lpnts(tmp,:);

% if the last value in tmp is the last segment, pick the minimum of the
% distance to the last lpnt and the distance to the next to last lpnt
if (~isempty(tmp) & (tmp(end) == nseg) & (d(tmp(end)) > d(tmp(end)+1)))
    % set the last point of intp to be the end of the last segment
    intp(tmp(end),3) = d(tmp(end)+1);
    intp(tmp(end),1:2) = lpnts(tmp(end)+1,:);
end


%[newdist minind] = min(d');
%intp(tmp,3) = newdist';
% the new value for the point is the closest endpoint whose index is tmp if
% minind is 1 and tmp+1 if minind is 2;
%intp(tmp,1:2) = lpnts((tmp+minind'-1),:);

% sort the intersection points by the distance element
intp = sortrows(intp,3);

% return all of the values
newloc = intp;
