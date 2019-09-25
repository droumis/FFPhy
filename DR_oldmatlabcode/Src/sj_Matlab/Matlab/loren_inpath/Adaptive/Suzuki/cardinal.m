% points = CARDINAL(cpx, cpy, npoints)
%          Retrns the x and y coordinates of npoints evenly spaced points on the spline
function [p] = cardinal(cpx, cpy, npoints)

t = 0;
s = (1 - t) / 2;
Mc = [(-s)  (2-s)    (s-2)   s;
     (2*s)  (s-3)    (3-2*s) (-s);
     (-s)     0        s      0;
     0        1        0       0];


p = zeros(npoints,2);
p(:,1) = (cpx(2):((cpx(end-1)-cpx(2))/(npoints-1)):cpx(end-1))';

% for each point, find the segment it is in
for i = 1:npoints
	% find the segment of the spline containing the current point %
	cseg = 2;
	while (cpx(cseg+1) < p(i,1))
		cseg = cseg + 1;
	end
	x1 = (p(i,1) - cpx(cseg)) / (cpx(cseg+1) - cpx(cseg));
	x2 = x1 * x1;
	x3 = x2 * x1;
	p(i,2) = (cpy(cseg-1)*(-.5*x3 + x2 -.5*x1) + ...
			   cpy(cseg)*(1.5*x3 - 2.5*x2 + 1) + ...
			   cpy(cseg+1)*(-1.5*x3 + 2.0*x2 + .5*x1) + ...
			   cpy(cseg+2)*(.5*x3 - .5*x2));  
end



%for i = 1:npoints
%	% find the segment of the spline containing the current point %
%	cseg = 2;
%	while (cpx(cseg+1) < p(i,1))
%		cseg = cseg + 1;
%	end
%	x1 = (p(i,1) - cpx(cseg)) / (cpx(cseg+1) - cpx(cseg));
%	x2 = x1 * x1;
%	x3 = x2 * x1;
%	p(i,2) = exp(cpy(cseg-1)*(-.5*x3 + x2 -.5*x1) + ...
%			   cpy(cseg)*(1.5*x3 - 2.5*x2 + 1) + ...
%			   cpy(cseg+1)*(-1.5*x3 + 2.0*x2 + .5*x1) + ...
%			   cpy(cseg+2)*(.5*x3 - .5*x2));  
%end
