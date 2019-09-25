%	dist(x1, x2) computes the distance between points x1 and x2
%                if x1 and x2 are matrices, dist returns a column vector of
%                distances between corresponding points. If x1 is 1xN point
%                and x2 is MxN, dist returns the distance from the x1 to
%                each point in x2
function [d] = dist(x1, x2)

if ((size(x1,1) == 1) & (size(x2,1) > 1))
	tmpx1 = zeros(size(x2,1), size(x2,2));
	tmpx1(:,1) = x1(1);
	tmpx1(:,2) = x1(2);
	x1 = tmpx1;	
end

d = sqrt(sum(((x1 - x2).^2)')'); 
