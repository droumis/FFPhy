%   m = ANGLEMEAN(angles)
%       returns the mean of the angles in the angle vector. 
%   [m r] = ANGLEMEAN(angles)
%       returns the mean of the angles in the angle vector and the length of
%       the mean resultant vector
%    

function [m, r] = anglemean(a)

% transform all of the angles into unit vectors
x = cos(a);
y = sin(a);

% sum all of the vectors and take the arctangent of the resultant vector

if (nargout == 2)
	r = dist([0 0], [sum(x) sum(y)]) / length(a);
end
m = atan2(sum(y), sum(x));

