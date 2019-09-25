%   s = ANGLESTD(angle)
%       returns the standard deviation of the angles in the angle vector. 

function [s] = anglestd(a)

[m r] = anglemean(a);

s = sqrt(2 * (1 - r));
