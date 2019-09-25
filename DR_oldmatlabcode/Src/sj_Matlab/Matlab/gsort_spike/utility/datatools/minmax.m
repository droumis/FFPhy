function [extrema,inds] = minmax(X)
%MINMAX            Simultaneous overall smallest and largest components.
%   Y = MINMAX(X) returns the minimum and maximum values, of the array X
%   such that Y(1) = MIN(X) and Y(2) = MAX(X).  For N-D arrays X,
%   MINMAX(X) is equivalent to MINMAX(X(:)).
%
%   [Y,I] = MINMAX(X) also returns the linear indices of the extrema such
%   that, Y(1) == X(I(1)) and Y(2) == X(I(2)).  When X has more than one
%   extremal element, the index of the first is returned.
%
%   X must be a real, double array.  NaN's are ignored.

%%%%% Argument checking.
if (~strcmp(class(X), 'double') || ~isreal(X)), 
    error('X must be a real, double array.');  
end;

[extrema(1),extrema(2),inds(1),inds(2)] = CORE_minmax(X(:));

%[extrema(1),inds(1)] = min(X(:));
%[extrema(2),inds(2)] = max(X(:));


