function ip= findPeaks(y,steps)
% function ip= findPeaks(y,steps)
%
% Given a vector y, return the indices, at which local maxima occur.

sz= size(y);

if(sz(1)>1) y= y'; end

n= length(y);
yl= [inf, y(1:n-1)];
yr= [y(2:n), inf];

ip= find(yl < y & yr < y);
