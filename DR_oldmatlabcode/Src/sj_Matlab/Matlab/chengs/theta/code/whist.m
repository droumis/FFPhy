function [h, nh]= whist(x, weights, edges)
% function whist(x, weights, edges)
%
% weighted histograms
% x and weights should be vectors

if size(x,1) > 1 & size(x,2) > 1
    error('can only take vectors');
else
    if size(x,2) > 1
        x= x';
    end
end
if size(weights,1) > 1 & size(weights,2) > 1
    error('can only take vectors');
else
    if size(weights,2) > 1
        weights= weights';
    end
end
if size(x,1)~= size(weights,1)
    error('x and weights have to be same length');
end

X= sort([x weights]);

nbins= length(edges);
nx= size(X,1);
h= zeros(nbins,1);
nh= zeros(nbins,1);
ix= 1;
for ih= 1:nbins
    while ix <= nx & ~(X(ix,1)>edges(ih));
        nh(ih)= nh(ih)+1;
        h(ih)= h(ih)+X(ix,2);
        ix= ix+1;
    end
end
nh(find(~nh))=nan;
h= h./nh;
h= h(2:nbins);
nh= nh(2:nbins);
