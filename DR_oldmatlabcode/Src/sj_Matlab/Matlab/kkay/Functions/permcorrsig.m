function [r p] = permcorrsig(x, y, nperm)
% function [r p] = permcorrsig(x, y, nperm)
% 	returns the correlation and the p value associated with correlation 
% 	between vectors x and y based on a permutation test where nperm 
% 	correlations of x and a scrambled version of y are calculated.  
%	Note that this correlation is two-tailed, so it determines the p-value
%	for the absolute value of the correlation
c = corrcoef(x,y);
r = c(1,2);


n = length(y);
cdist = zeros(nperm,1);
% for each permulation we generate a random vector of indeces and recompute the
% correlation
for i = 1:nperm
    tmpind = randperm(n);
    ctmp = corrcoef(x,y(tmpind));
    cdist(i) = ctmp(1,2);
end
% sort the items so we can look up the index of the actual correlation
cdist = sort(abs(cdist));

% find the index of the element closest to the actual correlation.
cind = lookup(r, cdist);

% the p value is the distance from that index to the maximum index
p = 1 - cind / nperm;
