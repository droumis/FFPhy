% [pvalue] = regress_permutation(x,y,nboot)
% Computes a permutation test to determine whether the regression
% is significantly different from null hypothesis (slope == 0).
% returns the pvalue.
% nperm specifies the number of resamples, Default = 1000

function [pvalue] = regress_permutation(x,y,nperm)

if (nargin < 3)
    nperm = 1000;
end

%Check size of x and y
if size(x,2) ~= 1 && size(x,1) == 1
    x = x';
elseif size(x,2) ~= 1
	error('x must be a vector with size Nx1');
end
if size(y,2) ~=1 && size(y,1) == 1
    y = y';
elseif size(y,2) ~=1
    error('y must be a vector with size Nx1');
end

[b,bint,r,rint,stats] = regress(y,[ones(length(x),1) x],0.05);
Q = stats(1);

%Run permutation test nperm times
qperm = nan(nperm,1);

for s = 1:nperm
    %Reassign x terms
    i = randperm(length(x));
    xperm = x(i);
    
    %Compute regression
    [b,bint,r,rint,stats] = regress(y,[ones(length(xperm),1) xperm],0.05);
    
    %permutation statistic is R2
    qperm(s) = stats(1);
end

%Determine if regression is significant by comparing R2 to distribution of
%R2 obtained by permuting x values.

pvalue = mean(qperm > Q);

end