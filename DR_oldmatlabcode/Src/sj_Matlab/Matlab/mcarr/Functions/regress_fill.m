% [L U slope intercept] = regress_fill(x,y,bin,nboot)
% Computes a bootstrap estimate of the regression coefficients and returns
% the lower and upper bounds at each bin.
% nboot specifies the number of resamples, Default = 1000
function [L U slope intercept] = regress_fill(x,y,bin,nboot)

if (nargin < 4)
    nboot = 1000;
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

%Run bootstrap resampling and regression
qboot = nan(nboot,2);
for s = 1:nboot
    boot = ceil(length(y)*rand(length(y),1));
    xboot = x(boot);
    yboot = y(boot);
    
    b = regress(yboot,[ones(length(xboot),1) xboot],0.05);
    qboot(s,1) = b(1);
    qboot(s,2) = b(2);
end

%Determine 95% confidence intervals at each bin location
L = nan(size(bin));
U = nan(size(bin));
for s = 1:length(bin)
    L(s) = prctile(qboot(:,2).*bin(s) + qboot(:,1),2.5);
    U(s) = prctile(qboot(:,2).*bin(s) + qboot(:,1),97.5);
end

slope = mean(qboot(:,2));
intercept = mean(qboot(:,1));
end