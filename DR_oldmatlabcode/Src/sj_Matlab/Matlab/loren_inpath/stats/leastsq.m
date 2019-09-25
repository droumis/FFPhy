% [b1 b2 sigb1 r sigr] = leastsq(y,x)
%           Returns the result of a least squares fit of a line through
%           the points (x,y): y = b1*x + b2
%           b1 is the slope
%           b2 is the intercept
%           sigb1 is the p value for b1 (as compared to zero)
%           r is the correlation coefficient
%           sigr is the p value for r (assumes two tailed test)
function [b1, b2, sigb1, r, sigr] = leastsq(y,x)

n = length(x);
if (n ~= length(y))
	error('Error in leastsq: x and y must be vectors of the same length');
elseif (n == 1)
	b1 = 0;
	b2 = 0;
	sigb1 = 1;
	r = 0;
	sigr = 1;
	return
end

xmean = sum(x) / n;
ymean = sum(y) / n;

b1 = (sum(x.*y) - n * xmean * ymean) / (sum(x.^2) - n*xmean^2);
b2 = ymean - b1 * xmean;

% get the significance of b1
yhat = b1 * x + b2;
Syx = sqrt(sum((y - yhat).^2) / (n - 2));
tval = abs(b1 * std(x) * sqrt(n - 1) / Syx);
sigb1 = 1 - tcdf(tval, n-2);

tmpr = corrcoef(x,y);
r = tmpr(1,2);

sigr = corrsig(r, n) / 2;


