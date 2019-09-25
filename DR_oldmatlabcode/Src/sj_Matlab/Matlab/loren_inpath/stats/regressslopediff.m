% [sig] = regressslopediff(x1, y1, x2, y2)
%         Returns the significance value of the comparison between the
%         slopes of the best fit linear regression to the two lines (x1,y1) and
%         (x2, y2). 
%         x1 and y1 must be the same size, and x2 and y2 must be the same size

function [sig] = regressslopediff(x1, y1, x2, y2)

% check the arguments
if (sum(size(x1) == size(y1)) ~= 2)
    error('x1 and y1 must be the same size');
end
if (sum(size(x2) == size(y2)) ~= 2)
    error('x2 and y2 must be the same size');
end


% to run regress we need to add a column of ones to x
% make sure the x's are n x 1
if (size(x1,1) == 1)
    x1 = x1';
    y1 = y1';
end
if (size(x2,1) == 1)
    x2 = x2';
    y2 = y2';
end

x1 = [ones(size(x1)) x1];
x2 = [ones(size(x2)) x2];


[b1, bint1, r, rint1, stats] = regress(y1, x1, .05);
[b2, bint2, r, rint1, stats] = regress(y2, x2, .05);

% The formula for the t value of the difference of the two slopes is
% t = (b1 - b2) / sqrt(sb1^2 + sb2^2) where
%	b1 is the slope of the first regression line
%	b2 is the slope of the second regression line
%	sb1 is the standard deviation of the estimate of b1
%	sb2 is the standard deviation of the estimate of b2
% sb1 = Syx / (Sx * sqrt(N-1)) where
% 	Syx = sqrt(sum(y - yest)^2) / (N-2)) (the standard error of the
% 		estimate)
%       Sx = sqrt(sum(x - xmean)^2 / (N - 1))
% From Statistical Methods in Psychology by David C. Howell, pg. 250

y1est = (b1' * x1')';
n1 = length(y1);
Syx1 = sqrt(sum((y1 - y1est).^2) / (n1 - 2));
Sx1 = sqrt(sum((x1(:,2) - mean(x1(:,2))).^2) / (n1 - 1));
sb1 = Syx1 / (Sx1 * sqrt(n1 - 1));

y2est = (b2' * x2')';
n2 = length(y2);
Syx2 = sqrt(sum((y2 - y2est).^2) / (n2 - 2));
Sx2 = sqrt(sum((x2(:,2) - mean(x2(:,2))).^2) / (n2 - 1));
sb2 = Syx2 / (Sx2 * sqrt(n2 - 1));

t = abs(b1(2) - b2(2)) / sqrt(sb1^2 + sb2^2) 


% lookup the p value for this t on N1 + N2 - 4 df. The value must be subtracted
% from one because we are looking on the positive side of the distribution, and
% multiple by 2 to get a two tailed comparison
sig = 2 * (1 - tcdf(t, n1 + n2 - 4));
