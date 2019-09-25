% [p] = ztestprop2(prop1, prop2)
%       Returns the p value for the comparison of the proportions in prop1 and
%       prop2.
%       prop1 and prop2 should be 1 x 2 vectors with the number of positives in
%       the first column and the total number in the second column
% Originally, ztestprop2, ztestn2
function [p, stderrp, stderrp1, stderrp2] = sj_ztestprop2(n1, n2)

phat = (n1(1) + n2(1)) / (n1(2) + n2(2));
p1hat = n1(1) / n1(2);
p2hat = n2(1) / n2(2);
qhat = 1 - phat;

stderrp = sqrt(phat * qhat * (1 / n1(2) + 1 / n2(2)));

z = abs(p1hat - p2hat) / stderrp;

p = (1 - normcdf(z)) * 2;  % this is two-tailed. For one-tailed, dont multiply by 2

stderrp1 = sqrt(p1hat*(1-p1hat)*(1/n1(2)));
stderrp2 = sqrt(p2hat*(1-p2hat)*(1/n2(2)));
