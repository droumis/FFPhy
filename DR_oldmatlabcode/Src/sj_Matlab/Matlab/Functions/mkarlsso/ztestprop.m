% [p] = ztestprop(prop1, expected_proportion)
%       Returns the p value for the comparison of the proportion in prop1 to
%       the expected proportion.
%       prop1 should be a 1 x 2 vector with the number of positives in
%       the first column and the total number in the second column
function [p, z] = ztestprop(n1, expprop)

phat = n1(1) / n1(2);
qhat = 1 - phat;

stderrp = sqrt(phat * qhat / n1(2));

z = abs(expprop - phat) / stderrp;

p = (1 - normcdf(z)) * 2;
