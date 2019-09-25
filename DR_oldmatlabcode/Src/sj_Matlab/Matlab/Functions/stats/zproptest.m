function [P, z] = zproptest(x,n)
%[P, z] = zproptest(x,n)
% Computes z-test for two proportions
%x: [x1 x2]
%n: [n1 n2]


phat = (x(1) + x(2)) / (n(1) + n(2));
p1hat = x(1) / n(1);
p2hat = x(2) / n(2);
qhat = 1 - phat;

z = (p1hat-p2hat)/(sqrt(phat*qhat*( (1/(n(1))+(1/n(2)) ))));
P = (1 - normcdf(z)) * 2;