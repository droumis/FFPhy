function [P,Z] = comparervalues(Rsquare,N)
%significance test for two Rsquare values 


R = sqrt(Rsquare);

z1 = .5*log((1+R(1))/(1-R(1)));
z2 = .5*log((1+R(2))/(1-R(2)));

Z = (z2-z1)/sqrt((1/(N(1)-3))+(1/(N(2)-3)));
P = (1 - normcdf(abs(Z))) * 2;

