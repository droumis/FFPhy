function [a,r,X,Y]= circstat(phi)
%function [a,r,X,Y]= circstat(phi)
%
% circular descriptive statistics
%   If phi is matrix, work down the columns
%   a     : circular mean 
%   r     : dispersion
%   X     : mean cos(phi)
%   Y     : mean sin(phi)

x= cos(phi); y= sin(phi);
X= nanmean(x);  Y= nanmean(y);
r= sqrt(X.^2+Y.^2);
a= atan2(Y,X);
