 % script nlineg1 example
% generate some data
t=[0.:.01:1]';
p=[4 5];
y=exp(-p(1)*t).*sin(2*pi*p(2)*t);
% rand('normal');
y=y+.01*randn(size(y));  %add noise

beta0=[2	3];	  %initial guess
[ypred,beta] = nlin(t,y,'mdlnlin',beta0);
