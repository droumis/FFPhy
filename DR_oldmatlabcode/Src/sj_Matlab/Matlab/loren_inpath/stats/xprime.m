function xdot =xprime(t,x,flag,k)
%function xdot =xprime(t,x,flag,k)
% k is extra parameter passed from ode45
% flag variable not used but must be in arg list

xdot(1)=-k(1).* x(1) - k(2).*x(2);
xdot(2)= k(1)*x(1) +k(2).*x(2);
% make sure output is a col vector
xdot=xdot(:);
