function y=modode5(k,tspan,x0)
%function y=modode5(k,tspan,x0)
%% example  to fit ode parameters
% x0(init condition) is extra passed parameter

[t,x12]=ode45('xprime',tspan,x0,[],k);
% xprime must be in separate file

%% take second x12 for output y, this is the one being fitted
y=x12(:,2);




