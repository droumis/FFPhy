function [theta,dev,dl,d2l] = levenb(y,x,theta)
%LEVENB	Outline of Levenberg modified Newton-Raphson optimisation algorithm.
%	Assume m files lik.m and derivs.m which calculate the likelihood
%	and its first and second derivates with respect to theta.
%       dev - minus twice the log-likihood
%       dl  - first derivative of log-likelihood
%       d2l - second derivative of log-likelihood

%	Gordon K Smyth, U of Queensland, Australia, gks@maths.uq.oz.au
%	Mar 23, 1995

% likelihood and derivatives at starting values
dev=lik(y,x,theta);
[dl,d2l]=derivs(y,x,theta);
epsilon=std(d2l(:))/1000;

% maximize likelihood using Levenberg modified Newton's method
iter=0;
while abs(dl'*(d2l\dl)/length(dl)) > tol,
   iter=iter+1;
   thetaold=theta;
   devold=dev;
   theta=thetaold-d2l\dl;
   dev=lik(y,x,theta);
   if (dev-devold)/(dl'*(theta-thetaold)) < 0,
      epsilon=epsilon/decr;
   else;
      while (dev-devold)/(dl'*(theta-thetaold)) > 0,
         epsilon=epsilon*incr;
         if epsilon>1e+15,
            disp('epsilon too large');
            return;
         end;
         theta=thetaold-(d2l-epsilon*eye(d2l))\dl;
         dev=lik(y,x,theta);
         disp('Epsilon'); disp(epsilon);
      end;
   end;
   [dl,d2l]=derivs(y,x,theta);
   disp('Iteration'); disp(iter);
   disp('Deviance'); disp(dev);
   disp('First derivative'); disp(dl');
   disp('Eigenvalues of second derivative'); disp(eig(d2l)');
end;
