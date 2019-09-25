function [beta,theta,dev,dl,d2l,p] = ordinal(y,x,print,beta,theta)
%ORDINAL Ordinal logistic regression.  ORDINAL(Y,X,1); regresses the ordinal
%	response  Y  on the design matrix  X  and prints summary results.
%
%	Suppose  Y  takes values in  k  ordered categories, and let  p_ij
%	be the cumulative probability that  Y(i)  falls in the j'th category
%	or higher.  The ordinal logistic regression model is
%
%	 logit(p_ij) = theta(j) + x_i'beta , i = 1,..,length(Y), j = 1,..,k-1,
%
%	where  x_i  is the i'th row of  X .  The number of ordinal
%	categories  k  is taken to be the number of distinct values of
%	round(Y) .  If  k = 2  the model is ordinary logistic regression.
%	X  is assumed to have full column rank.
%
%	ORDINAL(Y)  fits the model with baseline logit odds only.  The full
%	form is  [BETA,THETA,DEV,DL,D2L,P] = ORDINAL(Y,X,PRINT,BETA,THETA)
%	in which all output arguments and all input arguments except  Y  are
%	optional.  PRINT = 1  requests summary information about the fitted
%	model to be displayed.  PRINT = 2  requests information about
%	convergence at each iteration.  Other values request no information
%	to be displayed.  Input arguments BETA  and  THETA  give initial
%	estimates for beta and theta.  DEV holds minus twice the
%	log-likelihood,  DL  the log-likelihood derivative vector with
%	respect to theta and beta,  D2L  the second derivative matrix, and
%	P  the estimates of the cell probabilities, p_{ij}-p_{i,j+1}.

%	Gordon K Smyth, U of Queensland, Australia, gks@maths.uq.oz.au
%	Nov 19, 1990.  Last revision Aug 29, 1995.

%	Calls LIK, DERIVS and DMULT.

% check input
if nargin<2, x=[]; end;
y=round(y(:)); [my ny]=size(y); [mx nx]=size(x);
if (mx>0)&(mx~=my), disp('row dimension of x does not equal the length of y'); return; end;

% initial calculations
if mx>0, xstd=std(x); x=-x./(ones(mx,1)*xstd); end;
tol=1e-6; incr=10; decr=2;
ymin=min(y); ymax=max(y); yrange=ymax-ymin;
z =( y*ones(1,yrange) )==( ones(my,ny)*( ymin   :(ymax-1)) );
z1=( y*ones(1,yrange) )==( ones(my,ny)*((ymin+1): ymax   ) );
z=z(:,any(z)); z1=z1(:,any(z1)); [mz nz]=size(z);

% starting values
if nargin<3, print=0; end;
if nargin<4, g=cumsum(sum(z))'./my; theta=log(g./(1-g)); end;
if nargin<5, beta=zeros(nx,1); else beta=beta.*(xstd'); end;
tb=[theta; beta];

% likelihood and derivatives at starting values
[g,g1,p,dev]=lik(y,x,tb,z,z1);
[dl,d2l]=derivs(x,z,z1,g,g1,p);
epsilon=std(d2l(:))/1000;

% maximize likelihood using Levenberg modified Newton's method
iter=0;
while abs(dl'*(d2l\dl)/length(dl)) > tol,
   iter=iter+1;
   tbold=tb;
   devold=dev;
   tb=tbold-d2l\dl;
   [g,g1,p,dev]=lik(y,x,tb,z,z1);
   if (dev-devold)/(dl'*(tb-tbold)) < 0,
      epsilon=epsilon/decr;
   else;
      while (dev-devold)/(dl'*(tb-tbold)) > 0,
         epsilon=epsilon*incr;
         if epsilon>1e+15,
            disp('epsilon too large');
            return;
         end;
         tb=tbold-(d2l-epsilon*eye(size(d2l)))\dl;
         [g,g1,p,dev]=lik(y,x,tb,z,z1);
         disp('Epsilon'); disp(epsilon);
      end;
   end;
   [dl,d2l]=derivs(x,z,z1,g,g1,p);
   if print==2,
      disp('Iteration'); disp(iter);
      disp('Deviance'); disp(dev);
      disp('First derivative'); disp(dl');
      disp('Eigenvalues of second derivative'); disp(eig(d2l)');
   end;
end;

% tidy up output

theta=tb(1:nz,1);
if mx>0, beta=tb((nz+1):(nz+nx),1)./(xstd'); end;

if print>=1,
   disp('Number of iterations'); disp(iter);
   disp('  Deviance'); disp(dev);
   se=sqrt(diag(inv(-d2l)));
   disp('     Theta      SE');
   disp([theta se(1:nz,1)]);
   if mx>0,
      disp('     Beta       SE');
      disp([beta se((nz+1):(nz+nx),1)./(xstd')]);
   end;
end;

if nargout==6,
   if nx>0,
      e=( (x*tb((nz+1):(nz+nx),1))*ones(1,nz) )+( ones(my,ny)*theta' );
   else
      e=ones(my,ny)*theta';
   end;
   p=diff([zeros(my,ny) exp(e)./(1+exp(e)) ones(my,ny)]')';
end;
