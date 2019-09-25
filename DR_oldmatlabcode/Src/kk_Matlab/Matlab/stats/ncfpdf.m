function Y = ncfpdf(x,nu1,nu2,delta)
%NCFPDF Noncentral F probability density function (pdf).
%   Y = NCFPDF(X,NU1,NU2,DELTA) returns the noncentral F pdf with numerator 
%   degrees of freedom (df) NU1, denominator df NU2, and noncentrality
%   parameter DELTA, at the values in X.
%
%   The size of Y is the common size of the input arguments. A scalar input  
%   functions as a constant matrix of the same size as the other inputs.     
%
%   See also NCFCDF, NCFINV, NCFRND, NCFSTAT, FPDF, PDF.

%   Reference:
%      Johnson, Kotz, and Balakrishnan, "Continuous Univariate
%        Distributions, Vol. 2" (2nd edition), Wiley, 1995, eq. 30.7.

%   Copyright 1993-2004 The MathWorks, Inc. 
%   $Revision: 2.15.4.7 $  $Date: 2004/12/06 16:37:49 $

if nargin < 4, 
    error('stats:ncfpdf:TooFewInputs','Requires four input arguments.');
end

[errorcode, x, nu1, nu2, delta] = distchck(4,x,nu1,nu2,delta);

if errorcode > 0
    error('stats:ncfpdf:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

[m,n] = size(x);

% Initialize Y to zero.
if isa(x,'single') || isa(nu1,'single') || isa(nu2,'single') || isa(delta,'single')
   y = zeros(m,n,'single');
   yeps = eps('single');
else
   y = zeros(m,n);
   yeps = eps;
end
Y = y;

% Take care of special cases:  invalid parameters, central distribution
% (delta=0), and edge case (x=0)
k = (nu1 <= 0 | nu2 <= 0 | x < 0 | delta < 0);
Y(k) = NaN;
Y(x==0 & nu1<=2 & ~k) = Inf;
k2 = (x==0 & nu1==2 & ~k);
if any(k2)
    Y(k2) = exp(-delta(k2)/2) ./ beta(nu1(k2)/2,nu2(k2)/2);
end
central = (delta==0) & ~k;
if any(central)
    Y(central) = fpdf(x(central),nu1(central),nu2(central));
end
k1 = ~(k | k2 | (x==0) | central);
if ~any(k1(:)), return; end

% Use good indices only, and make everything a vector
y = y(k1);
nu1 = nu1(k1);
nu2 = nu2(k1);
x = x(k1);
delta = delta(k1);

% to simply computation, pre-divide nu1,nu2 and delta
nu1 = nu1/2; nu2 = nu2/2; delta = delta/2;

% Sum an infinite series each way from a decent starting point
j0 = floor(delta);
g = nu1.*x./nu2;
olddelta = 0;
b = beta(nu1+j0,nu2);
todo = 1:numel(delta);
N1 = nu1;
N2 = nu2;
D = delta;
P = exp(-D + j0.*log(D) - gammaln(j0+1));
G = exp((j0+N1-1).*log(g) - (j0+N1+N2).*log(1+g));
P0 = P;
G0 = G;
g0 = g;
b0 = b;
j = j0;
niter = 0;
while(true)     % sum from j0 upward until contributions are small
   % In principle we are going to sum tmp/b over j, where:
   %   b      = beta(j + nu1,nu2);
   %   tmp    = poisspdf(j,delta).*(g.^(j-1+nu1))./((1+g).^(j + nu1 + nu2));
   % However, we computed the beta, poisson, and g parts for j=j0 above
   % the loop, and we will use recurrence relations to update them each time.

   tmp    = P .* G;
   deltay = tmp./b;
   y(todo) = y(todo) + deltay;
   
   % Find elements that have not yet converged
   mask = ((deltay >= yeps * (y(todo)+yeps)) | ...
           (deltay >= olddelta)) & (deltay>0);
   if ~any(mask(:)), break; end
   niter = niter + 1;
   if (niter == 1000) 
      warning('stats:ncfpdf:NoConvergence',...
              'Failed to converge, use the result cautiously.');
      break;    
   end

   % Update for next iteration, possibly working on fewer array elements
   if ~all(mask(:))
      todo = todo(mask);
      b = b(mask);
      g = g(mask);
      deltay = deltay(mask);
      N1 = N1(mask);
      N2 = N2(mask);
      D = D(mask);
      j = j(mask);
   end
   olddelta = deltay;
   b = b .* (j+N1) ./ (j+N1+N2);
   P = P(mask) .* D ./ (j+1);
   G = G(mask) .* g ./ (1+g);
   j = j + 1;
end

% Now count down
todo = 1:numel(delta);
N1 = nu1;
N2 = nu2;
D = delta;
P = P0;
G = G0;
g = g0;
b = b0;
j = j0-1;
mask = (j>=0);
olddelta = zeros(size(mask));
while(true)
   if ~any(mask(:)), break; end
   if ~all(mask(:))
      todo = todo(mask);
      b = b(mask);
      g = g(mask);
      N1 = N1(mask);
      N2 = N2(mask);
      D = D(mask);
      j = j(mask);
      olddelta = olddelta(mask);
   end
   b = b .* (j+N1+N2) ./ (j+N1);
   P = P(mask) .* (j+1) ./ D;
   G = G(mask) .* (1+g) ./ g;
   tmp    = P .* G;
   deltay = tmp./b;
   y(todo) = y(todo) + deltay;
   
   % Find elements that have not yet converged
   j = j - 1;
   mask = (j>=0) & ((deltay >= yeps * (y(todo)+yeps)) | ...
                    (deltay >= olddelta)) & (deltay>0);
   
   olddelta = deltay;
end

Y(k1) = nu1.*y./nu2;

