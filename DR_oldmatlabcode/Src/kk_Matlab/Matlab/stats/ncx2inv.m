function x = ncx2inv(p,v,delta)
%NCX2INV Inverse of the noncentral chi-square cdf.
%   X = NCX2INV(P,V,DELTA)  returns the inverse of the noncentral chi-square
%   cdf with parameters V and DELTA, at the probabilities in P.
%
%   The size of X is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   NCX2INV uses Newton's method to converge to the solution.
%
%   See also NCX2CDF, NCX2PDF, NCX2RND, NCX2STAT, CHI2INV, ICDF.

%   Copyright 1993-2005 The MathWorks, Inc.
%   $Revision: 2.9.4.8 $  $Date: 2006/02/21 20:44:09 $

if nargin<3
    error('stats:ncx2inv:TooFewInputs','Requires three input arguments.');
end
[errorcode p v delta] = distchck(3,p,v,delta);

if errorcode > 0
    error('stats:ncx2inv:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

% Initialize x to return NaN for arguments are outside their respective limits.
if isa(p,'single') || isa(v,'single') || isa(delta,'single')
   x = NaN(size(p),'single');
   crit = sqrt(eps('single'));
else
   x = NaN(size(p));
   crit = sqrt(eps);
end

% If some delta==0, call chi2inv for those entries, call ncx2inv for other entries.
f = (delta == 0);
if any(f(:))
   x(f) = chi2inv(p(f),v(f));
   f = ~f;
   if any(f(:)), x(f) = ncx2inv(p(f),v(f),delta(f)); end
   return
end

okParams = (v > 0 & delta >= 0);

% The inverse cdf of 0 is 0, and the inverse cdf of 1 is Inf.
x(p == 0 & okParams) = 0;
x(p == 1 & okParams) = Inf;

k = find(p > 0 & p < 1 & okParams);
pk = p(k);

% Newton's Method
% Permit no more than count_limit interations.
count_limit = 100;
count = 0;

% Supply a starting guess for the iteration.
%   Use a method of moments fit to the lognormal distribution.
mn = v(k) + delta(k);
variance = 2*(v(k) + 2*delta(k));
temp = log(variance + mn .^ 2);
mu = 2 * log(mn) - 0.5 * temp;
sigma = -2 * log(mn) + temp;
xk = exp(norminv(pk,mu,sigma));

h = ones(size(xk),class(xk));

% Break out of the iteration loop for three reasons:
%  1) the last update is very small (compared to x)
%  2) the last update is very small (compared to sqrt(eps))
%  3) There are more than 100 iterations.
F = ncx2cdf(xk,v(k),delta(k));
while(count < count_limit),

    count = count + 1;
    f = ncx2pdf(xk,v(k),delta(k));
    h = (F - pk) ./ f;
    
    % Avoid stepping too far
    xnew = max(xk/5, min(5*xk, xk - h));

    % Back off if the step gives a worse result
    newF = ncx2cdf(xnew,v(k),delta(k));
    while(true)
       worse = (abs(newF-pk) > abs(F-pk)*(1+crit)) & ...
               (abs(xk-xnew)> crit*xk);
       if ~any(worse), break; end
       xnew(worse) = 0.5 * (xnew(worse) + xk(worse));
       newF(worse) = ncx2cdf(xnew(worse),v(k(worse)),delta(k(worse)));
    end
    h = xk-xnew;
    x(k) = xnew;
    mask = (abs(h)>crit*abs(xk)) & (abs(h)>crit);
    if ~any(mask), break; end
    k = k(mask);
    xk = xnew(mask);
    F = newF(mask);
    pk = pk(mask);
end

% Store the converged value in the correct place
x(k) = xk;

if count == count_limit
    fprintf('\nWarning: NCX2INV did not converge.\n');
    str = 'The last step was:  ';
    outstr = sprintf([str,'%13.8f'],h);
    fprintf(outstr);
end
