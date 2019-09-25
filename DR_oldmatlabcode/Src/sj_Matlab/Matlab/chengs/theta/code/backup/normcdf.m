function p = normcdf(x,mu,sigma)
%NORMCDF Normal cumulative distribution function (cdf).
%   P = NORMCDF(X,MU,SIGMA) computes the normal cdf with mean MU and
%   standard deviation SIGMA at the values in X.
%
%   The size of P is the common size of X, MU and SIGMA. A scalar input  
%   functions as a constant matrix of the same size as the other inputs.    
%
%   Default values for MU and SIGMA are 0 and 1 respectively.
%
%   See also NORMINV, ERF, ERFC, ERFINV, ERFCINV.

%   Copyright 1993-2002 The MathWorks, Inc. 
%   $Revision: 2.12 $  $Date: 2002/01/17 21:31:30 $

if nargin < 2, mu = 0; end
if nargin < 3, sigma = 1; end

[errorcode x mu sigma] = distchck(3,x,mu,sigma);
if errorcode > 0
    error('Requires non-scalar arguments to match in size.');
end

% It is numerically preferable to use the complementary error function
% and normcdf(x) = 0.5*erfc(-x/sqrt(2)) to produce accurate results
% approaching zero for large negative x.

p = 0.5 * erfc(-(x-mu)./(sqrt(2)*sigma));
