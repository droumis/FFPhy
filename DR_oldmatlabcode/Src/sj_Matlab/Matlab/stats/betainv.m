function x = betainv(p,a,b);
%BETAINV Inverse of the beta cumulative distribution function (cdf).
%   X = BETAINV(P,A,B) returns the inverse of the beta cdf with 
%   parameters A and B at the values in P.
%
%   The size of X is the common size of the input arguments. A scalar input  
%   functions as a constant matrix of the same size as the other inputs.    
%
%   BETAINV uses Newton's method to converge to the solution.
%
%   See also BETACDF, BETAFIT, BETALIKE, BETAPDF, BETARND, BETASTAT, ICDF.

%   Reference:
%      [1]     M. Abramowitz and I. A. Stegun, "Handbook of Mathematical
%      Functions", Government Printing Office, 1964.

%   Copyright 1993-2005 The MathWorks, Inc.
%   $Revision: 2.11.2.8 $  $Date: 2006/12/15 19:29:53 $

if nargin < 3
    error('stats:betainv:TooFewInputs',...
          'Requires at least three input arguments.');
end

[errorcode, p, a, b] = distchck(3,p,a,b);
if errorcode > 0
    error('stats:betainv:InputSizeMismatch',...
          'Non-scalar arguments must match in size.');
end

% Weed out any out of range parameters or edge/bad probabilities.
okAB = (0 < a & a < Inf) & (0 < b & b < Inf);
k = (okAB & (0 < p & p < 1));
allOK = all(k(:));

% Fill in NaNs for out of range cases, fill in edges cases when P is 0 or 1.
if ~allOK
    if isa(p,'single') || isa(a,'single') || isa(b,'single')
       x = NaN(size(k),'single');
    else
       x = NaN(size(k));
    end
    x(okAB & p == 0) = 0;
    x(okAB & p == 1) = 1;

    % Remove the bad/edge cases, leaving the easy cases.  If there's
    % nothing remaining, return.
    if any(k(:))
        if numel(p) > 1, p = p(k); end
        if numel(a) > 1, a = a(k); end
        if numel(b) > 1, b = b(k); end
    else
        return;
    end
end

% ==== Newton's Method to find a root of GAMCDF(X,A,B) = P ====

% Use the mean as a starting guess for q.
q = a ./ (a + b);
if isa(p,'single'), q = single(q); end

% Limit the number of iterations.
maxiter = 500;
reltol = eps(class(q)).^(3/4);

% Move starting values away from the boundaries.
q = max(reltol, min(1-reltol, q));

F = betacdf(q,a,b);
dF = F - p;
for iter = 1:maxiter
    % Compute the Newton step, but limit its size to prevent stepping out
    % of the unit interval.
    f = betapdf(q,a,b);
    h = dF ./ f;
    qNew = max(q/10, min(1-(1-q)/10, q - h)); % Avoid taking too large of a step
    
    % Break out of the iteration loop when the relative size of the last step
    % is small for all elements of q.
    done = (abs(h) <= reltol*q);
    if all(done(:))
        q = qNew;
        break
    end
    
    % Check the error, and backstep for those elements whose error did not
    % decrease.  The direction of h is always correct, because f > 0
    dFold = dF;
    F = betacdf(qNew,a,b);
    for j = 1:25 % If this fails, the outer loop may too
        dF = F - p;
        worse = (abs(dF) > abs(dFold)) & ~done;
        if ~any(worse(:))
            break
        end
        qNew(worse) = (q(worse) + qNew(worse))/2;
        F(worse) = betacdf(qNew(worse),a(worse),b(worse));
    end
    q = qNew;
end

badcdf = (abs(dF./F) > reltol.^(2/3));
if iter>maxiter || any(badcdf(:))   % too many iterations or cdf is too far off
    didnt = find(~done | badcdf, 1, 'first');
    if numel(a) == 1, abad = a; else abad = a(didnt); end
    if numel(b) == 1, bbad = b; else bbad = b(didnt); end
    if numel(p) == 1, pbad = p; else pbad = p(didnt); end
    warning('stats:betainv:NoConvergence',...
            'BETAINV did not converge for a = %g, b = %g, p = %g.',...
            abad,bbad,pbad);
end

% Broadcast the values to the correct place if need be.
if allOK
    x = q;
else
    x(k) = q;
end
