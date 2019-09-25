function y = mnpdf(x,p)
%MNPDF Multinomial probability density function (pdf).
%   Y = MNPDF(X,PROB) returns the pdf for the multinomial distribution with
%   probabilities PROB, evaluated at each row of X. X and PROB are M-by-K
%   matrices or 1-by-K vectors, where K is the number of multinomial bins or
%   categories.  Each row of PROB must sum to one, and the sample sizes for
%   each observation (rows of X) are given by the row sums SUM(X,2).  Y is an
%   M-by-K matrix, and MNPDF computes each row of Y using the corresponding
%   rows of the inputs, or replicates them if needed.
%
%   See also MNRFIT, MNRVAL, MNRND

%   Copyright 2006 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2006/06/20 20:51:13 $

if nargin < 2
    error('stats:mnpdf:TooFewInputs', ...
          'Requires two input arguments.');
end

[m,k] = size(x);
if k < 1
    error('stats:mnpdf:NoCategories', ...
          'P must have at least one column.');
end
n = sum(x,2);

[mm,kk] = size(p);
if kk ~= k
    error('stats:mnpdf:InputSizeMismatch', ...
          'P and X must have the same number of columns.');
elseif mm == 1 % when m > 1
    p = repmat(p,m,1);
elseif m == 1
    m = mm;
    x = repmat(x,m,1);
    n = repmat(n,m,1);
elseif mm ~= m
    error('stats:mnpdf:InputSizeMismatch', ...
          'P and X must have the same number of rows, or either can be a row vector.');
end

outClass = superiorfloat(n,p);

xBad = any(x < 0 | x ~= round(x), 2);
pBad = any(p < 0 | 1 < p, 2) & abs(sum(p,2)-1) > eps(class(p));
nBad = n < 0 | round(n) ~= n;

xPos = (x > 0); % avoid 0*log(0), but let (pi==0) & (y>0) happen
xlogp = zeros(m,k,outClass);
xlogp(xPos) = x(xPos).*log(p(xPos)); % may be bad if p was, those won't be used
xlogp = sum(xlogp, 2);

% Initialize to return zero for noninteger or negative x
y = zeros(m,1,outClass); % y(xBad) = 0;

t = ~(xBad | pBad | nBad);
y(t) = exp(gammaln(n(t,:) + 1) - sum(gammaln(x(t,:) + 1), 2) + xlogp(t,:));

% Return NaN for invalid parameter values
y(pBad | nBad) = NaN;
