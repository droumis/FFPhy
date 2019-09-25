function r = mvtrnd(C,df,cases)
%MVTRND Random matrices from the multivariate t distribution.
%   R = MVTRND(C,DF,N) returns an N-by-D matrix R of random numbers from
%   the multivariate t distribution with correlation parameters C and
%   degrees of freedom DF.
%
%   C is a D-by-D symmetric non-negative definite matrix, typically a
%   correlation matrix.  If its diagonal elements are not 1, MVTRND scales
%   C to correlation form.
%
%   Each row of R is generated as a multivariate normal random vector with
%   mean 0 and covariance C, divided by a scaled chi-square random variate
%   with DF degrees of freedom.
% 
%   Example:
%
%      C = [1 .4; .4 1]; df = 2;
%      r = mvtrnd(C, df, 500);
%      plot(r(:,1),r(:,2),'.');
%
%   See also MVNRND, MVTPDF, MVTCDF, TRND.

%   Copyright 1993-2006 The MathWorks, Inc. 
%   $Revision: 1.6.4.4 $  $Date: 2006/11/11 22:55:37 $

if (nargin < 2)
   error('stats:mvtrnd:TooFewInputs','MVTRND requires two arguments.');
end
if nargin<3 || isempty(cases)
    cases = 1;
end

[m,n] = size(C);
if (m ~= n)
   error('stats:mvtrnd:BadCorrelation',...
         'The correlation matrix C must be square');
end

df = df(:);
if ~(isscalar(df) || (isvector(df) && length(df) == cases))
   error('stats:mvtrnd:InputSizeMismatch',...
         'DF must be a scalar or a vector with CASES elements.');
elseif any(df <= 0)
   error('stats:mvtrnd:BadDF',...
         'DF must be positive.');
end

% Make sure C is a correlation matrix, then get Cholesky factor
s = diag(C);
if (any(s~=1))
   C = C ./ sqrt(s * s');
end
[T,err] = cholcov(C);
if (err ~= 0)
   error('stats:mvtrnd:BadCorrelation',...
         'C must be a symmetric positive semi-definite matrix.');
end

% Generate normal and sqrt(normalized chi-square), then divide
r = randn(cases, size(T,1)) * T;
x = sqrt(gamrnd(df./2, 2, cases, 1) ./ df);
r = r ./ x(:,ones(n,1));
