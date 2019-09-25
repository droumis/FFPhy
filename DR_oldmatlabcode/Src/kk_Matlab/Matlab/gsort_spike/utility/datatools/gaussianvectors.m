function data = gaussianvectors(numvectors, meanvector, datacov)
%GAUSSIANVECTORS   Returns vectors whose components are Gaussian RV's.
%   DATA = GAUSSIANVECTORS(NUMVECTORS, MEANVECTOR, DATACOVARIANCE) returns
%   an (M x N) matrix DATA, where M is NUMVECTORS and N is the number of
%   columns in MEANVECTOR.  The columns of DATA are each drawn from a
%   Gaussian distribution with mean given by the corresponding column in
%   MEANVECTOR.
%
%   The covariance of the various columns is specified by the
%   DATACOVARIANCE input, which must be a symmetric positive definite
%   matrix with dimension matching the number of columns in MEANVECTOR.

numdims = size(meanvector, 2);

% error checking
if ((size(datacov, 1) ~= size(datacov, 2)) || (size(datacov, 1) ~= numdims))
	error('Covariance matrix dimensions inconsistent with mean vector.');
end
if (all(all((datacov ~= datacov'))) || (det(datacov) <= 0))
	error('Covariance matrix must be symmetric positive definite.');
end

% We start with random vectors W, where W_k = N(0,1) and W_i and W_j
% are indpendent random variables for all i ~= j.  These vectors can
% be made to have the desired covariance by premultiplying them with
% the Cholesky factorization of the (symm pos def) covariance matrix.
% The desired mean is then added in to give Gaussian vectors with the
% desired statistics.
data = randn(numvectors, numdims) * chol(datacov) + repmat(meanvector, [numvectors, 1]);