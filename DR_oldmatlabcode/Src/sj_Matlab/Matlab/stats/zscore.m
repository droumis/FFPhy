function [z,mu,sigma] = zscore(x)
%ZSCORE Standardized z score.
%   Z = ZSCORE(X) returns a centered, scaled version of X, 
%   the same size as X. For vector input X, output is the 
%   vector of z-scores Z = (X-mean(X))./std(X). For matrix 
%   input X, z-scores are computed using the mean and standard 
%   deviation along each column of X. For higher-dimensional 
%   arrays, z-scores are computed using the mean and standard 
%   deviation along the first non-singleton dimension.
%
%   The columns of Z have sample mean zero and sample standard
%   deviation one (unless a column of X is constant, in which
%   case that column of Z is constant at 0).
%
%   [Z,MU,SIGMA] = ZSCORE(X) also returns mean(X) in MU and 
%   std(X) in SIGMA.
%
%   See also MEAN, STD.

%   Copyright 1993-2006 The MathWorks, Inc. 
%   $Revision: 1.7.2.4 $  $Date: 2006/11/11 22:56:00 $

% [] is a special case for std and mean, just handle it out here.
if isequal(x,[]), z = []; return; end

% Compute X's mean and sd, and standardize it
mu = mean(x);
sigma = std(x);
sigma0 = sigma;
sigma0(sigma0==0) = 1;
z = bsxfun(@minus,x, mu);
z = bsxfun(@rdivide, z,sigma0);

