function L = poisson_log_likelihood(counts,lambda,timestep)
%POISSON_LOG_LIKELIHOOD Compute log likelihood of intensity function for a discretized Poisson process
%
%   L = POISSON_LOG_LIKELIHOOD(COUNTS,LAMBDA,TIMESTEP) computes the
%   log-likelihood of the intensity function LAMBDA of an inhomogeneous Poisson
%   process, given observed COUNTS in bins of duration TIMESTEP. COUNTS and
%   LAMBDA must be vectors of the same length. COUNTS must contain non-negative
%   real integers, and LAMBDA must contain non-negative real finite
%   floating-point values. TIMESTEP must be a real positive finite
%   floating-point scalar. TIMESTEP and LAMBDA must be expressed in
%   commensurable units.
%
%Written by SMK, 2010 January 29.
%

if ~isvector(counts) || ~isreal(counts) || ~all(counts >= 0) || ...
    ~all(isfinite(counts)) || ~all(round(counts) == counts)
  error('COUNTS must be a vector of non-negative real integers');
end

if ~isvector(lambda) || ~isfloat(lambda) || ~isreal(lambda) || ...
    ~all(isfinite(lambda)) || ~all(lambda > 0) 
  error('LAMBDA must be a vector of positive finite real floating-point values');
end
if (numel(counts) ~= numel(lambda))
  error('COUNTS and LAMBDA must have the same number of elements');
end
if ~isscalar(timestep) || ~isreal(timestep) || ~isfinite(timestep) || ...
    ~(timestep > 0)
  error('TIMESTEP must be a positive finite real scalar');
end

L = sum(count .* log(lambda)) + log(timestep)*sum(count) - ...
    sum(log(factorial(count))) - timestep*sum(lambda);

