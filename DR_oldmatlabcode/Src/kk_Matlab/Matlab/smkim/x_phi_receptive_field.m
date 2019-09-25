function ratemap = x_phi_receptive_field(unit,linearized,lfp,timerange,direction,params)
%X_PHI_RECEPTIVE_FIELD Estimate position-phase place field by local Poisson likelihod regression.
%
%   RATEMAP =
%   X_PHI_RECEPTIVE_FIELD(UNIT,LINEARIZED,LFP,TIMERANGE,DIRECTION,PARAMS)
%   produces a two-dimensional place field estimate on the cylinder defined by
%   linearized position and LPF oscillation phase, given matching single-unit
%   spike data UNIT, linearized position data LINEARIZED, filtered
%   Hilbert-transformed continuous data LFP, TIMERANGE, running DIRECTION and
%   PARAMS.
%
%   UNIT must be a scalar struct of single-unit spike data which validates with
%   IS_UNIT.
%
%   LINEARIZED must be a scalar struct of linearized position data which
%   validates with IS_LINEARIZED, whose subject, day, epoch, and environment
%   fields all match those of UNIT.
%
%   LFP must be a scalar struct of continuous data which validates with
%   IS_CONTINUOUS, whose subject, day, epoch, and environment fields all match
%   those of UNIT. Furthermore, the 'phase' field of LFP must contain phase
%   values that lie on the interval [-pi,+pi].
%   
%   TIMERANGE must be a cell vector (or 1-element cell array), where each cell
%   corresponds to a single trial. Each cell must contain a Nx2 matrix of uint32
%   timestamp intervals (of the type which validates with IS_TIMERANGE), and
%   these intervals must fall entirely within the time ranges of UNIT,
%   LINEARIZED and LFP, and intervals corresponding to different trials must not
%   overlap.
%
%   DIRECTION must assume the value either +1 or -1, indicating the direction of
%   running. This value is irrelevant to the computation, but it is included in
%   the output as meta-data. It is the caller's responsiblity to set this value
%   correctly.
%
%   PARAMS is a scalar struct with the following fields:
%     'discretization_timestep': must be a uint32 scalar, which specifies the
%         time discretization resolution. Smaller values yield better estimates
%         at greater computational cost. This parameter should be substantially
%         smaller than the period of the LFP oscillation.
%     'grid': a scalar struct with fields 'x' and 'phi'
%       grid.x is a double vector of monotonically-increasing grid points along
%         the linearized position dimension where the firing rate is to be
%         estimated, expressed in the same units as
%         LINEARIZED.position.
%       grid.phi is a double vector of monotonically-increasing grid points
%         along the LFP oscillatory phase dimension where the firing rate is to
%         be estimated, expressed in radians.
%     'halfwidth': a scalar struct with fields 'x' and 'phi'
%       halfwidth.x is a real finite positive double scalar, specifying the
%         halfwidth of the Epanechnikov kernel (for local regression) along the
%         linearized position dimension, expressed in the same units as
%         LINEARIZED.position.
%       halfwidth.phi is real positive double scalar no greater than pi,
%         specifying the halfwidth of the Epanechnikov kernel (for local
%         regression) along the LFP oscillatory phase dimension, expressed in
%         radians.
%
%   The output RATEMAP is a scalar struct array. RATEMAP has the following
%   fields:
%     'uid': inherited from UNIT
%     'region': inherited from UNIT
%     'hemisphere': inherited from UNIT
%     'subject': inherited from UNIT and LINEARIZED
%     'tetrode': inherited from UNIT
%     'depth': inherited from UNIT
%     'reference': inherited from UNIT
%     'clustnum': inherited from UNIT
%     'cluster_quality': inherited from UNIT
%     'day': inherited from UNIT and LINEARIZED
%     'epoch': inherited from UNIT and LINEARIZED
%     'environment': inherited from UNIT and LINEARIZED
%     'lfp': struct containing various metadata inherited from LFP
%     'timerange': inherited from TIMERANGE
%     'discretization_timestep': inherited from PARAMS
%     'grid': inherited from PARAMS
%     'halfwidth': inherited from PARAMS
%     'kernel': the string 'epanechnikov'
%     'position_units': inherited from LINEARIZED
%     'direction': inherited from DIRECTION
%     'spikes': a struct with fields 'timestamp', 'x', 'phi', 'trial', with the
%         time, linearized position, instantaneous LFP phase, and trial number
%         of each spike
%     'firing_rate': estimated firing rate matrix; each row correspnds to a
%         linearized position bin, and each columns corresponds to a phase bin
%     'coefficients': model coefficient estimates
%     'jackknife_coefficients': leave-one-out jackknife estimates
%     'shrinkage': parameter between 0 and 1, indicating the data-adaptive
%         shrinkage applied to the estimated coefficients of the quadratic terms
%
%   References:
%   [1] Tibshirani R., Hastie H. (1987) Local likelihood estimation. _Journal of
%       the American Statistical Association_ 82: 559-567.
%   [2] Cleveland W.S., Loader C.L. (1996) Smoothing by local regression:
%       principles and methods. In Haerdle W., Schimek M.G. (eds.), _Statistical
%       Theory and Computational Aspects of Smoothing_, pages 10-49. New York:
%       Springer.
%   [3] Itskov V., Curto C., Harris K.D. (2008) Valuations for spike train
%       prediction. _Neural Computation_ 20:644-667.
%
%Depends on:
%   IS_UNIT (written by smk)
%   IS_LINEARIZED (written by smk)
%   IS_CONTINUOUS (written by smk)
%   STRUCT_CMP (written by SMK)
%   IS_TIMERANGE (written by smk)
%   DIFF_INTERVALS (written by smk)
%   ISMEMBER_INTERVALS (written by smk)
%   TRICUBE (written by SMK)
%   MINUS_ANGLE (written by SMK)
%   INTERP1_ANGLE (written by SMK)
%
%Written by SMK, 2010 January 28.
%

TS_PER_SEC = 1e4;

if (exist('is_unit') ~= 2)
  error(['X_PHI_RECEPTIVE_FIELD depends on m-file IS_UNIT ' ...
      '(written by smk)']);
end
if (exist('is_linearized') ~= 2)
  error(['X_PHI_RECEPTIVE_FIELD depends on m-file IS_LINEARIZED ' ...
      '(written by smk)']);
end
if (exist('is_continuous') ~= 2)
  error(['X_PHI_RECEPTIVE_FIELD depends on m-file IS_CONTINUOUS ' ...
      '(written by smk)']);
end
if (exist('struct_cmp') ~= 2)
  error(['X_PHI_RECEPTIVE_FIELD depends on m-file STRUCT_CMP ' ...
      '(written by smk)']);
end
if (exist('is_timerange') ~= 2)
  error(['X_PHI_RECEPTIVE_FIELD depends on m-file IS_TIMERANGE ' ...
      '(written by smk)']);
end
if (exist('diff_intervals') ~= 2)
  error(['X_PHI_RECEPTIVE_FIELD depends on m-file DIFF_INTERVALS ' ...
      '(written by smk)']);
end
if (exist('ismember_intervals') ~= 2)
  error(['X_PHI_RECEPTIVE_FIELD depends on m-file ISMEMBER_INTERVALS ' ...
      '(written by smk)']);
end
if (exist('epanechnikov') ~= 2)
  error(['X_PHI_RECEPTIVE_FIELD depends on m-file EPANECHNIKOV ' ...
      '(written by smk)']);
end
if (exist('minus_angle') ~= 2)
  error(['X_PHI_RECEPTIVE_FIELD depends on m-file MINUS_ANGLE ' ...
      '(written by smk)']);
end
if (exist('interp1_angle') ~= 2)
  error(['X_PHI_RECEPTIVE_FIELD depends on m-file INTERP1_ANGLE ' ...
      '(written by smk)']);
end

if ~is_unit(unit) || ~isscalar(unit)
  error(['UNIT must be a scalar struct of singe-unit spike data']);
end
if ~is_linearized(linearized) || ~isscalar(linearized)
  error(['LINEARIZED must be a scalar struct of linearized ' ...
      'position data']);
end
if ~is_continuous(lfp) || ~isscalar(lfp) || ~isfield(lfp,'phase')
  error(['LFP must be a scalar struct of continuous LFP data ' ...
        'with a ''phase'' field']);
end
if ~struct_cmp(unit,linearized, ...
    {'subject','day','epoch','environment'})
  error(['UNIT and LINEARIZED must have ' ...
      'matching subject, day, epoch, environment metadata']);
end
if ~struct_cmp(unit,lfp, ...
    {'subject','day','epoch','environment'})
  error(['UNIT and LFP must have ' ...
      'matching subject, day, epoch, environment metadata']);
end

if ~iscell(timerange) || ~isvector(timerange) || ...
    ~all(cellfun(@is_timerange,timerange))
  error(['TIMERANGE must be a cell vector whose cells contain ' ...
      'Nx2 matrix of uint32 timestamp intervals.']);
end
% Check that the contents of timerange are non-overlapping
t_concat = cell2mat(timerange);
if ~is_timerange(t_concat)
  error(['Time intervals corresponding to individual trials/passes ' ...
      'must not be overlapping']);
end
if ~isempty(diff_intervals(t_concat, ...
    [linearized.timestamp(1), linearized.timestamp(end)])) || ...
    ~isempty(diff_intervals(t_concat,unit.timerange)) || ...
    ~isempty(diff_intervals(t_concat, ...
    [lfp.timestamp(1), lfp.timestamp(end)]))
  error(['All time intervals in TIMERANGE must fall within time ranges ' ...
      'of UNIT, LINEARIZED and LFP']);
end

if ~isscalar(direction) || ~isreal(direction) || ...
    ~any(direction == [+1, -1])
  error('DIRECTION must equal +1 (rightbound) or -1 (leftbound)');
end

if ~isstruct(params) || ~isscalar(params)
  error('PARAMS must be a scalar struct array');
end
if ~all(isfield(params,{'discretization_timestep','grid', ...
    'halfwidth'}))
  error('PARAMS is missing one or more required fields');
end
if ~isscalar(params.discretization_timestep) || ...
    ~isa(params.discretization_timestep,'uint32') || ...
    ~(params.discretization_timestep > 0)
  error(['''discretization_timestep'' field of PARAMS must be non-zero ' ...
      'uint32 scalar']);
end
if ~isstruct(params.grid) || ~isscalar(params.grid) || ...
    ~all(isfield(params.grid,{'x','phi'}))
  error(['''grid'' field of params must be a scalar struct with fields ' ...
      '''x'' and ''phi''']);
end
if ~isvector(params.grid.x) || ~isreal(params.grid.x) || ...
    ~isfloat(params.grid.x) || any(diff(params.grid.x) <= 0)
  error(['PARAMS.grid.x must be a real floating-point vector whose ' ...
      'elements are monotonically increasing']);
end
if ~isvector(params.grid.phi) || ~isreal(params.grid.phi) || ...
    ~isfloat(params.grid.phi) || (numel(params.grid.phi) < 3) || ...
    any((params.grid.phi < -pi) | (params.grid.phi > +pi)) || ...
    any(diff(params.grid.phi) <= 0) || ...
    (params.grid.phi(1) + 2*pi <= params.grid.phi(end))
  error(['PARAMS.grid.phi must be a real floating-point vector whose ' ...
      'elements all lie within the interval [-pi,+pi] and are ' ...
      'monotonically increasing and span less than one full cycle']);
end

if ~isstruct(params.halfwidth) || ~isscalar(params.halfwidth) || ...
    ~all(isfield(params.halfwidth,{'x','phi'}))
  error(['''halfwidth'' field of params must be a scalar struct with ' ...
      'fields ''x'' and ''phi''']);
end
if ~isscalar(params.halfwidth.x) || ~isreal(params.halfwidth.x) || ...
    ~isfinite(params.halfwidth.x) || ~(params.halfwidth.x > 0)
  error(['PARAMS.halfwidth.x must be a real positive finite scalar']);
end
if ~isscalar(params.halfwidth.phi) || ~isreal(params.halfwidth.phi) || ...
    ~isfinite(params.halfwidth.phi) || ~(params.halfwidth.phi > 0) || ...
    (params.halfwidth.phi > +pi)
  error(['PARAMS.halfwidth.phi must be a real positive finite scalar ' ...
      'no greater than +pi']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DATA PREPARATION
num_trials = 0;
binned_data = struct( ...
    'x'      , [], ...
    'phi'    , [], ...
    'counts' , [], ...
    'trial'  , [] );
% Get time, position, phase, trial number for each spike (useful for generating
% scatterplots)
spikes = struct( ...
    'timestamp' , [], ...
    'x'         , [], ...
    'phi'       , [], ...
    'trial'     , [] );
bin_edges = (unit.timerange(1):params.discretization_timestep: ...
    unit.timerange(end))';
bin_centers = 0.5*(bin_edges(1:end-1) + bin_edges(2:end));
binned_data.counts = [ binned_data.counts; ...
    histc(double(unit.timestamp),bin_edges) ];
% histc returns an extra element at the end which we don't want
binned_data.counts(end) = [];
% Resampled linearized position
binned_data.x = [ binned_data.x; interp1(double(linearized.timestamp), ...
    double(linearized.position),double(bin_centers)) ];
% Resampled LFP oscillation phase
binned_data.phi = [ binned_data.phi; interp1_angle(double(lfp.timestamp), ...
    double(lfp.phase),double(bin_centers),pi) ];
% Label trials (allocate zeros here, and populate nonzero elements in the
% loop below)
binned_data.trial = [ binned_data.trial; zeros(size(bin_centers)) ];
for j = 1:numel(timerange)
  % Trials are numbered cumulatively over *all* elements of UNIT
  num_trials = num_trials + 1;
  spiketimes = unit.timestamp(ismember_intervals(unit.timestamp, ...
      timerange{j}));
  spikes.timestamp = [ spikes.timestamp; spiketimes ];
  spikes.x = [ spikes.x; ...
      interp1(double(linearized.timestamp), ...
      double(linearized.position), ...
      double(spiketimes),'spline') ];
  spikes.phi = [ spikes.phi; ...
      interp1_angle(double(lfp.timestamp), ...
      double(lfp.phase),double(spiketimes),pi,'spline') ];
  spikes.trial = [ spikes.trial; num_trials*ones(size(spiketimes)) ];
  binned_data.trial(numel(binned_data.trial) - numel(bin_centers) + ...
      find(ismember_intervals(bin_centers,timerange{j}))) = ...
      num_trials;
end

% Exclude bins whose trial number is zero
include_idx = find(binned_data.trial > 0);
binned_data.x = binned_data.x(include_idx);
binned_data.phi = binned_data.phi(include_idx);
binned_data.counts = binned_data.counts(include_idx);
binned_data.trial = binned_data.trial(include_idx);
% Sanity checks
assert(all(isfinite(binned_data.x)));
assert(all(isfinite(binned_data.phi)));
assert(all(isfinite(binned_data.counts) & (binned_data.counts >= 0)));
assert(all(isfinite(binned_data.trial) & (binned_data.trial > 0) & ...
    (round(binned_data.trial) == binned_data.trial)));
if (num_trials < 2)
  error('Must have at least two trials for cross-validation');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SET UP FOR LOCAL LIKELIHOOD REGRESSION
% Generate plaid grid; the order of dimensions is [phi, x] because this
% guarantees that elements with adjacent linear indices lie nearby on the
% cylinder. This is convenient because we can use the estimate at one grid point
% as a starting seed for estimating at the next grid point.
[plaid.phi, plaid.x] = ndgrid(params.grid.phi,params.grid.x);
grid_sz = size(plaid.phi);
num_grid_points = prod(grid_sz);
% For each grid point, find the points that lie within bandwidth.
bandwidth_select.phi = sparse(abs(bsxfun(@minus_angle,binned_data.phi(:), ...
    params.grid.phi(:)')) < params.halfwidth.phi);
bandwidth_select.x = sparse(abs(bsxfun(@minus,binned_data.x(:), ...
    params.grid.x(:)')) < params.halfwidth.x);
% Parameter for scaling from counts/bin to count/second (this is called the
% "exposure" in the Poisson GLM literature)
timestep = double(params.discretization_timestep)/TS_PER_SEC;
log_timestep = log(timestep);
% We consider the IRLS to have converged if the estimated firing rate for every
% time bin does not change more than this amount
tol = 1e-2;
% We require at least this many iterations for burn-in before deciding that IRLS
% diverged; this is important because the algorithm can overshoot and recover
% during the first few iterations
iter_burnin = 10; 
% We quit? raise an error? if we still don't have convergence after this many
% iterations
iter_max = 25;
% Impose a minimum on counts so that log-counts do not diverge to -Inf
count_min = 1e-6;
binned_data.counts = max(binned_data.counts,count_min);
% Put bounds on the inverse link function so that the estimated firing rate
% does not diverge. log_intensity_min must be matched to count_min and timestep
% so that the IRLS algorithm remains stable
log_intensity_min = log(count_min/timestep);
log_intensity_max = log(1e+9);

% Nested function for quickly computing IRLS. We assume that predictors(:,1) is
% a column of ones. This condition is not directly tested, so be careful!
function [b, nonconv] = irls(predictors,counts,weights,offset,b_seed)
  num_obs = numel(counts);
  num_pred = size(predictors,2);
  iter = 0;
  b = b_seed;
  log_intensity_old = nan(size(counts));
  while true
    iter = iter + 1;
    % Compute expected spike count in each time bin, with bounds to prevent
    % overflow/underflow
    log_intensity = predictors*b;
    if (iter >= iter_burnin) && any(log_intensity > log_intensity_max)
      %warning('IRLS diverged');
      nonconv = true;
      break;
    end
    log_mu = offset + max(log_intensity_min,log_intensity);
    mu = exp(log_mu);
    % Working dependent variable
    z = log_mu + (counts - mu)./mu;
    % Weights
    w_sqrt = sqrt(weights.*mu);
    % Update estimate
    [Q,R] = qr(bsxfun(@times,predictors,w_sqrt),0);
    b = R \ (Q' * ((z - offset).*w_sqrt));
    % Check stopping condition
    if ~any(abs(b - b_seed) > tol) || ...
        ((iter >= iter_burnin) && ...
        ~any(abs(exp(log_intensity_old) - exp(log_intensity)) > tol))
      nonconv = false;
      break; 
    end
    log_intensity_old = log_intensity;
    if (iter >= iter_max)
      %warning('IRLS reached max iterations (%d)',iter_max);
      nonconv = true;
      break;
    end
    b_seed = b;
  end
end % end nested function irls

% Cross-validation optimization options
optim_options = optimset('TolX',1e-2,'Display','notify');
% Nested objective function which computes mixture coefficients and returns
% weighted negative log likelihood for the resulting Poisson intensity estimate.
% Note that this function depends on local_data being defined and in scope.
function obj = cv_objective(d)
  % We assume that d is on the interval [0, 1]. To improve performance, this
  % condition is not checked, so handle with care.
  % Linear combination of coefficients
  b = (1-d).*local_data.cv.linear + d.*local_data.cv.quadratic;
  % Dot product of coeffs by predictors gives the instantaneous log-intensity
  % estimate for each bin. We impose an upper bound so that the firing rate
  % estimate is not +Inf.
  mu_hat = timestep*exp(min(log_intensity_max,sum(local_data.basis.*b,2)));
  % Weighted quadratic valuation (Itskov et al., 2008). Smaller values imply
  % better fit. We use this quadratic valuation it is more robust than
  % likelihood/deviance.
  obj = nansum(local_data.weights .* (local_data.counts - mu_hat).^2) / ...
      sum(local_data.weights(~isnan(mu_hat)));
end % end nested function cv_objective

% Struct to store all fitted model coefficients. The linear coefficient vectors
% are padded with trailing zeros to match the number of coefficients in the
% quadratic model.
fits = struct( ...
    'constant'            , zeros([6, num_grid_points])           , ...
    'linear'              , zeros([6, num_grid_points])           , ...
    'linear_jackknife'    , nan([6, num_grid_points, num_trials]) , ...
    'quadratic'           , zeros([6, num_grid_points])           , ...
    'quadratic_jackknife' , nan([6, num_grid_points, num_trials]) , ...
    'shrinkage'           , zeros([1, num_grid_points])           , ...
    'linear_nonconv'      , true([1, num_grid_points])            , ...
    'quadratic_nonconv'   , true([1, num_grid_points])            );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ESTIMATE AT EACH GRID POINT
for n = 1:num_grid_points
  %disp(n)
  % Find bins that fall within bandwidth intersection
  [i,j] = ind2sub(grid_sz,n);
  idx = find(bandwidth_select.phi(:,i) & bandwidth_select.x(:,j));
  local_data.counts = binned_data.counts(idx);
  local_data.trial = binned_data.trial(idx);
  
  % Skip the estimation if there are no spikes
  if all(local_data.counts <= count_min)
    % Trivially assign the pure linear model
    fits.shrinkage(n) = 0;
    %fits.linear_nonconv(n) = true;
    %fits.quadratic_nonconv(n) = true;
    fits.constant(1,n) = log_intensity_min;
    %fits.constant(2:end,n) = 0;
    fits.linear(:,n) = fits.constant(:,n);
    fits.quadratic(:,n) = fits.constant(:,n);
    t = find(ismember(1:num_trials,unique(local_data.trial)));
    fits.linear_jackknife(1,n,t) = log_intensity_min;
    fits.linear_jackknife(2:end,n,t) = 0;
    fits.quadratic_jackknife(1,n,t) = log_intensity_min;
    fits.quadratic_jackknife(2:end,n,t) = 0;
    % Skip the rest of the computation...
    continue;
  end

  % Recenter predictors around the current grid point
  local_data.x = (binned_data.x(idx) - plaid.x(n));
  local_data.phi = minus_angle(binned_data.phi(idx),plaid.phi(n));
  % Tensor-product local weighting
  local_data.weights = ...
      epanechnikov(local_data.x ./ params.halfwidth.x) .* ...
      epanechnikov(local_data.phi ./ params.halfwidth.phi);
  % Local 2nd-order Taylor expansion. The coefficients approximate the partial
  % derivatives of the log-intensity surface.
  local_data.basis = [ ones(size(idx)), ...
      local_data.x, local_data.phi, ...
      local_data.x.*local_data.phi, ...
      0.5*local_data.x.^2, 0.5*local_data.phi.^2 ];
  % Nadaraya-Watson estimate. We use this as a seed value because it has a
  % stable, bounded value
  fits.constant(1,n) = max(log_intensity_min, ...
      log(sum(local_data.counts.*local_data.weights) / ...
      sum(local_data.weights) / timestep));

  % Local linear regression
  [fits.linear(1:3,n), fits.linear_nonconv(n)] = ...
      irls(local_data.basis(:,1:3),local_data.counts,local_data.weights, ...
      log_timestep,fits.constant(1:3,n));
  %fits.linear(4:end,n) = 0;
  if fits.linear_nonconv(n)
    % Fallback onto Nadaraya-Watson estimate
    fits.linear(:,n) = fits.constant(:,n);
    % No need to fit the quadratic model if the linear model did not converge
    fits.shrinkage(n) = 0;
    fits.quadratic_nonconv(n) = true;
    fits.quadratic(:,n) = fits.linear(:,n);
  else
    % Local quadratic regression
    [fits.quadratic(:,n), fits.quadratic_nonconv(n)] = ...
        irls(local_data.basis,local_data.counts,local_data.weights, ...
        log_timestep,fits.linear(:,n));
    % Default to linear model if quadratic model did not converge
    if fits.quadratic_nonconv(n)
      fits.quadratic(:,n) = fits.linear(:,n);
      fits.shrinkage(n) = 0;
    end
  end

  % Leave-one-out cross-validation.
  if ~fits.quadratic_nonconv(n)
    local_data.cv.linear = nan([numel(idx), 6]);
    local_data.cv.quadratic = nan([numel(idx), 6]);
  end
  for t = 1:num_trials
    % Select all bins that do not belong to this trial
    jackknife_in = find(local_data.trial ~= t);
    jackknife_out = find(local_data.trial == t);
    % No jackknife estimate if the trial did not pass through the local window
    % around the current grid point; this leaves in place the original NaN
    % values
    if isempty(jackknife_out)
      continue;
    end
    if fits.linear_nonconv(n)
      % If linear model was not fit, then linear and quadratic jackknife
      % estimates are both identical to the Nadaraya-Watson estimate
      fits.linear_jackknife(1,n,t) = max(log_intensity_min, ...
          log(sum(local_data.counts(jackknife_in) .* ...
          local_data.weights(jackknife_in)) / ...
          sum(local_data.weights(jackknife_in)) / timestep));
      fits.linear_jackknife(2:end,n,t) = 0;
      fits.quadratic_jackknife(:,n,t) = fits.linear_jackknife(:,n,t);
    elseif fits.quadratic_nonconv(n)
      % If quadratic model was not fit, then quadratic jacknife estimate is
      % identical to the linear estimate
      fits.linear_jackknife(1:3,n,t) = irls( ...
          local_data.basis(jackknife_in,1:3),local_data.counts(jackknife_in), ...
          local_data.weights(jackknife_in),log_timestep,fits.linear(1:3,n));
      fits.linear_jackknife(4:end,n,t) = 0;
      fits.quadratic_jackknife(:,n,t) = fits.linear_jackknife(:,n,t);
    else
      fits.linear_jackknife(1:3,n,t) = irls( ...
          local_data.basis(jackknife_in,1:3),local_data.counts(jackknife_in), ...
          local_data.weights(jackknife_in),log_timestep,fits.linear(1:3,n));
      fits.linear_jackknife(4:end,n,t) = 0;
      fits.quadratic_jackknife(:,n,t) = irls( ...
          local_data.basis(jackknife_in,:),local_data.counts(jackknife_in), ...
          local_data.weights(jackknife_in),log_timestep,fits.quadratic(:,n));
    end
    % Repmat so that each time bin is matched to its complementary leave-one-out
    % jackknife estimate for cross-validation
    if ~fits.quadratic_nonconv(n)
      local_data.cv.linear(jackknife_out,:) = ...
          repmat(fits.linear_jackknife(:,n,t)',size(jackknife_out));
      local_data.cv.quadratic(jackknife_out,:) = ...
          repmat(fits.quadratic_jackknife(:,n,t)',size(jackknife_out));
    end
  end
  % Shrinkage of quadratic terms to optimize cross-validation score
  if ~fits.quadratic_nonconv(n)
    fits.shrinkage(n) = fminbnd(@cv_objective,0,1,optim_options);
  end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONSTRUCT OUTPUT

%DEBUG
%{
ratemap.jackknife_linear_coefficients = ...
    permute(reshape(fits.linear_jackknife,[6 grid_sz num_trials]),[3 2 1 4]);
ratemap.jackknife_quadratic_coefficients = ...
    permute(reshape(fits.quadratic_jackknife,[6 grid_sz num_trials]),[3 2 1 4]);
ratemap.linear_coefficients = ...
    permute(reshape(fits.linear,[6 grid_sz]),[3 2 1]);
ratemap.quadratic_coefficients = ...
    permute(reshape(fits.quadratic,[6 grid_sz]),[3 2 1]);
ratemap.linear_nonconv = ...
    permute(reshape(fits.linear_nonconv,grid_sz),[2 1]);
ratemap.quadratic_nonconv = ...
    permute(reshape(fits.quadratic_nonconv,grid_sz),[2 1]);
%}

ratemap.uid = unit.uid;
ratemap.region = unit.region;
ratemap.hemisphere = unit.hemisphere;
ratemap.subject = unit.subject;
ratemap.tetrode = unit.tetrode;
ratemap.depth = unit.depth;
ratemap.reference = unit.reference;
ratemap.clustnum = unit.clustnum;
ratemap.cluster_quality = unit.cluster_quality;
ratemap.day = unit.day;
ratemap.epoch = unit.epoch;
ratemap.environment = unit.environment;
ratemap.direction = direction;
LFP_REQUIRED_FIELDS = { ...
    'subject'         , ...
    'day'             , ...
    'epoch'           , ...
    'environment'     , ...
    'electrode'       , ...
    'channel'         , ...
    'depth'           , ...
    'hemisphere'      , ...
    'region'          , ...
    'reference'       , ...
    'passband'        , ...
    'Fs'              , ...
    'sources'          };
lfp_fn = fieldnames(lfp);
ratemap.lfp = rmfield(lfp,lfp_fn(~ismember(lfp_fn,LFP_REQUIRED_FIELDS)));
ratemap.timerange = timerange;
ratemap.discretization_timestep = params.discretization_timestep;
ratemap.grid = params.grid;
ratemap.kernel = 'epanechnikov';
ratemap.halfwidth = params.halfwidth;
ratemap.position_units = linearized(1).units;
ratemap.spikes = spikes;

% Reshape so that dimensions are ordered [x, phi, coeff, trial]
ratemap.shrinkage = permute(reshape(fits.shrinkage,grid_sz),[2 1]);
% Polynomial mixing of linear and quadratic models, weighted according to the
% shrinkage parameter
ratemap.coefficients = permute(reshape( ...
    bsxfun(@times,1-fits.shrinkage,fits.linear) + ...
    bsxfun(@times,fits.shrinkage,fits.quadratic), ...
    [6 grid_sz]),[3 2 1]);
ratemap.jackknife_coefficients = permute(reshape( ...
    bsxfun(@times,1-fits.shrinkage,fits.linear_jackknife) + ...
    bsxfun(@times,fits.shrinkage,fits.quadratic_jackknife), ...
    [6, grid_sz, num_trials]),[3 2 1 4]);
ratemap.firing_rate = exp(ratemap.coefficients(:,:,1));

end % end main function (this is necessary because there are nested functions)

