function cv = spike_xcovariogram(reference,target,timerange,params)
%   CV = SPIKE_XCOVARIOGRAM(REFERENCE,TARGET,TIMERANGE,PARAMS) computes the
%   cross-covariogram between the spike trains in the single-unit spike data
%   structs REFERENCE and TARGET. REFERENCE and TARGET must be structs of the
%   same size which both validate with IS_UNIT (written by SMK) and which have
%   consistent meta-data. If REFERENCE and TARGET are identical, then this gives
%   the auto-covariogram.
%
%   TIMERANGES must be a cell vector of the same size as REFERENCE and TARGET.
%   Each cell must contain a Nx2 matrix of uint32 timestamp intervals (of the
%   type which validates with IS_TIMERANGE), and these intervals must fall
%   entirely within the time ranges of REFERENCE and TARGET. If REFERENCE and
%   TARGET are scalar structs, then TIMERANGES may be an Nx2 uint32 matrix
%   without the packaging cell. TIMERANGES specifies the times when spikes in
%   REFERENCE are used as reference triggers for populating the histogram bins.
%
%   PARAMS is a struct which specifies how the histogram is to be computed,
%   with the following fields:
%     bin_width: width of the bin, expressed in seconds. This must be a real
%       positive non-Nan, non-Inf scalar.
%     lags: vector of bin centers, expressed in seconds. An event at lag u is
%       counted in the ith bin if 
%       (u >= lags(i)-bin_width/2) & (u <lags(i)+bin_width/2)
%       Note that the bin centers do not need to be uniformly spaced; however,
%       a warning will be issued if this is the case.
%     normalize: a logical scalar, indicating whether to normalize the result
%       with respect to asymptotic normal distribution under the null hypothesis
%       that the two spike trains are independent.
%
%   The output CV is a scalar struct array with the following fields:
%     'covariogram': the (cross/auto)-covariogram
%     'lags': inherited from PARAMS
%     'bin_width': inherited from PARAMS
%     'num_reference_spikes': number of reference spikes used to compute the
%       histogram
%     'normalize': inherited from PARAMS
%     'timerange': inherited from TIMERANGES
%     'subject': inherited from REFERENCE and TARGET
%     'environment': inherited from REFERENCE and TARGET
%     'days': inherited from REFERENCE and TARGET
%     'epochs': inherited from REFERENCE and TARGET
%     'reference': struct containing various metadata inherited from REFERENCE
%     'target': struct containing various metadata inherited from TARGET
%
%   The cross-covariogram is the raw cross-correlogram minus the product of the
%   sample mean rates of the two spike trains. If the spike trains are realized
%   from a stationary bivariate point process, then the cross-covariogram is an
%   estimator of the cross-covariance function whose asymptotic distribution is
%   approximately normal with unbiased mean, although in practice this
%   asymptotic normality is approaches only when the mean intensities are high
%   (in the same way that Poisson counts are approximately normal when the
%   intensity is high). If the two processes are independent, then the estimate
%   is normally distributed with mean zero and variance equal to the product of
%   the intensities of the two processes.
%
%   References:
%   [1] Brilllinger D.R. (1992) Nerve cell spike train data analysis: a
%       progression of technique. _Journal of the American Statistical
%       Association_ 87:260-271.
%   [2] Brody C.D. (1999) Correlations without synchrony. _Neural Computation_
%       11:1537-1551.
%   [3] Siapas AG, Lubenov EV, Wilson MA. (2005) Prefrontal phase locking to
%       hippocampal theta oscillations. _Neuron_ 46:141-151.
% 
%Depends on:
%   IS_UNIT (written by SMK)
%   STRUCT_CMP (written by SMK)
%   IS_TIMERANGE (written by smk)
%   DIFF_INTERVALS (written by smk)
%   LENGTH_INTERVALS (written by SMK)
%   ISMEMBER_INTERVALS (written by smk)
%   PP_XCORRELOGRAM_MEX (written by SMK)
%
%Written 2009 October 1 by smk.
%   


TS_PER_SEC = 1e4;

if (exist('is_unit') ~= 2)
  error(['SPIKE_XCORRELOGRAM depends on m-file IS_UNIT ' ...
      '(written by smk)']);
end
if (exist('struct_cmp') ~= 2)
  error(['SPIKE_XCORRELOGRAM depends on m-file STRUCT_CMP ' ...
      '(written by smk)']);
end
if (exist('is_timerange') ~= 2)
  error(['SPIKE_XCORRELOGRAM depends on m-file IS_TIMERANGE ' ...
      '(written by smk)']);
end
if (exist('diff_intervals') ~= 2)
  error(['SPIKE_XCORRELOGRAM depends on m-file DIFF_INTERVALS ' ...
      '(written by smk)']);
end
if exist('length_intervals') ~= 2
  error(['PP_XCOVARIOGRAM depends on m-file LENGTH_INTERVALS ' ...
      '(written by smk)']);
end
if (exist('ismember_intervals') ~= 2)
  error(['SPIKE_XCORRELOGRAM depends on m-file ISMEMBER_INTERVALS ' ...
      '(written by smk)']);
end
if exist('pp_xcorrelogram_mex') ~= 3
  error(['SPIKE_XCORRELOGRAM depends on mex-file PP_XCORRELOGRAM_MEX ' ...
      '(written by smk)']);
end

if ~is_unit(reference) || ~(isscalar(reference) || isvector(reference))
  error(['REFERENCE must be a scalar or vector struct of ' ...
      'singe-unit spike data']);
end
if ~all(struct_cmp(reference(1),reference, ...
    {'subject','uid','tetrode','reference','depth','region','hemisphere', ...
    'clustnum','environment'}))
  error(['Elements of REFERENCE must correspond to the same single unit ' ...
      'and environment']);
end

if ~is_unit(target) || ~(isscalar(target) || isvector(target))
  error(['TARGET must be a scalar or vector struct of ' ...
      'singe-unit spike data']);
end
if ~all(struct_cmp(target(1),target, ...
    {'subject','uid','tetrode','reference','depth','region','hemisphere', ...
    'clustnum','environment'}))
  error(['Elements of TARGET must correspond to the same single unit ' ...
      'and environment']);
end

if is_timerange(timerange) 
  timerange = {timerange};
end
if ~iscell(timerange) || ~all(cellfun(@is_timerange,timerange))
  error(['TIMERANGES must be a Nx2 matrix of uint32 timestamp intervals, ' ...
      'or a cell array of such timestamp interval matrices']);
end

if ~isequal(size(reference),size(target),size(timerange))
  error('REFERENCE, TARGET and TIMERANGES must match in size');
end
if ~all(arrayfun(@(r,t) struct_cmp(r,t, ...
    {'subject','day','epoch','environment'}),reference,target))
  error(['Corresponding elements of REFERENCE and TARGET must have ' ...
      'matching subject, day, epoch, environment metadata']);
end
for i = 1:numel(timerange)
  if ~isempty(diff_intervals(timerange{i},reference(i).timerange)) || ...
      ~isempty(diff_intervals(timerange{i},target(i).timerange))
    error(['All time intervals in TIMERANGES must fall within time range ' ...
        'of corresponding elements of REFERENCE and TARGET']);
  end
end

if ~isstruct(params) || ~all(isfield(params,{'bin_width','lags'}))
  error('PARAMS argument is not a struct or is missing field(s)');
end
if ~isnumeric(params.bin_width) || ~isscalar(params.bin_width) || ...
    ~isreal(params.bin_width) || isnan(params.bin_width) || ...
    isinf(params.bin_width) || ~isfloat(params.bin_width) || ...
    (params.bin_width <= 0)
  error('PARAMS.bin_width must be a positive floating-point real scalar');
end
if ~isnumeric(params.lags) || ~isvector(params.lags) || ...
    ~isreal(params.lags) || any(isnan(params.lags)) || ...
    any(isinf(params.lags)) || ~isfloat(params.lags) || ...
    any(diff(params.lags) <= 0) || ...
    (numel(unique(params.lags)) ~= numel(params.lags))
  error(['PARAMS.lags must be a vector of strictly ' ...
      'monotonically-increasing floating-point values']);
end
if any(abs(diff(params.lags) - params.bin_width) > ...
    sqrt(eps(params.bin_width)))
  warning('PARAMS.lags are spaced apart greater than PARAMS.bin_width');
end
% coerce lags to be a row vector
if size(params.lags,1) > size(params.lags,2)
  params.lags = params.lags';
end
if ~islogical(params.normalize) || ~isscalar(params.normalize)
  error('PARAMS.normalize must be a logical scalar');
end

raw_cc = zeros(size(params.lags));
num_reference_spikes = 0;
num_target_spikes = 0;
total_duration = 0;
for i = 1:numel(reference)
  % Select only those reference spikes which occur within specified time
  % interval(s)
  reference_times = double(reference(i).timestamp( ...
      ismember_intervals(reference(i).timestamp,timerange{i})))/TS_PER_SEC;
  target_times = double(target(i).timestamp)/TS_PER_SEC;

  % Accumulate counts in histogram bins (mex function)
  raw_cc = raw_cc + pp_xcorrelogram_mex(reference_times,target_times, ...
      double(params.lags),double(params.bin_width));

  % Update cumulative count of spikes that occur within the specified time
  % intervals and the cumulative duration of these time intervals
  num_reference_spikes = num_reference_spikes + numel(reference_times);
  num_target_spikes = num_target_spikes + ...
      nnz(ismember_intervals(target(i).timestamp,timerange{i}));
  total_duration = total_duration + ...
      sum(double(length_intervals(timerange{i}))/TS_PER_SEC);
end

% total duration
mean_rate1 = num_reference_spikes / total_duration;
mean_rate2 = num_target_spikes / total_duration;
covariogram = (raw_cc / params.bin_width / total_duration) - ...
    (mean_rate1 * mean_rate2);
if params.normalize
  var_est = mean_rate1 * mean_rate2 / total_duration / params.bin_width;
  covariogram = covariogram ./ sqrt(var_est);
end

% Construct output struct
cv.covariogram = covariogram;
cv.normalize = params.normalize;
cv.lags = params.lags;
cv.bin_width = params.bin_width;
cv.timerange = timerange;
cv.subject = reference(1).subject;
cv.environment = reference(1).environment;
cv.days = {reference(:).day}';
cv.epochs = {reference(:).epoch}';
UNIT_REQUIRED_FIELDS = { ...
    'uid'                     , ...
    'subject'                 , ...
    'day'                     , ...
    'epoch'                   , ...
    'environment'             , ...
    'timerange'               , ...
    'tetrode'                 , ...
    'depth'                   , ...
    'hemisphere'              , ...
    'region'                  , ...
    'reference'               , ...
    'passbands'               , ...
    'thresholds'              , ...
    'clustnum'                , ...
    'cluster_quality'         , ...
    'sources'                 };
reference_fn = fieldnames(reference);
cv.reference = rmfield(reference, ...
    reference_fn(~ismember(reference_fn,UNIT_REQUIRED_FIELDS)));
target_fn = fieldnames(target);
cv.target = rmfield(target, ...
    target_fn(~ismember(target_fn,UNIT_REQUIRED_FIELDS)));

