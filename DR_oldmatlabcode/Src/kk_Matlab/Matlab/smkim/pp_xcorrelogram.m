function cc = pp_xcorrelogram(reference,target,timerange,params)
%PP_XCORRELOGRAM Compute cross-correlogram of a bivariate point process realization.
%
%   CC = PP_XCORRELOGRAM(REFERENCE,TARGET,TIMERANGE,PARAMS) takes two
%   vectors of event times, TARGET and REFERENCE, which are sampled during the
%   uninterrupted time interval specified by TIMERANGE.
%
%   REFERENCE and TARGET must be vectors of monotonically-increasing real
%   uint32 timestamps (1e-4 seconds per timestamp unit). REFERENCE must be
%   non-empty. If REFERENCE and TARGET are identical, then this gives the
%   auto-correlogram.
%
%   TIMERANGE must be a 1x2 row vector of monotonically-increasing uint32
%   timestamps [start, end], such that all elements of TARGET and REFERENCE fall
%   within TIMERANGE. An error will be raised if TARGET or REFERENCE contain
%   times that lie outside of the intervals specified by TIMERANGE.
%
%   CC is the histogram of TARGET events at lags relative to REFERENCE events.
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
%     normalize: logical scalar, indicating whether to normalize counts by
%       (numel(REFERENCE) * bin_width)
%
%   References:
%    
%   Brilllinger D.R. (1992) Nerve cell spike train data analysis: a progression
%   of technique. _Journal of the American Statistical Association_ 87:260-271.
%
%   Brody C.D. (1999) Correlations without synchrony. _Neural Computation_
%   11:1537-1551.
% 
%Depends on:
%   IS_TIMERANGE (written by SMK)
%   ISMEMBER_INTERVALS (written by SMK)
%   PP_XCORRELOGRAM_MEX (written by SMK)
%
%Written 2009 October 1 by smk.
%   

TS_PER_SEC = 1e4;

if exist('is_timerange') ~= 2
  error('PP_XCORR depends on m-file IS_TIMERANGE (written by smk)');
end
if exist('ismember_intervals') ~= 2
  error('PP_XCORR depends on m-file ISMEMBER_INTERVALS (written by smk)');
end
if exist('pp_xcorrelogram_mex') ~= 3
  error(['PP_XCORR depends on mex-file PP_XCORRELOGRAM_MEX ' ...
      '(written by smk)']);
end

if ~isnumeric(reference) || ~isvector(reference) || ...
    ~isa(reference,'uint32') || ~isreal(reference) || ...
    isempty(reference) || any(diff(reference) < 0)
  error(['REFERENCE must be a non-empty monotonically-increasing vector ' ...
      'of uint32 timestamps']);
end
if ~isnumeric(target) || ~isvector(target) || ...
    ~isa(target,'uint32') || ~isreal(target) || any(diff(target) < 0)
  error(['TARGET must be a monotonically-increasing vector of ' ...
      'uint32 timestamps']);
end
if ~is_timerange(timerange) || (size(timerange,1) ~= 1)
  error(['TIMERANGE must be a 1-by-2 vector of monotonically-increasing ' ...
      'uint32 timestamps, corresponding to a single uninterrupted time ' ...
      'interval']);
end
if ~all(ismember_intervals(target,timerange))
  error('event times in TARGET must lie within TIMERANGE');
end
if ~all(ismember_intervals(reference,timerange))
  error('event times in REFERENCE must lie within TIMERANGE');
end
if ~isstruct(params) || ~all(isfield(params,{'bin_width','lags','normalize'}))
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
if any(diff(params.lags) > params.bin_width)
  warning('PARAMS.lags are spaced apart greater than PARAMS.bin_width');
end
if ~islogical(params.normalize) || ~isscalar(params.normalize)
  error('PARAMS.normalize must be a logical scalar');
end

% coerce event time vectors to be column vectors
if size(target,1) < size(target,2)
  target = target';
end
if size(reference,1) < size(reference,2)
  reference = reference;
end
% coerce lags to be a row vector
if size(params.lags,1) > size(params.lags,2)
  params.lags = params.lags';
end
% call the mex function
cc = pp_xcorrelogram_mex( ...
    double(reference)/TS_PER_SEC,double(target)/TS_PER_SEC, ...
    double(params.lags),double(params.bin_width));
% normalize by the number of reference events and bin_width
if params.normalize
  cc = cc / numel(reference) / params.bin_width;
end

