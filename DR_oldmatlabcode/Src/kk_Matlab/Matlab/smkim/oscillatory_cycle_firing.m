function [rate, mean_phase] = oscillatory_cycle_firing(unit,lfp,ti,phase_range)
%OSCILLATORY_CYCLE_FIRING Compute firing rate per cycle(s) of LFP oscillation
%
%   RATE = OSCILLATORY_CYCLE_FIRING(UNIT,LFP,TI,PHASE_RANGE) computes
%   firing rate in a boxcar window whose variable size is defined in terms of
%   adjacent oscillatry cycles. UNIT is a scalar struct containing single-unit
%   spike data. LFP is a scalar struct containing bandpass-filtered LFP data
%   with a valid ''phase'' field, whose meta-data match that of UNIT. TI is an
%   array of real uint32 timestamps at which the estimate is to be computed,
%   such that 
%
%       all(ismember_intervals(TI,UNIT.timerange)) && ...
%       all(ismember_intervals(TI,LFP.timestamp([1 end])))
%
%   PHASE_RANGE is a 2-element row vector of real finite floating-point values
%   such that 
%
%       sign(PHASE_RANGE(1)) < sign(PHASE_RANGE(2))
%
%   For each element of TI, this function counts all spikes within PHASE_RANGE
%   radians (as unwrapped phase) and divides by the duration of this window. If
%   this window is truncated at the start or end of the LFP record, the part of
%   the window that is available will be used. Intervals on which phase is not
%   strictly increasing will be censored; if the filtered LFP has many phase
%   slips, this will give poor results.
%
%   [RATE, MEAN_PHASE] = OSCILLATORY_CYCLE_FIRING(...) also computes the
%   circular mean LFP oscillation phase in each window. If there are no spikes
%   in a window, the corresponding phase value is NaN.
%
%   References:
%   [1] Harris K.D., Henze D.A., Hirase H., Leinekugel X., Dragoi G., Czurko A.,
%       Buzsaki G. (2002) Spike train dynamics predicts theta-related phase
%       precession in hippocampal pyramidal cells. _Nature_ 417: 738-741.
%   [2] Huxter J., Burgess N., O'Keefe J. (2003) Independent rate and temporal
%       coding in hippocampal pyramidal cells. _Nature_ 425: 838-832.
%
%Depends on:
%   IS_UNIT (written by SMK)
%   IS_CONTINUOUS (written by SMK)
%   STRUCT_CMP (written by SMK)
%   ISMEMBER_INTERVALS (written by SMK)
%   DIFF_INTERVALS (written by SMK)
%   MEAN_ANGLE (written by SMK)
%   INTERP1_ANGLE (written by SMK)
%
%Written by SMK, 2010 January 18.
%

TS_PER_SEC = 1e4;

if (exist('is_unit') ~= 2)
  error(['OSCILLATORY_CYCLE_FIRING depends on m-file IS_UNIT ' ...
      '(written by smk)']);
end
if (exist('is_continuous') ~= 2)
  error(['OSCILLATORY_CYCLE_FIRING depends on m-file IS_CONTINUOUS ' ...
      '(written by smk)']);
end
if (exist('struct_cmp') ~= 2)
  error(['OSCILLATORY_CYCLE_FIRING depends on m-file STRUCT_CMP ' ...
      '(written by smk)']);
end
if (exist('ismember_intervals') ~= 2)
  error(['OSCILLATORY_CYCLE_FIRING depends on m-file ISMEMBER_INTERVALS ' ...
      '(written by smk)']);
end
if (exist('diff_intervals') ~= 2)
  error(['OSCILLATORY_CYCLE_FIRING depends on m-file DIFF_INTERVALS ' ...
      '(written by smk)']);
end
if (exist('mean_angle') ~= 2)
  error(['OSCILLATORY_CYCLE_FIRING depends on m-file MEAN_ANGLE ' ...
      '(written by smk)']);
end
if (exist('interp1_angle') ~= 2)
  error(['OSCILLATORY_CYCLE_FIRING depends on m-file INTERP1_ANGLE ' ...
      '(written by smk)']);
end

if ~is_unit(unit) || ~isscalar(unit)
  error(['UNIT must be a scalar struct of singe-unit spike data']);
end
if ~is_continuous(lfp) || ~isscalar(lfp) || ~isfield(lfp,'phase')
  error(['LFP must be a scalar struct of continuous LFP data ' ...  
      'with a ''phase'' field']);
end
if ~struct_cmp(unit,lfp,{'subject','day','epoch','environment','hemisphere'})
  error(['UNIT and LFP must have matching subject, day, epoch, ' ...
      'environment, hemisphere metadata']);
end
if ~struct_cmp(unit,lfp,{'subject','day','epoch','environment','hemisphere'})
  error(['UNIT and LFP must have matching subject, day, epoch, ' ...
      'environment, hemisphere metadata']);
end
if ~isempty(diff_intervals(unit.timerange, ...
    [lfp.timestamp(1), lfp.timestamp(end)]))
  error('UNIT.timerange must fall within the range of timestamps in LFP');
end

if ~isreal(ti) || ~isa(ti,'uint32')
  error('TI must be an array of real uint32 timestamps');
end
if ~all(ismember_intervals(ti(:),unit.timerange))
  error('All elements of TI must fall within UNIT.timerange');
end

if ~isvector(phase_range) || ~isequal(size(phase_range),[1 2]) || ...
    ~isreal(phase_range) || ~isfloat(phase_range) || ...
    ~all(isfinite(phase_range(:))) || (diff(phase_range) <= 0)
  error(['PHASE_RANGE must be a 2-element row vector that defines a real ' ...
      'finite floating-point interval']);
end

% Remove all negative phase differences (also, note that t_censored is cast as
% double this allows floating-point computations)
t_censored = double(lfp.timestamp);
p_censored = unwrap(double(lfp.phase));
while true
  exclude_flag = (diff(p_censored) <= eps(pi));
  if ~any(exclude_flag)
    break;
  else
    t_censored([false; exclude_flag]) = [];
    p_censored([false; exclude_flag]) = [];
 end
end

% Interpolate desired phase interval around each time point in TI. We use the
% original, *uncensored* unwrapped LFP phase
p_target = bsxfun(@plus,double(phase_range), ...
    interp1(double(lfp.timestamp),unwrap(lfp.phase),double(ti(:))));
assert(isequal(size(p_target),[numel(ti), 2]));

% Now, for each element ti(j), we look for time interval when the unwrapped
% phase spans (p_target(j) + phase_range). 
t_ranges = interp1(p_censored,t_censored,p_target);
% Interpolation will yield NaN values where the target phase exceeds the
% unwrapped phase range. We replace these NaN values with nearest-neighbor
% extrapolation, which is guaranteed to be either t_censored(1) or
% t_censored(end) because monotonicity was enforced by censoring.
nan_flag = isnan(t_ranges);
t_ranges(nan_flag) = interp1(p_censored,t_censored,p_target(nan_flag), ...
    'nearest','extrap');
% Get durations of these time intervals
durations = diff(t_ranges,1,2) ./ TS_PER_SEC;
% Reformat t_ranges to be a cell column vector in which each cell contains a
% uint32 timestamp interval
t_ranges = mat2cell(uint32(t_ranges),ones([numel(ti) 1]),2);

% Cell array of lookup indices for spikes that fall within each window.
spike_idx = cellfun(@(arg) find(ismember_intervals(unit.timestamp,arg)), ...
    t_ranges,'UniformOutput',false);

rate = cellfun(@numel,spike_idx,'UniformOutput',true) ./ durations;
rate = reshape(rate,size(ti));

if (nargout >= 2)
  mean_phase = nan(size(ti));
  spike_phases = interp1_angle(double(lfp.timestamp),lfp.phase, ...
      double(ti),pi);
  assert(all(isfinite(spike_phases)));
  for i = 1:numel(ti)
    if ~isempty(spike_idx{i})
      mean_phase(i) = mean_angle(spike_phases(spike_idx{i}));
    end
  end
end

