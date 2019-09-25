function [x, phi] = x_phi_prettify(linearized,lfp,timerange,timestep,center,tol)
%X_PHI_PRETTIFY Resample trajectory defined by linearized position and LFP oscillation phase for pretty plotting.
%
%   [X,PHI] = X_PHI_PRETTIFY(LINEARIZED,LFP,TIMERANGE,TIMESTEP) takes a
%   linearized position data struct LINEARIZED, a filtered LFP data struct LFP,
%   an N-by-2 matrix of nonoverlapping uint32 timestamp intervals TIMERANGE, and
%   a real nonzero uint32 scalar TIMESTEP. LINEARIZED and LFP must be scalar
%   structs with matching meta-data fields, and the time intervals defined in
%   TIMERANGE must fall within the common range of timestamps in LINEARIZED and
%   LFP. The ''linearized_position'' field of LINEARIZED and the ''phase'' field
%   of LFP are resampled on a discrete time grid with a spacing greater than or
%   equal to TIMESTEP.
%
%   The outputs X and PHI are column vectors which contain linearized position
%   data and LFP phase data, respectively, which contain NaN padding values to
%   indicate line breaks for pretty plotting.
%
%   If TIMERANGE is a cell array, each of which contains an N-by-2 matrix of
%   nonoverlappig uint32 timestamp intervals, then X and PHI are cell arrays of
%   matching size whose contents are the resampled data for each cell of
%   TIMERANGE.
%
%   [X,PHI] = X_PHI_PRETTIFY(LINEARIZED,LFP,TIMERANGE,CENTER) takes a real
%   finite floating-point scalar argument CENTER. Phase values are mapped onto
%   the interval [CENTER-pi, CENTER+pi]. By default, CENTER = 0.
%
%   [X,PHI] = X_PHI_PRETTIFY(LINEARIZED,LFP,TIMERANGE,CENTER,TOL) takes a real
%   finite floating-point scalar argument TOL, which specifies the criterion for
%   detecting phase wrapping jumps. By default, TOL = pi. See the MATLAB
%   function UNWRAP.
%
%Depends on:
%   IS_LINEARIZED (written by SMK)
%   IS_CONTINUOUS (written by SMK)
%   IS_TIMERANGE (written by SMK)
%   STRUCT_CMP (written by SMK)
%   DIFF_INTERVALS (written by SMK)
%   INTERSECT_INTERVALS (written by SMK)
%   RESAMPLE_INTERVALS (written by SMK)
%   PRETTIFY_ANGLE_WRAP (written by SMK)
%
%Written by SMK, 2010 January 28.
%

if (exist('is_linearized') ~= 2)
  error(['X_PHI_PRETTIFY depends on m-file IS_LINEARIZED ' ...
      '(written by SMK)']);
end
if (exist('is_continuous') ~= 2)
  error(['X_PHI_PRETTIFY depends on m-file IS_CONTINUOUS ' ...
      '(written by SMK)']);
end
if (exist('is_timerange') ~= 2)
  error(['X_PHI_PRETTIFY depends on m-file IS_TIMERANGE ' ...
      '(written by SMK)']);
end
if (exist('struct_cmp') ~= 2)
  error(['X_PHI_PRETTIFY depends on m-file STRUCT_CMP ' ...
      '(written by SMK)']);
end
if (exist('diff_intervals') ~= 2)
  error(['X_PHI_PRETTIFY depends on m-file DIFF_INTERVALS ' ...
      '(written by SMK)']);
end
if (exist('intersect_intervals') ~= 2)
  error(['X_PHI_PRETTIFY depends on m-file INTERSECT_INTERVALS ' ...
      '(written by SMK)']);
end
if (exist('resample_intervals') ~= 2)
  error(['X_PHI_PRETTIFY depends on m-file RESAMPLE_INTERVALS ' ...
      '(written by SMK)']);
end
if (exist('prettify_angle_wrap') ~= 2)
  error(['X_PHI_PRETTIFY depends on m-file PRETTIFY_ANGLE_WRAP ' ...
      '(written by SMK)']);
end

if ~is_linearized(linearized) || ~isscalar(linearized)
  error('LINEARIZED must be a scalar linearized data struct');
end
if ~is_continuous(lfp) || ~isscalar(lfp) || ~isfield(lfp,'phase')
  error(['LFP must be a scalar LFP data struct with a ''phase'' field']);
end
if ~struct_cmp(linearized,lfp,{'subject','day','epoch','environment'})
  error('LINEARIZED and LFP must have concordant metadata');
end
if ~(is_timerange(timerange) || (iscell(timerange) && ...
    all(cellfun(@is_timerange,timerange))))
  error(['TIMERANGE must be a N-by-2 matrix that specifies ' ...
      'non-overlapping uint32 timestamp intervals, or a cell array of ' ...
      'such matrices']);
end
cell_flag = iscell(timerange);
if ~cell_flag
  % For computational convenience, package into a cell array. The output will
  % later be unpackaged 
  timerange = {timerange};
end
if ~all(cellfun(@(tr) isempty(diff_intervals(tr, ...
    intersect_intervals( ...
    [linearized.timestamp(1), linearized.timestamp(end)], ...
    [lfp.timestamp(1), lfp.timestamp(end)]))),timerange))
  error(['All intervals specified in TIMERANGE must fall within the ' ...
      'range of timestamps in LINEARIZED and LFP']);
end
if ~isreal(timestep) || ~isscalar(timestep) || ~isa(timestep,'uint32') || ...
    ~(timestep > 0)
  error('TIMESTEP must be a positive uint32 real scalar');
end
if (nargin < 3)
  center = 0;
end
if (nargin < 4)
  tol = pi;
end
if ~isscalar(center) || ~isreal(center) || ~isfinite(center) || ...
    ~isfloat(center)
  error('CENTER must be a real finite floating-point scalar');
end
if ~isscalar(tol) || ~isreal(tol) || ~isfinite(tol) || ~isfloat(tol)
  error('TOL must be a real finite floating-point scalar');
end

x = cell(size(timerange));
phi = cell(size(timerange));

for i = 1:numel(timerange)
  try
    [t, x_resampled] = resample_intervals(linearized.timestamp, ...
        double(linearized.linearized_position),timerange{i},timestep);
  catch
    error('Error calling RESAMPLE_INTERVALS on linearized position data');
  end
  try
    [t, phi_resampled] = resample_intervals(lfp.timestamp, ...
        unwrap(lfp.phase),timerange{i},timestep);
  catch
    error('Error calling RESAMPLE_INTERVALS on LFP phase data');
  end
  assert(isequal(size(x_resampled),size(phi_resampled)));
  x{i} = cell(size(x_resampled));
  phi{i} = cell(size(phi_resampled));
  for j = 1:numel(x_resampled)
    try
      [x{i}{j}, phi{i}{j}] = prettify_angle_wrap(x_resampled{j}, ...
          phi_resampled{j},center);
    catch
      error('Error calling PRETTIFY_ANGLE_WRAP');
    end
  end
  x{i} = cell2mat(x{i});
  phi{i} = cell2mat(phi{i});
end

if ~cell_flag
  x = x{1};
  phi = phi{1};
end


