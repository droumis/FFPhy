function y = align_waveforms(x,Fs,t_align,t_grid)
%ALIGN_WAVEFORMS Interpolate waveforms and align at time of interest.
%
%   Y = ALIGN_WAVEFORMS(X,FS,T_ALIGN,T_GRID) takes the waveforms X sampled in
%   windows of uniform sampling rate FS and aligns them at the times specified
%   by the vector T_ALIGN, interpolated over the time vector T_GRID. This
%   function is useful when the waveforms were originally acquired in windows
%   that were not synchronously locked to an event time of interest, so that
%   different waveforms must be shifted in time by different lags to achieve
%   the desired alignment.
%
%   X must be an array of floating-point values (real or complex).
%   X((:,i,j,...,n) describes the (i,j,...,n)th waveform sampled over a
%   size(X,1)-point window at frequency FS. The (ndims(X)-1)-dimensional
%   slice X(:,:,:,...,n) is regarded as a stack of waveforms acquired
%   simultaneously in the nth sampling window (e.g. as recorded with a
%   multi-channel electrode); all such waveforms in a given stack are aligned
%   together.
%
%   FS must be a positive real floating-point scalar.
%
%   T_ALIGN must be a column vector of length size(SAMPLES,ndims(SAMPLES)) real
%   floating-point values, where each element corresponds to a
%   simultaneously-acquired stack of waveforms in X. The elements of T_ALIGN
%   are interpreted as lags relative to the time of the first ssample in the
%   window, such that the time of the first sample is 0 and the time of the
%   last sample is (size(SAMPLES,1)-1)/FS.
%
%   T_GRID must be a column vector of real floating-point values, specifying
%   the lags relative to alignment times T_ALIGN at which the waveforms are to
%   interpolated. Thus a value of zero in T_GRID corresponds to the alignment
%   time. Y is the array of interpolated waveforms, with the same
%   floating-point class and dimensionality as X. size(Y,1) == numel(T_GRID)
%   and size(X,d) == size(Y,d) for d >= 2. For elements of T_GRID which lie
%   outside of the original sample times of X, the waveform is not extrapolated
%   but rather Y is padded with NaN values.
%
%Written by SMK, 2009 September 25.
%     

if ~isnumeric(x) || ~isfloat(x) || ~(ndims(x) >= 2)
  error('X must be a non-vector array of real floating-point values');
end
if any(isnan(x(:))) || any(isinf(x(:)))
  warning('X contains one or more NaN or Inf values');
end

if ~isnumeric(Fs) || ~isscalar(Fs) || ~isreal(Fs) || ~isfloat(Fs) || ~(Fs > 0)
  error('FS must be a positive floating-point scalar');
end

if ~isnumeric(t_align) || ~isfloat(t_align) || ~isreal(t_align) || ...
    any(isnan(t_align(:))) || any(isinf(t_align(:))) || ...
    ~isvector(t_align) || (size(t_align,1) < size(t_align,2)) || ...
     (numel(t_align) ~= size(x,ndims(x)))
  error(['T_ALIGN must be a column vector whose length matches the size of ' ...
      'the last dimension of X, containing real finite non-NaN ' ...
      'floating-point values']);
end
if any(t_align < 0) || any(t_align > (size(x,1)-1)/Fs)
  warning(['elements of T_ALIGN lie outside of the range of original ' ...
      'sample times; is this what you intended?']);
end

if ~isnumeric(t_grid) || ~isfloat(t_grid) || ~isreal(t_grid) || ...
    any(isnan(t_grid(:))) || any(isinf(t_grid(:))) || ...
    ~isvector(t_grid) || (size(t_grid,1) < size(t_grid,2)) || ...
    (numel(unique(t_grid)) ~= numel(t_grid)) || ...
    ~strcmp(class(t_grid),class(t_align))
  error(['T_GRID must be a column vector of unique real finite non-NaN ' ...
      'flaoting-point values of the same class as T_ALIGN']);
end

% X and Y differ in size only in their first dimension
out_sz = size(x);
out_sz(1) = numel(t_grid);
y = nan(out_sz,class(x));
% cell array of subscripts for accessing multidimensional slices of x and y
x_subs = arrayfun(@(s)1:s,size(x),'UniformOutput',false);
y_subs = arrayfun(@(s)1:s,size(y),'UniformOutput',false);
t_original = (0:(size(x,1)-1))'/Fs;
for i = 1:numel(t_align)
  % slice along the last dimension
  x_subs{ndims(x)} = i;
  y_subs{ndims(y)} = i;
  y(y_subs{:}) = interp1(t_original - t_align(i),x(x_subs{:}), ...
      t_grid,'spline',NaN);
end


