function m = running_median(t,x,hw)
%RUNNING_MEDIAN Compute running median of a time series.
%
%   M = RUNNING_MEDIAN(T,X,HW) take a time series X sampled at T and computes a
%   running median M whose half-width is HW. T must be a vector of unique real
%   finite non-NaN double values. X must be a numeric array such that 
%   size(X,1) == numel(T). HW must be a positive real non-NaN finite double
%   scalar. 
%
%   The output M is an array of the same size and numeric type as X. M(i,:)
%   contains the median of X(j,:) for j such that (T(j) >= T(i) - HW) & (T(j) <
%   T(i) - HW). Note the cadlag assymmetry of this inclusion condition; if T is
%   a vector of evenly-spaced integers, you may wish to add a small fudge factor
%   in HW to include the rightmost sample.
%
%   NaN elements in X are ignored in the manner of NANMEDIAN (MATLAB Statistics
%   Toolbox).
%
%   References: 
%
%   Tukey J.W. (1977) _Exploratory Data Analysis_. Addison-Wesley.
%
%   Gallagher N.C., Wise G.L. (1981)  A theoretical analysis of the properties
%   of median filters. _IEEE Transactions on Acoustics, Speech, and Signal
%   Processing_ 29:1136-1141.
%
%Depends on:
%   FIND_NEARBY_MEX (written by SMK)
%   NANMEDIAN (MATLAB Statistics Toolbox)
%
%Written by SMK, 2009 November 19.
%
%!!!TODO: Generalize this to a block median filter (where T is no longer a
%vector but can be a multidimensional array)
%

if (exist('find_nearby_mex') ~= 3)
  error(['RUNNING_MEDIAN depends on the mex-file FIND_NEARBY_MEX ' ...
      '(written by SMK)']);
end
if (exist('nanmedian') ~= 2)
  error(['RUNNING_MEDIAN depends on the m-file NANMEDIAN ' ...
      '(MATLAB Statistics Toolbox)']);
end

if ~isvector(t) || ~isa(t,'double') || ~isreal(t) || ~all(isfinite(t)) || ...
    (numel(unique(t)) ~= numel(t))
  error('T must be a vector of unique real finite non-NaN double values');
end
if ~isnumeric(x) || isempty(x) || (size(x,1) ~= numel(t))
  error(['X must be a non-empty array whose length along the first ' ...
      'dimension equals the length of T']);
end
if ~isscalar(hw) || ~isa(hw,'double') || ~isreal(hw) || ~isfinite(hw) || ...
    ~(hw > 0)
  error('HW must be a positive real non-NaN finite double scalar');
end

% Force t to be a column vector
if (size(t,1) < size(t,2))
  t = t';
end

% j is a cell array, such that m(i,:) = nanmedian(x(j{i},:)
j = find_nearby_mex(t,t,-hw,+hw);

% Determine the size of each "row" of X
sz = size(x);
slice_ind = arrayfun(@(n) 1:n,sz(2:end),'UniformOutput',false);

% Compute median in each row and concatenate
try
  m = cell2mat(cellfun(@(i) nanmedian(x(i,slice_ind{:}),1),j, ...
      'UniformOutput',false));
  assert(isequal(size(m),size(x)) & isequal(class(m),class(x)));
catch
  error('Something went wrong. Sorry!');
end


