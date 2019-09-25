function [occ,d] = density_map(x,y,Fs,params)
%DENSITY_MAP Multidimensional kernel density estimate.
%
%   [OCC,D] = DENSITY_MAP(X,Y,FS,PARAMS) computes approximate kernel
%   density estimates of events that occur upon a multidimensional space (useful
%   for estimating neural receptive fields from spike trains and stimulus). 
%
%   X is an M-by-N matrix which describes a real-valued time-varying
%   N-dimensional stimulus that is sampled at a uniform sampling rate. X(i,:) is
%   the ith discrete sample of a trajectory through stimulus space.
%
%   Y is an M-by-P matrix of event counts (integer values) in bins sampled at
%   the same sampling rate as X. Each column of Y represents a stream of events,
%   for which a separate density map will be estimated. X and Y must have the
%   same length along their first dimension. 
%
%   FS is the (scalar) sampling rate of X and Y; 1/FS is the time step between
%   consecutive rows in X and Y. This parameter is necessary so that the output
%   D is correctly expressed as the rate of events per unit time. The default
%   value is 1.
%
%   PARAMS is a struct array of length N, with the following fields:
%     'periodic': either true or false
%     'grid': a vector of evenly-spaced grid points at which the kernel
%         density is to be estimated; if 'periodic' is true, then grid(1)
%         coincides with grid(end)
%     'kernel': a string identifying the smoothing kernel. Available options are
%         'epanechnikov', 'gaussian' (only if 'periodic' is false), 'vonmises'
%         (only if 'periodic' is true).
%     'dispersion': a real positive scalar that parametrizes the smoothing width
%         of the kernel. For 'epanechnikov', this is the half-width of the
%         kernel support, expressed in grid units. For 'gaussian', this is the
%         standard deviation, expressed in grid units. For 'vonmises', this is
%         the concentration parameter (kappa) of a von Mises distribution on the
%         periodic interval [grid(1) grid(end)].If the empty vector [] is given,
%         then no smoothing is done.
%
%   OCC is a N-dimensional array of size CELLFUN(@LENGTH,{PARAMS(:).GRID})-1,
%   which is the estimated occupancy density of the stimulus.
%
%   D is an (N+1)-dimensional array of size
%   [cellfun(@length,{PARAMS(:).GRID})-1, P] which contains the density
%   estimates corresponding to the event streams Y(:,1),Y(:,2),...,Y(:,P).
%   Because of computational performance, this function computes a discretized
%   *approximation* over voxels of stimulus space, using a discrete (*not*
%   continuous) smoothing kernel. This approximation comes close to the true
%   kernel density estimate only when the grid is finely spaced.
%   
%Depends on:
%   HISTCN (written by smk)
%   VONMISESPDF (written by smk)
%   EPANECHNIKOV (written by smk)
%   NORMPDF (MATLAB Statistics Toolbox)
%   RANGE (MATLAB Statistics Toolbox)
%
%Written by smk, 6 November 2008.
%

if (exist('histcn') ~= 2)
  error('DENSITY_MAP depends on m-file HISTCN (written by smk)');
end
if (exist('epanechnikov') ~= 2)
  error('DENSITY_MAP depends on m-file EPANECHNIKOV (written by smk)');
end
if (exist('vonmisespdf') ~= 2)
  error('DENSITY_MAP depends on m-file VONMISESPDF (written by smk)');
end
if (exist('normpdf') ~= 2)
  error('DENSITY_MAP depends on m-file NORMPDF (MATLAB Statistics toolbox');
end
if (exist('range') ~= 2)
  error('DENSITY_MAP depends on m-file RANGE (MATLAB Statistics toolbox');
end

if ~isnumeric(x(:)) || ~isreal(x(:)) || ~all(isfinite(x(:)))
  error('X must contain only real finite numbers');
end
if ~isnumeric(y) || ~isreal(y) || ~all(isfinite(y(:))) || ...
    ~all(y(:) == round(y(:))) || ~all(y(:) >= 0)
  error('All elements of Y must be non-negative real integers');
end
if (ndims(x) ~= 2) || (ndims(y) ~= 2) || (size(x,1) ~= size(y,1))
  error('X and Y must be matrices with the same number of rows');
end
% M is the number of samples/time bins
% N is the dimensionality of the stimulus X
% P is the number of concurrent event processes
[m, n] = size(x);
[m, p] = size(y);
if ~isstruct(params) || (length(params) ~= n)
  error('PARAMS must be a struct array of length SIZE(X,2)');
end
if ~all(isfield(params,{'periodic','grid','kernel','dispersion'}))
  error('PARAMS is missing one or more required fields');
end

for i = 1:n
  if ~isscalar(params(i).periodic) || ~islogical(params(i).periodic)
    error('PARAMS.periodic must be a logical scalar (true or false)');
  end
  if ~isnumeric(params(i).grid) || ~isvector(params(i).grid) || ...
    (numel(params(i).grid) < 3)
    % Avoid collapse of singleton dimensions of N-dimensional histogram. 3 grid
    % points corresponds to 2 histogram bins, which avoids the singleton problem
    error('PARAMS must contain numeric grid vectors with at least 3 elements');
  end
  if ~all(x(:,i) <= params(i).grid(end)) || ~all(x(:,i) >= params(i).grid(1)) 
    warning('some values of X are outside the bounds defined by PARAMS');
    % Wrap values if the support is periodic
    if params(i).periodic
      warning('Values will be wrapped along the periodic %dth dimension',i);
      x(:,i) = mod(double(x(:,i)) - params(i).grid(1), ...
          params(i).grid(end) - params(i).grid(1)) + params(i).grid(1);
    end
  end
  if any(diff(params(i).grid) <= 0) || ...
      any(abs(diff(diff(params(i).grid))) > sqrt(max(eps(params(i).grid))))
    error('grid spacings in PARAMS must be uniform monotonic increasing');
  end
  if ~ischar(params(i).kernel) || ~any(strcmp(params(i).kernel, ...
      {'epanechnikov','gaussian','vonmises'}))
    error('PARAMS.kernel is not a recognized string value');
  end
  if params(i).periodic && strcmp(params(i).kernel,'gaussian')
    error('''gaussian'' kernel can only be used on a non-periodic support');
  end
  if ~params(i).periodic && strcmp(params(i).kernel,'vonmises')
    error('''vonmises'' kernel can only be used on a periodic support');
  end
  if ~(isscalar(params(i).dispersion) && isreal(params(i).dispersion) && ...
      isfinite(params(i).dispersion) && (params(i).dispersion > 0)) && ...
      ~isequal(params(i).dispersion,[])
    error(['PARAMS dispersion field must be either a real positive finite ' ...
      'scalar or the empty matrix []']);
  end
  if params(i).periodic && strcmp(params(i).kernel,'epanechnikov') && ...
      (params(i).dispersion >= (params(i).grid(end) - params(i).grid(1))/2)
    error(['halfwidth of Epanechnikov kernel must be less than one half ' ...
        'cycle of the periodic support']);
  end
end
if nargin==3
  Fs = 1;
end
% define the grid
edges = {params(:).grid};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mini-tutorial: How does this function work? Kernel density estimation replaces
% each data point (which can be represented as a Dirac delta function) with a
% smooth "bump" of density, and then summing these bumps over all data points.
% Equivalently, one can take the sum of Dirac delta functions at the data points
% and then convolve with the smoothing kernel. Here, the Dirac delta function is
% approximated by a voxel indicator function, and the smoothing kernel is
% discrete, not continuous. This approximation approaches the kernel density
% estimate as the grid spacing is made finer.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Construct a discrete approximation to the sum of Dirac delta indicators of the
% data points in stimulus space
[xcounts, bin] = histcn(x,edges{:});
nzsubs = find(all(bin,2));
% For each column of Y, compute spike counts in the same bins as xcounts
for j = 1:p
  ycounts{j} = accumarray(bin(nzsubs,:),y(nzsubs,j),size(xcounts),@sum,0);
end

% Wait, something isn't quite ready! HISTCN returns extra counts along the
% ends. For each non-periodic dimension, we can simply throw away these extra
% values. For each periodic dimension, we need to add these counts to the
% first bin (edge wrap-around).
for i = 1:n
  % keep track of the current size of the empirical density matrix. it's tricky
  % because it changes with every iteration of this for loop!
  sz = size(xcounts);
  [subs{1:n}] = ind2sub(sz,1:prod(sz));
  % find the indices corresponding to the "extra" elements at the end of this
  % dimension 
  lastbinidx = (subs{i} == length(edges{i}));
  % new size which we want to enforce
  sz(i) = sz(i) - 1;
  if params(i).periodic
    firstbinidx = (subs{i} == 1);
    xcounts(firstbinidx) = xcounts(firstbinidx) + xcounts(lastbinidx);
    xcounts(lastbinidx) = [];
    xcounts = reshape(xcounts,sz);
    for j = 1:p
        ycounts{j}(firstbinidx) = ycounts{j}(firstbinidx) + ...
            ycounts{j}(lastbinidx);
        ycounts{j}(lastbinidx) = [];
        ycounts{j} = reshape(ycounts{j},sz);
    end
  else
    xcounts(lastbinidx) = [];
    xcounts = reshape(xcounts,sz);
    for j = 1:p
        ycounts{j}(lastbinidx) = [];
        ycounts{j} = reshape(ycounts{j},sz);
    end
  end
end
% sanity check
sz = size(xcounts);
if ~all(cellfun(@(c) isequal(sz,size(c)),ycounts))
  error('bug in this function');
end

% The multidimensional smoothing kernel is a tensor product of orthogonal
% components.
for i = 1:n
  switch params(i).kernel
  case 'epanechnikov'
    if ~isempty(params(i).dispersion)
      % Kernel spans from -halfwidth to +halfwidth
      hw = ceil(params(i).dispersion/mean(diff(params(i).grid)));
      kernsz(i) = 1 + 2*hw;
      temp{i} = epanechnikov(linspace(-hw,+hw,kernsz(i))' * ...
          mean(diff(params(i).grid))/params(i).dispersion) / ...
          mean(diff(params(i).grid));
    else
      % These artificial kernel components are designed to work around the
      % annoying automatic suppression of singleton dimensions in MATLAB
      kernsz(i) = 3;
      temp{i} = [0 1 0];
    end
  case 'gaussian'
    assert(~params(i).periodic);
    if ~isempty(params(i).dispersion)
      % Kernel spans from -4*sigma to +4*sigma
      hw = ceil(4*params(i).dispersion/mean(diff(params(i).grid)));
      kernsz(i) = 1 + 2*hw;
      temp{i} = normpdf(linspace(-hw,+hw,kernsz(i))',0, ...
          params(i).dispersion/mean(diff(params(i).grid)));
    else
      % These artificial kernel components are designed to work around the
      % annoying automatic suppression of singleton dimensions in MATLAB
      kernsz(i) = 3;
      temp{i} = [0 1 0];
    end
  case 'vonmises'
    assert(params(i).periodic);
    if ~isempty(params(i).dispersion)
      % along each circular dimension, kernel extends over a full cycle
      kernsz(i) = sz(i);
      temp{i} = vonmisespdf(linspace(-pi,+pi,kernsz(i))',0,params(i).dispersion);
    else
      % These artificial kernel components are designed to work around the
      % annoying automatic suppression of singleton dimensions in MATLAB
      kernsz(i) = 3;
      temp{i} = [0 1 0];
    end
  otherwise
    error('kernel type is not recognized');
  end
end
% This test condition works around the annoying fact that every vector in MATLAB
% is really a 2-dimensional array
if (n > 1)
  [kerncomp{1:n}] = ndgrid(temp{:});
  % multiply these orthogonal components together
  kernel = ones(kernsz);
  for i = 1:n
    kernel = kernel .* kerncomp{i};
  end
else
  % If there is only one dimension, we want to make sure that kernel is a vector
  kernel = temp{1};
  assert(isvector(kernel));
end
% normalize kernel
kernel = kernel ./ sum(kernel(:));

% Repeatedly tile each empirical density function along its circular dimensions
% in order to correctly handle wrap-around at the edges. This excess will be
% trimmed at a later step.
repsz = ones([1 n]);
repsz([params(:).periodic]) = 3;
% Indices to recover the original array, which will be the middle third in
% between flanking repeats
[subs{1:n}] = ind2sub(sz,1:prod(sz));
for i = 1:n
  if repsz(i)==3
      subs{i} = subs{i} + sz(i);
  end
end
middle_idx = sub2ind(sz.*repsz,subs{:});
% Convolve the tiled empirical density functions with the smoothing kernel to
% obtain approximate kernel density estimate; note that 
% size(xsmooth) == sz .* repsz
xsmooth = convn(repmat(xcounts,repsz),kernel,'same');
% Recover the middle third section
occ = reshape(xsmooth(middle_idx),sz)/Fs;
d = nan([size(occ) size(y,2)]);
leading_subs = arrayfun(@(i) 1:i,size(occ),'UniformOutput',false);
for j = 1:p
  ysmooth{j} = convn(repmat(ycounts{j},repsz),kernel,'same');
  % remember that size(ysmooth{j}) = sz .* repsz
  d(leading_subs{:},j) = reshape(ysmooth{j}(middle_idx),sz) ./ occ;
end

