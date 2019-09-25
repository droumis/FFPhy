function [fout0,xout,u]=ksdensity(yData,varargin)
%KSDENSITY Compute kernel density or distribution estimate
%   [F,XI]=KSDENSITY(X) computes a probability density estimate of the sample
%   in the vector X.  KSDENSITY evaluates the density estimate at 100 points
%   covering the range of the data.  F is the vector of density values and XI
%   is the set of 100 points.  The estimate is based on a normal kernel
%   function, using a window parameter (bandwidth) that is a function of the
%   number of points in X.
%
%   F=KSDENSITY(X,XI) specifies the vector XI of values where the density
%   estimate is to be evaluated.
%
%   [F,XI,U]=KSDENSITY(...) also returns the bandwidth of the kernel smoothing
%   window.
%
%   KSDENSITY(...) without output arguments produces a plot of the results.
%
%   KSDENSITY(AX,...) plots into axes AX instead of GCA.
%
%   [...]=KSDENSITY(...,'PARAM1',val1,'PARAM2',val2,...) specifies parameter
%   name/value pairs to control the density estimation.  Valid parameters
%   are the following:
%
%      Parameter    Value
%      'censoring'  A logical vector of the same length of X, indicating which
%                   entries are censoring times (default is no censoring).
%      'kernel'     The type of kernel smoother to use, chosen from among
%                   'normal' (default), 'box', 'triangle', and
%                   'epanechnikov'.
%      'npoints'    The number of equally-spaced points in XI.
%      'support'    Either 'unbounded' (default) if the density can extend
%                   over the whole real line, or 'positive' to restrict it to
%                   positive values, or a two-element vector giving finite
%                   lower and upper limits for the support of the density.
%      'weights'    Vector of the same length as X, giving the weight to
%                   assign to each X value (default is equal weights).
%      'width'      The bandwidth of the kernel smoothing window.  The default
%                   is optimal for estimating normal densities, but you
%                   may want to choose a smaller value to reveal features
%                   such as multiple modes.
%      'function'   The function type to estimate, chosen from among 'pdf',
%                   'cdf', 'icdf', 'survivor', or 'cumhazard' for the density,
%                   cumulative probability, inverse cumulative probability,
%                   survivor, or cumulative hazard functions, respectively.
%
%   In place of the kernel functions listed above, you can specify another
%   kernel function by using @ (such as @normpdf) or quotes (such as 'normpdf').
%   The function must take a single argument that is an array of distances
%   between data values and places where the density is evaluated, and
%   return an array of the same size containing corresponding values of
%   the kernel function.  When the 'function' parameter value is 'pdf',
%   this kernel function should return density values, otherwise it should
%   return cumulative probability values.  Specifying a custom kernel when the
%   'function' parameter value is 'icdf' is an error.
%
%   If the 'support' parameter is 'positive', KSDENSITY transforms X using
%   a log function, estimates the density of the transformed values, and
%   transforms back to the original scale.  If 'support' is a vector [L U],
%   KSDENSITY uses the transformation log((X-L)/(U-X)).  The 'width' parameter
%   and U outputs are on the scale of the transformed values.
%
%   Example:
%      x = [randn(30,1); 5+randn(30,1)];
%      [f,xi] = ksdensity(x);
%      plot(xi,f);
%   This example generates a mixture of two normal distributions, and
%   plots the estimated density.
%
%   See also HIST, @.

% If there is any censoring, we would like to estimate the density up
% to the last non-censored observation.  Say this is XMAX.  Without
% censoring, the density estimate near XMAX would consist of contributions
% from kernels centered above and below XMAX.  We can't compute the
% contributions above XMAX, though, because we have no data.  Using only
% the kernels centered below XMAX makes the density estimate biased.
%
% In an attempt to reduce bias, we will compute the contributions
% from kernels centered below XMAX, and fold their values around XMAX.
% The result should be good if the density is nearly flat in this area.
% If the density is increasing then the estimate will still be biased
% downward, and if the density is decreasing it will still be biased
% upward, but the bias will be reduced.

% Reference:
%   A.W. Bowman and A. Azzalini (1997), "Applied Smoothing
%      Techniques for Data Analysis," Oxford University Press.

%   Copyright 1993-2006 The MathWorks, Inc.
%   $Revision: 1.9.6.11 $  $Date: 2006/12/15 19:30:08 $

if (nargin > 0) && isscalar(yData) && ishandle(yData) ...
                && isequal(get(yData,'type'),'axes')
    axarg = {yData};
    if nargin>1
       yData = varargin{1};
       varargin(1) = [];
    else
       yData = [];  % error to be dealt with below
    end
else
    axarg = {};
end

% Get y vector and its dimensions
if ~isvector(yData) || isempty(yData)
   error('stats:ksdensity:VectorRequired','X must be a non-empty vector.');
end
yData = yData(:);
yData(isnan(yData)) = [];
n = length(yData);
ymin = min(yData);
ymax = max(yData);

% Maybe xi was specified, or maybe not
xispecified = false;
if ~isempty(varargin)
   if ~ischar(varargin{1})
      xi = varargin{1};
      varargin(1) = [];
      xispecified = true;
   end
end

% Process additional name/value pair arguments
okargs={'width'     'npoints'   'kernel'   'support' ...
        'weights'   'censoring' 'cutoff'   'function'   };
defaults = {[]      []          'normal'   'unbounded' ...
            1/n     false(n,1)  []         'pdf'};
[eid,emsg,u,m,kernelname,support,weight,cens,cutoff,ftype] = ...
               statgetargs(okargs, defaults, varargin{:});
if ~isempty(eid)
   error(sprintf('stats:ksdensity:%s',eid),emsg);
end

if isnumeric(support)
   if numel(support)~=2
      error('stats:ksdensity:BadSupport',...
            'Value of ''support'' parameter must have two elements.');
   end
   if support(1)>=ymin || support(2)<=ymax
      error('stats:ksdensity:BadSupport',...
            'Data values must be between lower and upper ''support'' values.');
   end
   L = support(1);
   U = support(2);
elseif ischar(support) && length(support)>0
   okvals = {'unbounded' 'positive'};
   rownum = strmatch(support,okvals);
   if isempty(rownum)
      error('stats:ksdensity:BadSupport',...
            'Invalid value of ''support'' parameter.')
   end
   support = okvals{rownum};
   if isequal(support,'unbounded')
      L = -Inf;
      U = Inf;
   else
      L = 0;
      U = Inf;
   end
   if isequal(support,'positive') && ymin<=0
      error('stats:ksdensity:BadSupport',...
            'Cannot set support to ''positive'' with non-positive data.')
   end
else
   error('stats:ksdensity:BadSupport',...
         'Invalid value of ''support'' parameter.')
end
if isempty(weight)
   weight = ones(1,n);
elseif numel(weight)==1
   weight = repmat(weight,1,n);
elseif numel(weight)~=n || numel(weight)>length(weight)
   error('stats:ksdensity:InputSizeMismatch',...
         'Value of ''weight'' must be a vector of the same length as X.');
else
   weight = weight(:)';
end
weight = weight / sum(weight);
if isempty(cens)
   cens = false(1,n);
elseif ~all(ismember(cens(:),0:1))
   error('stats:ksdensity:BadCensoring',...
         'Value of ''censoring'' must be a logical vector.');
elseif numel(cens)~=n || numel(cens)>length(cens)
   error('stats:ksdensity:InputSizeMismatch',...
         'Value of ''censoring'' must be a vector of the same length as X.');
end

% Kernel can be the name of a function local to here, or an external function
% Kernel numbers are used below:
%              1        2               3               4        5
kernelnames = {'normal' 'epanechinikov' 'epanechnikov'  'box'    'triangle'};
kernelhndls = {@normal  @epanechnikov   @epanechnikov   @box     @triangle};
cdfhndls    = {@cdf_nl  @cdf_ep         @cdf_ep         @cdf_bx  @cdf_tr};
kernelcuts  = [4        sqrt(5)         sqrt(5)         sqrt(3)  sqrt(6)];

% Check function type
okvals = {'pdf' 'cdf' 'survivor' 'cumhazard' 'icdf'};
if ischar(ftype)
   ftype = lower(ftype);
end
rownum = strmatch(ftype,okvals);
if isempty(rownum)
   error('stats:ksdensity:BadFunction',...
         'Invalid value of ''function'' parameter.')
elseif length(rownum)>1
   error('stats:ksdensity:BadFunction',...
         'Ambiguous value of ''function'' parameter.')
end
ftype = okvals{rownum};

% Set a flag indicating we are to compute the cdf; later on
% we may transform to another function that is a transformation
% of the cdf
iscdf = isequal(ftype,'cdf') | isequal(ftype,'survivor') ...
                             | isequal(ftype,'cumhazard');
kernel = kernelname;
if isempty(kernelname)
   if iscdf
      kernel = cdfhndls{1};
   else
      kernel = kernelhndls{1};
   end
   kernelname = kernelnames{1};
   kernelcutoff = kernelcuts(1);
elseif ischar(kernelname)
   % If this is an abbreviation of our own methods, expand the name now.
   % If the string matches the start of both variants of the Epanechnikov
   % spelling, that is not an error so pretend it matches just one.
   knum = strmatch(lower(kernelname), kernelnames);
   if all(ismember(2:3,knum))   % kernel number used here
      knum(knum==3) = [];
   end
   if (length(knum) == 1)
      if iscdf
         kernel = cdfhndls{knum};
      else
         kernel = kernelhndls{knum};
      end
      kernelcutoff = kernelcuts(knum);
   else % custom kernel specified by name
      if isequal(ftype,'icdf')
         error('stats:ksdensity:IcdfNotAllowed',...
               'Cannot compute inverse cdf for a custom kernel.');
      end
      kernelcutoff = Inf;
   end
elseif (isa(kernelname,'function_handle') || isa(kernelname,'inline')) % custom kernel
   if isequal(ftype,'icdf')
      error('stats:ksdensity:IcdfNotAllowed',...
            'Cannot compute inverse cdf for a custom kernel.');
   else
      kernelcutoff = Inf;
   end
else
   error('stats:ksdensity:BadKernel',...
         'Smoothing kernel must be a function name or a function handle.');
end

if isempty(cutoff)
   cutoff = kernelcutoff;
end

% Inverse cdf is special, so deal with it here
if isequal(ftype,'icdf')
   if xispecified
      p = xi;
   else
      p = (1:99)/100;
   end
   
   if isequal(kernelname,'normal')
      % Truncation for the normal kernel creates small jumps in the CDF.
      % That's not a problem for the CDF, but it causes convergence problems
      % for ICDF calculation, so use a cutoff large enough to make the jumps
      % smaller than the convergence criterion.
      cutoff = max(cutoff,6);
   else
      % Other kernels have a fixed finite width.  Ignore any requested
      % truncation for these kernels; it would cause convergence problems if
      % smaller than the kernel width, and would have no effect if larger.
      cutoff = kernelcutoff;
   end

   % To get starting x values for the ICDF evaluated at p, first create a
   % grid xi of values spanning the data on which to evaluate the CDF
   sy = sort(yData);
   xi = linspace(sy(1), sy(end), 100);
   
   % Estimate the CDF on the grid
   [Fi,xi,u] = ksdensity(yData,xi, 'censoring',cens, 'kernel',kernelname, ...
                       'support',support, 'weights',weight, 'width',u, ...
                       'function','cdf');
   
   % If there are any gaps in the data wide enough to create regions of
   % exactly zero density, include points at the edges of those regions
   % in the grid, to make sure a linear interpolation smooth of the gridded
   % CDF captures them as constant
   halfwidth = cutoff*u;
   gap = find(diff(sy) > 2*halfwidth);
   if ~isempty(gap)
      sy = sy(:)';
      xi = sort([xi, sy(gap)+halfwidth, sy(gap+1)-halfwidth]);
      [Fi,xi,u] = ksdensity(yData,xi, 'censoring',cens, 'kernel',kernelname, ...
                          'support',support, 'weights',weight, 'width',u, ...
                          'function','cdf');
   end
   
   % Find any regions where the CDF is constant, these will cause problems
   % inverse interpolation for x at p
   t = (diff(Fi) == 0);
   if any(t)
       % Remove interior points in constant regions, they're unnecessary
       s = ([false t] & [t false]);
       Fi(s) = [];
       xi(s) = [];
       % To make Fi monotonic, nudge up the CDF value at the end of each
       % constant region by the smallest amount possible.
       t = 1 + find(diff(Fi) == 0);
       Fi(t) = Fi(t) + eps(Fi(t));
       % If the CDF at the point following is that same value, just remove
       % the nudgee.
       if (t(end) == length(Fi)), t(end) = []; end
       s = t(Fi(t) >= Fi(t+1));
       Fi(s) = [];
       xi(s) = [];
   end

   % Get starting values for ICDF(p) by inverse linear interpolation of
   % the gridded CDF, plus some clean-up
   x1 = interp1(Fi,xi,p);               % interpolate for p in a good range
   x1(isnan(x1) & p<min(Fi)) = min(xi); % use lowest x if p>0 too low
   x1(isnan(x1) & p>max(Fi)) = max(xi); % use highest x if p<1 too high
   x1(p<=0) = L;                        % use lower bound if p<=0
   x1(p>=1) = U;                        % and upper bound if p>=1
   
   % Now refine the ICDF using Newton's method for cases with 0<p<1
   notdone = find(p>0 & p<1);
   maxiter = 100;
   for iter = 1:maxiter
      x0 = x1(notdone);

      % Compute cdf and derivative (pdf) at this value
      F0 = ksdensity(yData,x0, 'censoring',cens, 'kernel',kernelname, ...
                       'support',support, 'weights',weight, 'width',u, ...
                       'function','cdf', 'cutoff',cutoff);
      dF0 = ksdensity(yData,x0, 'censoring',cens, 'kernel',kernelname, ...
                       'support',support, 'weights',weight, 'width',u, ...
                       'function','pdf', 'cutoff',cutoff);

      % Perform a Newton's step
      dp = p(notdone) - F0;
      dx = dp ./ dF0;
      x1(notdone) = x0 + dx;

      % Continue if the x and function (probability) change are large
      notdone = notdone(abs(dx) > 1e-6*abs(x0) & abs(dp) > 1e-8);
      if isempty(notdone), break; end
   end
   if iter>maxiter
      warning('stats:ksdensity:NoConvergence',...
              'Inverse CDF calculation did not converge for p = %g.', ...
              p(notdone(1)));
   end

   % Plot or return these values
   if nargout==0
       plot(axarg{:},x1,p);
   else
       fout0 = x1;
       xout = p;
   end
   return
end

% Compute transformed values of data
if isequal(support,'unbounded')
   ty = yData;
elseif isequal(support,'positive')
   ty = log(yData);
else % finite support [L, U]
   ty = log(yData-L) - log(U-yData);    % same as log((y-L)./(U-y))
end

% Deal with censoring
iscensored = any(cens);
if iscensored
   % Compute empirical cdf and create an equivalent weighted sample
   [F,XF] = ecdf(ty, 'censoring',cens, 'frequency',weight);
   weight = diff(F(:)');
   ty = XF(2:end);
   n = length(ty);
   N = sum(~cens);
   issubdist = (F(end)<1);  % sub-distribution, integrates to less than 1
   ymax = max(yData(~cens));
else
   N = n;
   issubdist = false;
end

% Get bandwidth if not already specified
if (isempty(u)),
   if ~iscensored
      % Get a robust estimate of sigma
      med = median(ty);
      sig = median(abs(ty-med)) / 0.6745;
   else
      % Estimate sigma using quantiles from the empirical cdf
      Xquant = interp1(F,XF,[.25 .5 .75]);
      if ~any(isnan(Xquant))
         % Use interquartile range to estimate sigma
         sig = (Xquant(3) - Xquant(1)) / (2*0.6745);
      elseif ~isnan(Xquant(2))
         % Use lower half only, if upper half is not available
         sig = (Xquant(2) - Xquant(1)) / 0.6745;
      else
         % Can't easily estimate sigma, just get some indication of spread
         sig = ty(end) - ty(1);
      end
   end
   if sig<=0, sig = max(ty)-min(ty); end
   if sig>0
      % Default window parameter is optimal for normal distribution
      u = sig * (4/(3*N))^(1/5);
   else
      u = 1;
   end
end

% Get XI values at which to evaluate the density
foldwidth = min(cutoff,3);
if ~xispecified
   % Compute untransformed values of lower and upper evaluation points
   ximin = min(ty) - foldwidth*u;
   if issubdist
      ximax = max(ty);
   else
      ximax = max(ty) + foldwidth*u;
   end

   if isequal(support,'positive')
      ximin = exp(ximin);
      ximax = exp(ximax);
   elseif ~isequal(support,'unbounded')
      ximin = (U*exp(ximin)+L) / (exp(ximin)+1);
      ximax = (U*exp(ximax)+L) / (exp(ximax)+1);
   end

   if isempty(m)
      m=100;
   end

   xi = linspace(ximin, ximax, m);

elseif ~isvector(xi)
   error('stats:ksdensity:VectorRequired','XI must be a vector');
end

% Compute transformed values of evaluation points that are in bounds
xisize = size(xi);
fout = zeros(xisize);
if iscdf && isfinite(U)
   fout(xi>=U) = 1;
end
xout = xi;
xi = xi(:);
if isequal(support,'unbounded')
   inbounds = true(size(xi));
   txi = xi;
   foldpoint = ymax;
elseif isequal(support,'positive')
   inbounds = (xi>0);
   xi = xi(inbounds);
   txi = log(xi);
   foldpoint = log(ymax);
else % finite support [L, U]
   inbounds = (xi>L) & (xi<U);
   xi = xi(inbounds);
   txi = log(xi-L) - log(U-xi);
   foldpoint = log(ymax-L) - log(U-ymax);
end
m = length(txi);


% If the density is censored at the end, add new points so that
% we can fold them back across the censoring point as a crude
% adjustment for bias
if issubdist && ~iscdf
   needfold = (txi >= foldpoint - foldwidth*u);
   nkeep = length(txi);
   nfold = sum(needfold);
   txifold = (2*foldpoint) - txi(needfold);
   txi(end+1:end+nfold) = txifold;
   m = length(txi);
else
   nkeep = length(txi);
   nfold = 0;
end

% Now compute density estimate at selected points
blocksize = 3e4;
if n*m<=blocksize && ~iscdf
   % For small problems, compute kernel density estimate in one operation
   z = (repmat(txi',n,1)-repmat(ty,1,m))/u;
   f = weight * feval(kernel, z);
else
   % For large problems, try more selective looping

   % First sort y and carry along weights
   [ty,idx] = sort(ty);
   weight = weight(idx);

   % Loop over evaluation points
   f = zeros(1,m);

   if isinf(cutoff)
      for k=1:m
         % Sum contributions from all
         z = (txi(k)-ty)/u;
         f(k) = weight * feval(kernel,z);
      end
   else
      % Sort evaluation points and remember their indices
      [stxi,idx] = sort(txi);

      jstart = 1;       % lowest nearby point
      jend = 1;         % highest nearby point
      for k=1:m
         % Find nearby data points for current evaluation point
         halfwidth = cutoff*u;
         lo = stxi(k) - halfwidth;
         while(ty(jstart)<lo && jstart<n)
            jstart = jstart+1;
         end
         hi = stxi(k) + halfwidth;
         jend = max(jend,jstart);
         while(ty(jend)<=hi && jend<n)
            jend = jend+1;
         end
         nearby = jstart:jend;

         % Sum contributions from these points
         z = (stxi(k)-ty(nearby))/u;
         fk = weight(nearby) * feval(kernel,z);
         if iscdf
            fk = fk + sum(weight(1:jstart-1));
         end
         f(k) = fk;
      end

      % Restore original x order
      f(idx) = f;
   end

end

% If we added extra points for folding, fold them now
if nfold>0
   % Fold back over the censoring point to give a crisp upper limit
   ffold = f(nkeep+1:end);
   f = f(1:nkeep);
   f(needfold) = f(needfold) + ffold;
   txi = txi(1:nkeep);
   f(txi>foldpoint) = 0;

   % Include a vertical line at the end
   if ~xispecified
      xi(end+1) = xi(end);
      f(end+1) = 0;
      inbounds(end+1) = true;
   end
end


if iscdf
   % Guard against roundoff.  Lower boundary of 0 should be no problem.
   f = min(1,f);
else
   % Apply reverse transformation and create return value of proper size
   f = f(:) ./ u;
   if isequal(support,'positive')
      f = f ./ xi;
   elseif isnumeric(support)
      f = f * (U-L) ./ ((xi-L) .* (U-xi));
   end
end
fout(inbounds) = f;
xout(inbounds) = xi;


% If another function based on the cdf, compute it now
if isequal(ftype,'survivor')
   fout = 1-fout;
elseif isequal(ftype,'cumhazard')
   fout = 1-fout;
   t = (fout>0);
   fout(~t) = NaN;
   fout(t) = -log(fout(t));
end

% Plot the results if they are not requested as return values
if nargout==0
    plot(axarg{:},xout,fout)
else
    fout0 = fout;
end

% -----------------------------
% The following are functions that define smoothing kernels k(z).
% Each function takes a single input Z and returns the value of
% the smoothing kernel.  These sample kernels are designed to
% produce outputs that are somewhat comparable (differences due
% to shape rather than scale), so they are all probability
% density functions with unit variance.
%
% The density estimate has the form
%    f(x;k,h) = mean over i=1:n of k((x-y(i))/h) / h

function f = normal(z)
%NORMAL Normal density kernel.
%f = normpdf(z);
f = exp(-0.5 * z .^2) ./ sqrt(2*pi);

function f = epanechnikov(z)
%EPANECHNIKOV Epanechnikov's asymptotically optimal kernel.
a = sqrt(5);
z = max(-a, min(z,a));
f = max(0,.75 * (1 - .2*z.^2) / a);

function f = box(z)
%BOX    Box-shaped kernel
a = sqrt(3);
f = (abs(z)<=a) ./ (2 * a);

function f = triangle(z)
%TRIANGLE Triangular kernel.
a = sqrt(6);
z = abs(z);
f = (z<=a) .* (1 - z/a) / a;

% -----------------------------
% The following are functions that define cdfs for smoothing kernels.

function f = cdf_nl(z)
%CDF_NL Normal kernel, cdf version
f = normcdf(z);

function f = cdf_ep(z)
%CDF_EP Epanechnikov's asymptotically optimal kernel, cdf version
a = sqrt(5);
z = max(-a, min(z,a));
f = ((z+a) - (z.^3+a.^3)/15) * 3 / (4*a);

function f = cdf_bx(z)
%CDF_BX Box-shaped kernel, cdf version
a = sqrt(3);
f = max(0, min(1,(z+a)/(2*a)));

function f = cdf_tr(z)
%CDF_TR Triangular kernel, cdf version
a = sqrt(6);
denom = 12;  % 2*a^2
f = zeros(size(z));                     % -Inf < z < -a
t = (z>-a & z<0);
f(t) = (a + z(t)).^2 / denom;           % -a < z < 0
t = (z>=0 & z<a);
f(t) = .5 + z(t).*(2*a-z(t)) / denom;   % 0 < z < a
t = (z>a);
f(t) = 1;                               % a < z < Inf
