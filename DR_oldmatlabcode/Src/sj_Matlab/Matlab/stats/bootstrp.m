function [bootstat, bootsam] = bootstrp(nboot,bootfun,varargin)
%BOOTSTRP Bootstrap statistics.
%   BOOTSTAT = BOOTSTRP(NBOOT,BOOTFUN,D1,...) draws NBOOT bootstrap data
%   samples, computes statistics on each sample using the function BOOTFUN,
%   and returns the results in the matrix BOOTSTATS.  NBOOT must be a
%   positive integer.  BOOTFUN is a function handle specified with @.
%   Each row of BOOTSTAT contains the results of applying BOOTFUN to one
%   bootstrap sample.  If BOOTFUN returns a matrix or array, then this
%   output is converted to a row vector for storage in BOOTSTAT.
%
%   The third and later input arguments (D1,...) are data (scalars,
%   column vectors, or matrices) that are used to create inputs to BOOTFUN.
%   BOOTSTRP creates each bootstrap sample by sampling with replacement
%   from the rows of the non-scalar data arguments (these must have the
%   same number of rows).  Scalar data are passed to BOOTFUN unchanged.
%
%   [BOOTSTAT,BOOTSAM] = BOOTSTRP(...) returns BOOTSAM, a matrix of indices
%   into the rows of the extra arguments.  To get the output samples BOOTSAM
%   without applying a function, set BOOTFUN to empty ([]).
%
%   Examples:
%
%   Compute a sample of 100 bootstrapped means of random samples taken from
%   the vector Y, and plot an estimate of the density of these bootstrapped
%   means:
%      y = exprnd(5,100,1);
%      m = bootstrp(100, @mean, y);
%      [fi,xi] = ksdensity(m);
%      plot(xi,fi);
%
%   Compute a sample of 100 bootstrapped means and standard deviations of
%   random samples taken from the vector Y, and plot the bootstrap estimate
%   pairs:
%      y = exprnd(5,100,1);
%      stats = bootstrp(100, @(x) [mean(x) std(x)], y);
%      plot(stats(:,1),stats(:,2),'o')
%
%   Estimate the standard errors for a coefficient vector in a linear
%   regression by bootstrapping residuals:
%      load hald ingredients heat
%      x = [ones(size(heat)), ingredients];
%      y = heat;
%      b = regress(y,x);
%      yfit = x*b;
%      resid = y - yfit;
%      se = std(bootstrp(1000, @(bootr) regress(yfit+bootr,x), resid));
%
%   Bootstrap a correlation coefficient standard error:
%      load lawdata gpa lsat
%      se = std(bootstrp(1000,@corr,gpa,lsat));
%
%   See also RANDOM, RANDSAMPLE, HIST, KSDENSITY.

%   Reference:
%      Efron, Bradley, & Tibshirani, Robert, J.
%      "An Introduction to the Bootstrap", 
%      Chapman and Hall, New York. 1993.

%   Copyright 1993-2005 The MathWorks, Inc. 
%   $Revision: 2.11.2.7 $  $Date: 2005/12/12 23:33:21 $

% Initialize matrix to identify scalar arguments to bootfun.
la = length(varargin);
scalard = zeros(la,1);

% find out the size information in varargin.
n = 1;
for k = 1:la
   [row,col] = size(varargin{k});
   if max(row,col) == 1
      scalard(k) = 1;
   end
   if row == 1 && col ~= 1
      row = col;
      varargin{k} = varargin{k}(:);
   end
   n = max(n,row);
end

% Create index matrix of bootstrap samples if necessary
haveallsamples = (nargout>=2);
if haveallsamples
   bootsam = ceil(n*rand(n,nboot));
end

if isempty(bootfun)
   bootstat = zeros(nboot,0);
   return
end

% Get result of bootfun on actual data and find its size.
bootstat = feval(bootfun,varargin{:});

% Initialize an array to contain the results of all the bootstrap
% calculations, preserving the output type
bootstat(nboot,1:numel(bootstat)) = bootstat(:)';

% Do bootfun - nboot times.
if la==1 && ~haveallsamples && ~any(scalard)
   % For special case of one non-scalar argument and one output, try to be fast
   X1 = varargin{1};
   for bootiter = 1:nboot
      onesample = ceil(n*rand(n,1));
      tmp = feval(bootfun,X1(onesample,:));
      bootstat(bootiter,:) = (tmp(:))';
   end
elseif la==2 && ~haveallsamples && ~any(scalard)
   % For two non-scalar arguments and one output, try to be fast
   X1 = varargin{1};
   X2 = varargin{2};
   for bootiter = 1:nboot
      onesample = ceil(n*rand(n,1));
      tmp = feval(bootfun,X1(onesample,:),X2(onesample,:));
      bootstat(bootiter,:) = (tmp(:))';
   end
else
   % General case
   db = cell(la,1);
   for bootiter = 1:nboot
      if haveallsamples
         onesample = bootsam(:,bootiter);
      else
         onesample = ceil(n*rand(n,1));
      end
      for k = 1:la
         if scalard(k) == 0
            db{k} = varargin{k}(onesample,:);
         else
            db{k} = varargin{k};
         end
      end
      tmp = feval(bootfun,db{:});
      bootstat(bootiter,:) = (tmp(:))';
   end
end