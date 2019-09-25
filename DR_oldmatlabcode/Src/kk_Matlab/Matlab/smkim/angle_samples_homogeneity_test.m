function stats = angle_samples_homogeneity_test(obs,labels,varargin)
%ANGLE_SAMPLES_HOMOGENEITY_TEST Test whether samples of angles are drawn from same circular distribution.
%
%   ANGLE_SAMPLES_HOMOGENEITY_TEST performs a nonparametric test on the
%   hypothesis that the angular observations in vector OBS, grouped into samples
%   by labels LABELS, are drawn from a common (possibly multimodal)
%   distribution.
%
%   OBS is an N-element vector. LABELS is an N-element cell array or a character
%   array whose N rows are string labels for grouping observations. LABELS must
%   contain at least two unique string labels.
%
%   STATS = ANGLE_SAMPLES_HOMOGENEITY_TEST(OBS,LABELS) returns a struct with the
%   following fields:
%   groups: cell array of group labels
%        n: vector of group sizes
%       df: degrees of freedom for the asymptotic chi-squared distribution under
%           the null hypothesis, or NaN if the asymptotic approximation was not
%           used
%       Wr: value of the test statistic
%     nsim: number of random permutations drawn for computing p-value. If
%           METHOD is not specified, the permutation test is used by default 
%           when any of the samples has less than 10 observations. If NSIM is
%           not provided, the default number of permutations is 1e5.
%        p: p-value
%   method: 'approximation' or 'permutation'
%
%   ANGLE_SAMPLES_HOMOGENEITY_TEST(OBS,LABELS,'approximation') uses the
%   asymptotic chi-squared distribution to compute the p-value. This argument
%   overrides the default behavior for small sample sizes.
%
%   ANGLE_SAMPLES_HOMOGENEITY_TEST(OBS,LABELS,'permutation') uses a permutation
%   test to compute the p-value, using randomized labels to obtain a null
%   distribution for the test statistic.
%
%   ANGLE_SAMPLES_HOMOGENEITY_TEST(OBS,LABELS,'permutation',NSIM) does the
%   permutation test with NSIM random permutations.
%
%   The test is described in section 8.3 (pages 130-152) Brunner E., Domhof S.,
%   Langer F. (2002) Nonparametric Analysis of Longitudinal Data in Factorial
%   Experiments. Wiley-Interscience
%
%Depends on:
%   TIEDRANK (MATLAB Statistics toolbox)
%   CHI2CDF (MATLAB Statistics toolbox)
%   GRP2IDX (MATLAB Statistics toolbox)
%
%Written by smk, 18 October 2008
%

if (exist('tiedrank') ~= 2)
  error(['ANGLE_SAMPLES_HOMOGENEITY_TEST depends on m-file TIEDRANK ' ...
      '(in MATLAB Statistics Toolbox)']);
end
if (exist('tcdf') ~= 2)
  error(['ANGLE_SAMPLES_HOMOGENEITY_TEST depends on m-file TCDF ' ...
      '(in MATLAB Statistics Toolbox)']);
end
if (exist('grp2idx') ~= 2)
  error(['ANGLE_SAMPLES_HOMOGENEITY_TEST depends on m-file GRP2IDX ' ...
      '(in MATLAB Statistics Toolbox)']);
end

% default number of permutation draws if nsim is not specified
DEFAULT_NSIM = 1e5;

% if labels are not provided, create a labels array whose elements are all the
% empty string ''
if (nargin == 1)
  labels = cellstr(char(ones([size(obs,1),1])));
end
if ~isfloat(obs) || ~isreal(obs) || ~isvector(obs) || ~all(isfinite(obs))
  error('observations must be real finite floating-point vector');
end
if any((obs < -pi) | (obs > +pi))
  warning('some observations are outside the interval [-pi,+pi]');
end
if ischar(labels) && ~isempty(labels)
  % convert character array to cell array
  labels = cellstr(labels);
elseif iscellstr(labels)
  % pass
else
  error('labels argument must be a character array or cell array of string');
end

% total number of observations
N = numel(obs);
% sample sizes
[groupidx, stats.labels] = grp2idx(labels);
% number of groups
r = numel(stats.labels);
if (r < 2)
  error('LABELS must contain at least two distinct values');
end
% count number of subjects in each group
ni = diff(find([1; diff([sort(groupidx); r+1])]));
if (N ~= sum(ni)) || (numel(ni) ~= r)
  error('bug: subjects are not counted correctly');
end
stats.n = ni;
lookup_idx = cell(size(ni));
for i = 1:r
  lookup_idx{i} = find(groupidx == i);
end

% pick default method according to sample size
if any(ni <= 10)
  stats.method = 'permutation';
  stats.nsim = DEFAULT_NSIM; 
else
  stats.method = 'approximation';
  stats.nsim = NaN;
end

if (length(varargin) == 0)
  % pass
elseif (length(varargin) == 1)
  if ischar(varargin{1}) && ...
      any(strcmp(varargin{1},{'permutation','approximation'}))
    stats.method = varargin{1};
    stats.nsim = DEFAULT_NSIM;
  else
    error('invalid METHOD argument specified');
  end
elseif (length(varargin) == 2)
  if ~strcmp(stats.method,'permutation')
    error('NSIM is a valid argument only with the permutation method');
  end
  if isnumeric(varargin{2}) && isscalar(varargin{2}) && ...
      isreal(varargin{2}) && (varargin{2} > 0) && ...
      (round(varargin{2}) == varargin{2})
    stats.nsim = varargin{2};
  else
    error('NSIM must be a positive integer');
  end
else
  error('too many arguments');
end

% compute circular ranks over all data
gamma = 2*pi*tiedrank(obs)/N;
z = nan(size(ni));
for i = 1:r
  z(i) = sum(complex(cos(gamma(lookup_idx{i})),sin(gamma(lookup_idx{i}))));
end
% compute test statistic
stats.Wr = 2*sum(z .* conj(z) ./ ni)/r;

switch stats.method
case 'approximation'
  stats.nsim = NaN;
  stats.df = 2*(r-1);
  stats.p = 1 - chi2cdf(stats.Wr,stats.df);
case 'permutation'
  Wr_randomized = nan([stats.nsim 1]);
  for i = 1:stats.nsim
    gamma_randomized = gamma(randperm(N));
    for j = 1:r
      z(j) = sum(complex(cos(gamma_randomized(lookup_idx{j})), ...
          sin(gamma_randomized(lookup_idx{j}))));
    end
    Wr_randomized(i) = 2*sum(z .* conj(z) ./ ni)/r;
  end
  stats.df = NaN;
  m = nnz(stats.Wr > Wr_randomized);
  stats.p = (stats.nsim - m)/stats.nsim;
otherwise
  error('unrecognized method');
end


