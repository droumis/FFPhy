function cluster_quality = measure_cluster_quality(features,labels,varargin)
%MEAUSRE_CLUSTER_QUALITY Compute cluster quality measures
%   CLUSTER_QUALITY = MEASURE_CLUSTER_QUALITY(FEATURES,LABELS,INGROUP_IDX)
%   computes cluster quality measures for assessing single-unit isolation.
%   Given N events that are measured in some D-dimensional (real-valued)
%   feature space, summarized in the NxD matrix DATA and labeled with the
%   D-element cell array of strings LABELS, we want to measure how well the
%   subset of events specified by INGROUP_IDX is separated from those evenets
%   that do not belong to the subset (e.g., the set complement).
%
%   MEASURE_CLUSTER_QUALITY(FEATURES,LABELS,INGROUP_IDX,OUTGROUP_IDX) measures
%   the separation of the cluster specified by INGROUP_IDX with respect to the
%   outgroup cluster specified by OUTGROUP_IDX. If OUTGROUP_IDX is omitted, then
%   it is assumed to equal setdiff(1:size(FEATURES,1),INGROUP_IDX).
%
%   This function returns a CLUSTER_QUALITY struct, with the following fields:
%     dimensions (copy of LABELS)
%     L_ratio 
%     isolation_distance
%
%   L_ratio and isolation_distance are described in Schmitzer-Torbert N.,
%   Jackson J., Henze D., Harris K., Redish A.D. (2005) Quantitative measures
%   of cluster quality for use in extracellular recordings. _Neuroscience_
%   131:1-11.
%
%   This function is adapted from code written by Ken Harris, who bears no
%   responsibility for any introduced errors by me (SMK).
%
%Depends on:
%   MAHAL (MATLAB Statistics Toolbox)
%   CHI2CDF (MATLAB Statistics Toolbox)
%
%Written by SMK, 2009 February 12.
%

if (exist('mahal') ~= 2)
  error(['MEASURE_CLUSTER_QUALITY depends on m-file MAHAL ' ...
      '(in MATLAB Statistics Toolbox)']);
end

if (exist('chi2cdf') ~= 2)
  error(['MEASURE_CLUSTER_QUALITY depends on m-file CHI2CDF ' ...
      '(in MATLAB Statistics Toolbox)']);
end

if isempty(features) || ~isnumeric(features) || (ndims(features) ~= 2) || ...
    ~isreal(features) || any(isinf(features(:))) || ~isfloat(features)
  error('FEATURES must be a two-dimensional finite real floating-point array');
end
if ~iscellstr(labels) || (numel(labels) ~= size(features,2)) || ...
    (numel(labels) ~= numel(unique(labels)))
  error(['LABELS must be a cell array of unique strings whose elements ' ...
      'correspond to the columns of FEATURES']);
end 

if length(varargin) > 2
  error('Too many arguments');
elseif length(varargin) == 2
  ingroup_idx = varargin{1};
  outgroup_idx = varargin{2};
elseif length(varargin) == 1
  ingroup_idx = varargin{1};
  outgroup_idx = setdiff(1:size(features,1),ingroup_idx);
else
  error('No cluster indices specified');
end

if isempty(ingroup_idx)
  warning('INGROUP_IDX is empty');
elseif isvector(ingroup_idx) && islogical(ingroup_idx) && ...
    (numel(ingroup_idx) == size(features,1))
  % If INGROUP_IDX is given as a logical bitmask vector, convert it to indices
  ingroup_idx = find(ingroup_idx);
elseif ~isvector(ingroup_idx) || ~isnumeric(ingroup_idx) || ...
      ~all(unique(ingroup_idx) == sort(ingroup_idx)) || ...
      ~all(round(ingroup_idx) == ingroup_idx) || ...
      ~all(ingroup_idx >= 1) || ~all(ingroup_idx <= size(features,1))
    error('INGROUP_IDX must be a vector of valid row indices into FEATURES');
end
if isempty(outgroup_idx)
  error('OUTGROUP_IDX is empty');
elseif isvector(outgroup_idx) && islogical(outgroup_idx) && ...
    (numel(outgroup_idx) == size(features,1))
  % If OUTGROUP_IDX is given as a logical bitmask vector, convert it to indices
  outgroup_idx = find(outgroup_idx);
elseif ~isvector(outgroup_idx) || ~isnumeric(outgroup_idx) || ...
      ~all(unique(outgroup_idx) == sort(outgroup_idx)) || ...
      ~all(round(outgroup_idx) == outgroup_idx) || ...
      ~all(outgroup_idx >= 1) || ~all(outgroup_idx <= size(features,1))
    error('OUTGROUP_IDX must be a vector of valid row indices into FEATURES');
end
%  Check for disjointness
if ~isempty(intersect(ingroup_idx,outgroup_idx))
  error('OUTGROUP_IDX and INGROUP_IDX overlap');
end
% Exclude events with NaN feature values
nan_idx = find(any(isnan(features),2));
if ~isempty(nan_idx)
  warning('%d events have NaN features and will be excluded',numel(nan_idx));
end
ingroup_idx = setdiff(ingroup_idx,nan_idx);
outgroup_idx = setdiff(outgroup_idx,nan_idx);

% Initialize output
cluster_quality = struct( ...
    'dimensions'        , {labels}, ...
    'L_ratio'           , {NaN}   , ...
    'isolation_distance', {NaN}   );
% Abort if cluster is empty
if isempty(ingroup_idx)
  warning('No valid events specified by INGROUP_IDX');
  return;
end
if isempty(outgroup_idx)
  warning('No valid events specified by OUTGROUP_IDX');
  return;
end

% compute 2-norm Mahalanobis distances from the outgroup to the ingroup (it
% matters which one we define as ingroup versus outgroup!)
outgroup2ingroup_D2 = mahal(features(outgroup_idx,:),features(ingroup_idx,:));

% L is computed from the chi-squared distribution of Mahalanobis distance, with
% degrees of freedom equal to the dimensionality of the feature space
L = sum(1 - chi2cdf(outgroup2ingroup_D2,size(features,2)));
cluster_quality.L_ratio = L/numel(ingroup_idx);

% isolation_distance can only be computed if the number of events in the
% cluster is less than the number of events in the outgroup and greater than or
% equal to the dimensionality of the feature space
if (numel(ingroup_idx) < size(features,2))
	warning(['isolation_distance could not be computed: ' ...
      'not enough events in cluster']);
elseif (numel(ingroup_idx) >= numel(outgroup_idx))
  warning(['isolation_distance could not be computed: ' ...
      'more than half of events belong to the cluster']);
else
  % isolation_distance is defined as the Mahalanobis distance from the cluster
  % which includes as many outgroup points as are inside the cluster
  sorted_D2 = sort(outgroup2ingroup_D2(:));
  cluster_quality.isolation_distance = sorted_D2(numel(ingroup_idx));
end

%{

% TODO: compute isolation_information
%
%   isolation_information is described in Fenton A.A., Kao H.-Y., Neymotin S.A.,
%   Olypher A., Vayntrub Y., Lytton W.W., Ludvig N. (2008) Unmasking the CA1
%   ensemble place code by exposures to small and large environments: more place
%   cells and multiple, irregularly arranged, and expanded place fields in the
%   larger space. _Journal of Neuroscience_ 28: 11250-11262.
%
% isolation information is the resistor-average symmetrized distance based on
% the Kullback-Leibler divergences between the probability distributions of the
% cluster versus non-cluster:
%
%cluster_quality.isolation_information = 1/ ...
%    (1/cluster2outgroup_KLdiv + 1/outgroup2ingroup_KLdiv);
%
% How do we estimate Kullback-Leibler divergence, given only the features points
% that are drawn from the distributions? One way is to histogram the features, and
% then use the K-T estimate [Krichevsky R.E., Trofimov V.K. (1981) The
% performance of universal encoding. _IEEE Transactions in Information Theory_
% 27:199-207]. Another way is to compute a smooth kernel density estimate, but
% principled selection of a multivariate kernel bandwidth is tricky. Yet
% another way to estimate KL divergence involves sampling the K nearest
% neighbors of each features point [Leonenko N., Pronzato L., Savani L. (2008)
% Estimation of entropies and divergences via nearest neighbors. _Tatra
% Mountains Mathematical Publications_ 39:265-273; Wang Q., Kulkarni S.R.,
% Verdu S. (2006) A nearest-neighbor approach to estimating divergence between
% continuous random vectors. _2006 IEEE Symposium on Information Theory_]. 
%

%}

