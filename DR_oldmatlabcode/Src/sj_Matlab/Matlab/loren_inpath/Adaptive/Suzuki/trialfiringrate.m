function frates = trialfiringrate (c, cluster, binsize, r, group)
% trialfiringrate.m
% This function returns trial-by-trial firing rates.
%
% Syntax1: FRATES = FIRINGRATE (COBJ, CLUSTER, BINSIZE, ROBJ, [GROUP])
% Syntax2: FRATES = FIRINGRATE (SPMAT, BINSIZE)
%
% where
% FRATES:    M-by-N matrix of spikes/second rates, M trials by N bins;
% COBJ:      the @cortex object;
% SPMAT:     a spike matrix returned by @cortex\spikematrix.m
% CLUSTER:   the spike number, as encoded in the CORTEX data file;
% BINSIZE:   binsize in ms;
% ROBJ:      one rule from a @rules object 
%            or the entire object, if GROUP is specified;
% GROUP:     the group to be analyzed (optional).
%
% Last modified: 28 May 98

% retrieve the spike array

switch nargin
case 5
   spmat = spikematrix (c, cluster, r, group);
case 4
   if sum(r.size)>1
      error ('''TRIALFIRINGRATE'' ERROR: multi-rule object passed using single-rule syntax');
   end
   spmat = spikematrix (c, cluster, r);
case 2
   spmat = c;
   binsize = cluster;
   if ndims (spmat) ~= 2
      error ('''TRIALFIRINGRATE'' ERROR: if you pass 2 args, the first should be a spike matrix');
   end
otherwise
   error ('''TRIALFIRINGRATE'' ERROR: wrong number of arguments');
end

if prod (size(binsize))>1
   error ('''TRIALFIRINGRATE'' ERROR: binsize should be a scalar');
end

% make sure that the length of the array is a multiple of binsize
[instances, winsize] = size(spmat);
extra = rem (winsize, binsize);
if extra
   warning ('''TRIALFIRINGRATE'' WARNING: time window is not a multiple of binsize;');
   disp (['rightmost ' num2str(extra) 'ms truncated from spike matrix.']);
   winsize = winsize - extra;
   spmat = spmat (:,1:winsize);
end

% compute the firing rates
numofbins = winsize / binsize;
spmat = reshape(spmat, instances, binsize, numofbins);
spptr =  double(sum (spmat, 2)) / binsize * 1000;
frates  = squeeze (spptr);