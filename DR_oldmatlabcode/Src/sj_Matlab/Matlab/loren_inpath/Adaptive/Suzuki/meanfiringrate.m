function [frates, instances] = meanfiringrate (c, cluster, binsize, r, group)
% meanfiringrate.m
% This function returns average firing rates.
%
% Syntax 1: [FRATES, INSTANCES] = MEANFIRINGRATE (COBJ, CLUSTER, BINSIZE, ROBJ, [GROUP])
% Syntax 2: [FRATES, INSTANCES] = MEANFIRINGRATE (SPMAT, BINSIZE)
%
% where
% FRATES:    array of spikes/second rates, one per bin;
% INSTANCES: number of spike arrays retrieved from COBJ by the RULES criteria
%            i.e., the number of values included in the averages;
% COBJ:      the @cortex object;
% SPMAT:     a spike matrix returned by @cortex\spikematrix.m
% CLUSTER:   the spike number, as encoded in the CORTEX data file;
% BINSIZE:   binsize in ms;
% ROBJ:      one rule from a @rules object 
%            or the entire object, if GROUP is specified;
% GROUP:     the group to be analyzed (optional).
%
% Last modified: 2 Dec 98

% retrieve the spike array

switch nargin
case 5
   spmat = spikematrix (c, cluster, r, group);
case 4
   if sum(r.size)>1
      error ('multi-rule object passed using single-rule syntax');
   end
   spmat = spikematrix (c, cluster, r);
case 2
   spmat = c;
   binsize = cluster;
   if ndims (spmat) ~= 2
      error ('with 2 args, the first should be a spike matrix');
   end
otherwise
   error ('wrong number of arguments');
end

if prod (size(binsize))>1
   error ('binsize should be a scalar');
end
% make sure that the length of the array is a multiple of binsize
[instances, winsize] = size(spmat);
extra = rem (winsize, binsize);
if extra
   warning ('time window is not a multiple of binsize;');
   disp (['rightmost ' num2str(extra) 'ms truncated from spike matrix.']);
   winsize = winsize - extra;
   spmat = spmat (:,1:winsize);
end

% compute the firing rates
numofbins = winsize / binsize;


%length(sum(spmat))
%binsize
%numofbins

spmat = reshape(sum(spmat), binsize, numofbins);
spptr = double(sum(spmat, 1))/instances;
frates  = spptr / binsize * 1000;