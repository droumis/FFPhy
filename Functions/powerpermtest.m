
function [powerzmask, MC_power_minmax] = powerpermtest(pwr, aIdx, bIdx, varargin)
% compute the power z mask for a two-condition power diff

% as : analytic signal. samples x events x freqs
% poweroutput : diff of mean power. samples x 1 x freqs
% aIdx : condition A indices into event array
% bIdx : condition B indices into event array
% options : n_permutes
% Author: Demetris Roumis 2019

% once this is working, call it from within the getPower nt loop so i don't
% have to load the individual AS's twice

% need to compute the significance mask for the individual settypes.. and
% then i can get to the diffs that I was doing previously.

n_permutes = 1000;
if ~isempty(varargin)
    assign(varargin{:});
end

for ian = 1:length(Fp.animals)
    animal = Fp.animals{ian};
    ntrodes = lfpstack(ian).ntrodes;
%     pwr = cell(1,length(ntrodes));
    for nti = 1:length(ntrodes)
%         dbpwr = cell(1,length(ripstate(ian).statesetsfields));
        nt = ntrodes(nti);
        as = loadAS(animal, nt, Fp.waveSet, 'AS');
        
        permdataInds = aIdx;
        
        
        %% create the shuffled distribution of indices for each pixel
        permdataInds = [aIdx; bIdx];
        [~, permsetfull] = sort(rand(n_permutes,length(permdataInds)),2);
        
        %% compute the cross-condition mean power diff for every permutation of indices
        powerperms = arrayfun(@(x) bsxfun(@minus, ...
            mean(abs( as(:,permdataInds(permsetfull(x,1:length(aIdx))),:)) .^2,2), ...
            mean(abs( as(:,permdataInds(permsetfull(x,length(aIdx)+1:end)),:)) .^2,2)), ...
            1:n_permutes, 'un', 0);
        
        powerpermCat = cat(2,powerperms{:});
        
        %% compute pixel-wise mean and standard deviation maps for perm null
        % pixel-wise Z-score the data against the perm null distribution for power
        powermean_h0(:,1,:) = mean(powerpermCat,2);
        powerstd_h0(:,1,:)  = std(powerpermCat,[],2);
        powerzmask(:,1,:) = bsxfun(@rdivide, bsxfun(@minus, pwr, powermean_h0),...
            powerstd_h0);
        
        %% pixel-based multiple comparisons null distribution of extreme values
        MC_power_minmax = sort(abs([min(min(powerpermCat,[],3),[],1) max(max(powerpermCat,[],3),[],1)]));
    end
end

end