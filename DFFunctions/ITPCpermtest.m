function [powerzmask, phasezmask, MCpower_minmax, MCphase_minmax] = ITPCpermtest(poweroutput, phaseoutput, phdata, as, condAInds, condBInds, n_permutes)

%% create the shuffled distribution of indices for each pixel
permdataInds = [condAInds; condBInds];
[~, permsetfull] = sort(rand(n_permutes,length(permdataInds)),2); % make permuted indices mat without loop

%% For every vec array of specifying permuted indices into the time-freq maps for this condition, compute the
% cross-condition ITPC and mean power for every permutation of indices
phaseperms = arrayfun(@(x) bsxfun(@minus, abs(mean(exp(1i*phdata(:,permdataInds(permsetfull(x,1:length(condAInds))),:)),2)), ...
    abs(mean(exp(1i*phdata(:,permdataInds(permsetfull(x,length(condAInds)+1:end)),:)),2))),[1:n_permutes], 'un', 0);
phasepermCat = cat(2,phaseperms{:});

powerperms = arrayfun(@(x) bsxfun(@minus, mean(abs(as(:,permdataInds(permsetfull(x,1:length(condAInds))),:)).^2,2), ...
    mean(abs(as(:,permdataInds(permsetfull(x,length(condAInds)+1:end)),:)).^2,2)), [1:n_permutes], 'un', 0);
powerpermCat = cat(2,powerperms{:});

%% compute pixel-wise mean and standard deviation maps for perm null
phasemean_h0(:,1,:) = mean(phasepermCat,2);
phasestd_h0(:,1,:)  = std(phasepermCat,[],2);
% pixel-wise Z-score the data against the perm null distribution for ITPC
phasezmask(:,1,:) = bsxfun(@rdivide, bsxfun(@minus, phaseoutput, phasemean_h0), ...
    phasestd_h0);
% % pixel-wise Z-score the data against the perm null distribution for power
powermean_h0(:,1,:) = mean(powerpermCat,2);
powerstd_h0(:,1,:)  = std(powerpermCat,[],2);
powerzmask(:,1,:) = bsxfun(@rdivide, bsxfun(@minus, poweroutput, powermean_h0),...
    powerstd_h0);

%% pixel-based multiple comparisons null distribution of extreme values
MCphase_minmax = sort(abs([min(min(phasepermCat,[],3),[],1) max(max(phasepermCat,[],3),[],1)]));
MCpower_minmax = sort(abs([min(min(powerpermCat,[],3),[],1) max(max(powerpermCat,[],3),[],1)]));

end