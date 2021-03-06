function out = ITPCpermtest(ixpc, phdata, as, ianimal, introde, outInd, condAInds, condBInds, n_permutes)


%% clear previous results
ixpc.phasemean_h0{ianimal}{introde}{outInd} = [];ixpc.phasestd_h0{ianimal}{introde}{outInd} = [];ixpc.phasezmap{ianimal}{introde}{outInd} = [];
ixpc.powermean_h0{ianimal}{introde}{outInd} = [];ixpc.powerstd_h0{ianimal}{introde}{outInd} = [];ixpc.powerzmap{ianimal}{introde}{outInd} = [];
ixpc.MCmax_clust_size{ianimal}{introde}{outInd} = [];ixpc.MCphase_minmax{ianimal}{introde}{outInd} = []; ixpc.MCpower_minmax{ianimal}{introde}{outInd} = [];

%% create the shuffled distribution of indices for each pixel
disp(sprintf('========== running perm test ============= nt%d perm x %d',introde, n_permutes))
permdataInds = [condAInds; condBInds];
[~, permsetfull] = sort(rand(n_permutes,length(permdataInds)),2); % make permuted indices mat without loop

%% For every vec array of specifying permuted indices into the time-freq maps for this condition, compute the
% cross-condition ITPC and mean power

testo = arrayfun(@(x) bsxfun(@minus, abs(mean(exp(1i*phdata(:,permdataInds(permsetfull(x,1:length(condAInds))),:)),2)), ...
    abs(mean(exp(1i*phdata(:,permdataInds(permsetfull(x,length(condAInds)+1:end)),:)),2))), [1:n_permutes], 'un', 0);
phasepermCat = cat(2,testo{:});

presto = arrayfun(@(x) bsxfun(@minus, mean(abs(as(:,permdataInds(permsetfull(x,1:length(condAInds))),:)).^2,2), ...
    mean(abs(as(:,permdataInds(permsetfull(x,length(condAInds)+1:end)),:)).^2,2)), [1:n_permutes], 'un', 0);
powerpermCat = cat(2,presto{:});

%% compute pixel-wise mean and standard deviation maps for perm null
ixpc.phasemean_h0{ianimal}{introde}{outInd}(:,1,:) = mean(phasepermCat,2);
ixpc.phasestd_h0{ianimal}{introde}{outInd}(:,1,:)  = std(phasepermCat,[],2);
% pixel-wise Z-score the data against the perm null distribution for ITPC
ixpc.phasezmask{ianimal}{introde}{outInd}(:,1,:) = bsxfun(@rdivide, bsxfun(@minus, ixpc.phaseoutput{ianimal}{introde}(:,outInd,:), ixpc.phasemean_h0{ianimal}{introde}{outInd}), ...
    ixpc.phasestd_h0{ianimal}{introde}{outInd});
% % pixel-wise Z-score the data against the perm null distribution for power
ixpc.powermean_h0{ianimal}{introde}{outInd}(:,1,:) = mean(powerpermCat,2);
ixpc.powerstd_h0{ianimal}{introde}{outInd}(:,1,:)  = std(powerpermCat,[],2);
ixpc.powerzmask{ianimal}{introde}{outInd}(:,1,:) = bsxfun(@rdivide, bsxfun(@minus, ixpc.poweroutput{ianimal}{introde}(:,outInd,:), ixpc.powermean_h0{ianimal}{introde}{outInd}),...
    ixpc.powerstd_h0{ianimal}{introde}{outInd});

%% pixel-based multiple comparisons null distribution of extreme values
ixpc.MCphase_minmax{ianimal}{introde}{outInd} = sort([min(min(phasepermCat,[],3),[],1) max(max(phasepermCat,[],3),[],1)]);
ixpc.MCpower_minmax{ianimal}{introde}{outInd} = sort([min(min(powerpermCat,[],3),[],1) max(max(powerpermCat,[],3),[],1)]);

out = ixpc;
end