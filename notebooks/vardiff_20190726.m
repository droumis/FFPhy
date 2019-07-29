

%{ 

now i need to make diffcors for power.
so start with outbound - inbound, rewarded-unrewarded

get all rips in the stack that are outbound, mean
get all rips in the stack that are inbound, mean
get diff 
permutes
i did diff with the original phase coherence stuff.
i also recently did the mean per outb, inb condition.. so i just need to
find that stuff and add diff and permuting to it

ripstates returns the outbound inbound rewarded unrewarded


--------------------------------
- need to make diffcorrs for phase
%}

load_stack = 0;

get_ripstate = 0;
load_ripstate = 0;
get_designmat = 0;
load_designmat = 0;
load_raw_pwr = 0;

get_tfdesign = 1;
resultname =  'trialType';
run_regression = 1;
load_regression = 0;

plot_regression = 1;
pausefigs = 0;
savefigs = 1;

run_TFVarcorr = 0;
plot_corrheatmap = 0;
run_multiple_regress = 0;

Fp.animals = {'D10'};
Fp.filtfunction = 'dfa_riptriglfp';
Fp.add_params = {'wtrackdays', 'excludeNoise','excludePriorFirstWell', '<4cm/s', ...
    'wavelets4-300Hz',  'excludeAfterLastWell'};

Fp = load_filter_params(Fp, 'add_params', Fp.add_params);
Fp.uselfptype = 'eeg';
Fp.useripstates = {'onlywdays','rewarded', 'unrewarded', 'inbound' , 'outbound', ...
    'proximalWell', 'distalWell'};
Fp.useDiffTrialTypes = {'diffRewardedUnrewarded', 'diffOutboundInbound'};
pconf = paramconfig;
me = animaldef('Demetris');

%% load lfpstack
if load_stack
    lfpstack = load_data(Fp.paths.resultsDirectory, 'riptriglfpstack_wtrack', Fp.animals);
end

%% load raw power
if load_raw_pwr
    wp = getWaveParams(Fp.waveSet);
    rawpwr = load_data(sprintf('%s/analyticSignal/', me{2}), ...
        sprintf('AS_waveSet-%s_%s_power', wp.waveSet, Fp.uselfptype), Fp.animals);
end

%% ripstates returns the outbound, inbound, rewarded, unrewarded
if get_ripstate
    ripstate = getStateFilters(lfpstack);
    save_data(ripstate, [pconf.andef{2}, 'ripstate/'], 'ripstate_wtrack');
end
if load_ripstate
    ripstate = load_data([pconf.andef{2}, 'ripstate/'], 'ripstate_wtrack', Fp.animals);
end

% 
% %% diffStates
% for ani = 1:length(Fp.animals)
%     animal = Fp.animals{ani};
%     diffTrialType.expvar = {};
%     diffTrialType(ani).dims = {'ripple', 'expvar'}; % dimensions of dm
%     for dtt = 1:length(Fp.useDiffTrialTypes)
%         ripstateidx = ismember(ripstate(ani).statesetsfields, Fp.useDiffTrialTypes{dtt});
%         dttA = 
%         diffTrialType(ani).dm(:,dtt) = % col design matrix
%         diffTrialType(ani).expvar{end+1} = Fp.useDiffTrialTypes{dtt};
%     end
% end

%% tfbox design matrix returns timefrex boxed pwr means, rip, gamma, theta, etc
if get_tfdesign
    tfboxDesign = makeTFBoxDesignMat(rawpwr);
end

%% Compute and Save MEAN Power, ITPC
% rely on the ripstates.. look up condA and condB to get their indices..
% then create a matrix of shuffled condA,condB indices x n_permutes
% then for each permute 
if calc_mean_PWR
    getDiffPower(ripstate, Fp, 'uselfptype', Fp.uselfptype, 'ripstatetypes', ...
        Fp.useDiffTrialTypes, 'savepower', 1);
end
if calc_ITPC
    getDiffITPC(ripstate, Fp, 'uselfptype', Fp.uselfptype, 'ripstatetypes', ...
        Fp.usediffripstates, 'saveresult', 1);
end


%% create shuffled indices 
condsIdx = [condAInds; condBInds];
[~, permsIdx] = sort(rand(wp.n_permutes,length(condsIdx)),2);

%% For every permute, cross-condition ITPC and power
% ITPCperms = arrayfun(@(x) bsxfun(@minus, abs(mean(exp(1i*phdata(:,permdataInds(permsetfull(x,1:length(condAInds))),:)),2)), ...
%     abs(mean(exp(1i*phdata(:,permdataInds(permsetfull(x,length(condAInds)+1:end)),:)),2))),[1:n_permutes], 'un', 0);
% ITPCpermCat = cat(2,ITPCperms{:});

nRips = size(rawpwr(ani).pwr,3);
nNTrodes = size(rawpwr(ani).pwr,1);
nfrex = size(rawpwr(ani).pwr,4);
nTimePoints = size(rawpwr(ani).pwr,2);
pwrPermMean = zeros(30, nrips, 
parfor p = 1:length(wp.n_permutes)
pwrPermMean(:,:,p,:) = squeeze(bsxfun(@minus, ...
    mean(rawpwr(ani).pwr(:,condsIdx(permsIdx(x,1:length(condsIdx))),:,:),3), ...
    mean(rawpwr(ani).pwr(:,condsIdx(permsIdx(x,length(condsIdx)+1:end)),:,:),3)));
end
powerpermCat = cat(2,pwrPermMean{:});

%% zmaps
ITPCmean_h0(:,1,:) = mean(ITPCpermCat,2);
ITPCstd_h0(:,1,:)  = std(ITPCpermCat,[],2);
% pixel-wise Z-score the data against the perm null distribution for ITPC
ITPCzmask(:,1,:) = bsxfun(@rdivide, bsxfun(@minus, ITPCoutput, ITPCmean_h0), ...
    ITPCstd_h0);
% % pixel-wise Z-score the data against the perm null distribution for power
powermean_h0(:,1,:) = mean(powerpermCat,2);
powerstd_h0(:,1,:)  = std(powerpermCat,[],2);
powerzmask(:,1,:) = bsxfun(@rdivide, bsxfun(@minus, poweroutput, powermean_h0),...
    powerstd_h0);

%% multiple comparisons correction
% pixel-based multiple comparisons null distribution of extreme values
MC_ITPC_minmax = sort(abs([min(min(ITPCpermCat,[],3),[],1) max(max(ITPCpermCat,[],3),[],1)]));
MC_power_minmax = sort(abs([min(min(powerpermCat,[],3),[],1) max(max(powerpermCat,[],3),[],1)]));







