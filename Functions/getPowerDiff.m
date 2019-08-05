
function out = getPowerDiff(expvarCat, rawpwr, Fp, varargin)
%
% Fp = filter params (see load_filter_params)
% wp = wave params (see getWaveParams)
pconf = paramconfig;
me = animaldef('Demetris');
invalid_rips = struct;
lfptype = 'eeg';
expvars = {{'rewarded', 'unrewarded'},{'outbound', 'inbound'},{'distalWell', 'proximalWell'},...
    {'rewarded_outbound', 'rewarded_inbound'}};
saveout = 1;
run_perm = 1;
outdir = 'expvarCatMeanPwrDiff';
if ~isempty(varargin)
    assign(varargin{:});
end

for ian = 1:length(Fp.animals)
    wp = getWaveParams(Fp.waveSet);
    animal = Fp.animals{ian};
    fprintf('animal %s \n', animal);
    andef = animaldef(animal);
    tetinfo = loaddatastruct(andef{2}, animal, 'tetinfo');
    den = cellfetch(tetinfo, '');
    ntrodes = unique(den.index(:,3));
    fprintf('lfptype %s \n', lfptype);
    meandbpowerDiff = cell(1,length(expvars));
    pwranidx = find(strcmp({rawpwr.animal}, animal));
    if ~isempty(invalid_rips)
        % exclude invalid rips
        noiseanidx = find(strcmp({invalid_rips.animal}, animal));
        invalidrips = invalid_rips(noiseanidx ).ripnums;
        userips = ones(length(rawpwr(pwranidx).day),1);
        userips(invalidrips) = 0;
    end
    % for each state/condition, compute dbpower, run timeshift permtest vs baseline
    for s = 1:length(expvars)
        stidx = find(strcmp(expvars{s}{1}, expvarCat(ian).expvars));
        Aidx = expvarCat(ian).dm(:,stidx);
        Aidx = find(all([Aidx, userips], 2));
        stidx = find(strcmp(expvars{s}{2}, expvarCat(ian).expvars));
        Bidx = expvarCat(ian).dm(:,stidx);
        Bidx = find(all([Bidx, userips], 2));
        meandbpowerDiff{s} = computePowerDiff(rawpwr(pwranidx).pwr, Aidx, Bidx, wp, 'dsamp',wp.dsamp, ...
            'run_permutation_test', run_perm);
    end
    
    out(ian).animal = animal;
    out(ian).dm = expvarCat;
    out(ian).expvars = expvars;
    out(ian).wp = wp;
    out(ian).Fp = Fp;
    out(ian).lfptype = lfptype;
    out(ian).meandbpowerDiff = meandbpowerDiff;
    
    if saveout
        save_data(out, [pconf.andef{2},outdir,'/'], [outdir, '_', Fp.epochEnvironment])
    end
end
end

function [pout] = computePowerDiff(pwr,Aidx,Bidx, wp,varargin)
fprintf('computing power diff\n');
tic
% dsamp = 10;
run_permutation_test = 1;
if ~isempty(varargin)
    assign(varargin{:});
end
pwrA = pwr(:,:,Aidx,:);
pwrB = pwr(:,:,Bidx,:);
% PWR DIMS ntrodes x samples x events x freqs
% pwr = abs(as(1:dsamp:end,:,:,:)).^2; %default dsamp is 2: 1.5kHz to .75kHz sampling
srate = wp.srate/wp.dsamp;
timeWin = wp.win(1):1/srate:wp.win(2);
baseind(1,1) = dsearchn(timeWin',wp.basewin(1));
baseind(1,2) = dsearchn(timeWin',wp.basewin(2));
% preind(1,1) = dsearchn(timeWin',wp.prewin(1));
% preind(1,2) = dsearchn(timeWin',wp.prewin(2));
% postind(1,1) = dsearchn(timeWin',wp.postwin(1));
% postind(1,2) = dsearchn(timeWin',wp.postwin(2));

fprintf('event time 0\n');
fprintf('win %.02f:%.02f sec\n',wp.win(1), wp.win(2));
% fprintf('baseline %.02f:%.02f sec \n', wp.basewin(1), wp.basewin(2));
% fprintf('prewin %.02f:%.02f sec \n', wp.prewin(1), wp.prewin(2));
% fprintf('postwin %.02f:%.02f sec \n', wp.postwin(1), wp.postwin(2));

% pre and post mean per ntrode per event per freq
% pout.premean = squeeze(nanmean(pwr(:,preind(1):preind(2),:,:),2));
% pout.postmean = squeeze(nanmean(pwr(:,postind(1):postind(2),:,:),2));

pout.bl_pwr_meanA = mean(mean(pwrA(:,baseind(1):baseind(2),:,:),2),3);
pwrmeanA = mean(pwrA,3); % mean across events
pout.pwr_mean_db_A = squeeze(10*log10(bsxfun(@rdivide, pwrmeanA, pout.bl_pwr_meanA)));

pout.bl_pwr_meanB = mean(mean(pwrB(:,baseind(1):baseind(2),:,:),2),3);
pwrmeanB = mean(pwrB,3); % mean across events
pout.pwr_mean_db_B = squeeze(10*log10(bsxfun(@rdivide, pwrmeanB, pout.bl_pwr_meanB)));

pout.pwr_meandb_diff = pout.pwr_mean_db_A - pout.pwr_mean_db_B;

pout.dsamp = dsamp;
pout.dims_pwr_mean_db = {'ntrode', 'sample', 'frequency'};

% baseline_std_power = std(mean(abs( as(baseind(1):baseind(2),:,:)) .^2,2),[],1);
% pctpower = 100 * bsxfun(@rdivide, bsxfun(@minus, power, basepower), basepower);
% meanZpower = (power-baseline_power) ./ baseline_std_power;
% medianpower = median(abs(as).^2,2);
% medianbasepower = median(median(abs( as(wp.baseind(1):wp.baseind(2),:,:)) .^2,2),1);
% mediandbpower = 10*log10(bsxfun(@rdivide, medianpower, medianbasepower));

pout.permt = struct;
if run_permutation_test
    pout.permt = run_pwr_perm_test(pwr,Aidx,Bidx,pout.pwr_meandb_diff,wp);
else
    fprintf('not running permtest\n');
end
fprintf('power took %.02f seconds \n', toc)
end

function permt = run_pwr_perm_test(pwr,Aidx,Bidx,dbmeandiff,wp)
fprintf('running %d permutes \n', wp.n_permutes);
tic
didx = [Aidx; Bidx];
[~, prms] = sort(rand(wp.n_permutes,length(didx)),2); % make permuted indices mat without loop
powerperms = cell(1,length(wp.n_permutes));
for p = 1:length(wp.n_permutes)
    powerperms{p} = bsxfun(@minus, mean(pwr(:,didx(prms(p,1:length(Aidx))),:),2), ...
    mean(pwr(:,didx(perms(p,length(Aidx)+1:end)),:),2));
end
powerpermCat = cat(2,powerperms{:});
% % pixel-wise Z-score the data against the perm null distribution for power
powermean_h0(:,1,:) = mean(powerpermCat,2);
powerstd_h0(:,1,:)  = std(powerpermCat,[],2);
zmap = squeeze(bsxfun(@rdivide, bsxfun(@minus, dbmeandiff, powermean_h0),...
    powerstd_h0));
threshmean = squeeze(dbmeandiff);
threshmean(abs(zmap)<norminv(1-wp.voxel_pval))=0;
permt.zmap = zmap;
permt.threshmean = threshmean;
fprintf('perm test took %.03f sec \n', toc)
end