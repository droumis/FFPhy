

function out = getPower(expvarCat, rawpwr, Fp, varargin)
%
% Fp = filter params (see load_filter_params)
% wp = wave params (see getWaveParams)
pconf = paramconfig;
me = animaldef('Demetris');

lfptype = 'eeg';
expvars = {'onlywdays','rewarded', 'unrewarded', 'inbound' , 'outbound', 'proximalWell', ...
    'distalWell'};
saveout = 1;
run_perm = 1;
outdir = 'expvarCatMeanPwr';
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
    meandbpower = cell(1,length(expvars));
    
    % for each state/condition, compute dbpower, run timeshift permtest vs baseline
    for stset = 1:length(expvars)
        fprintf('ripstate %s \n', expvars{stset});
        stidx = find(strcmp(expvars{stset}, expvarCat(ian).expvars));
            sripidx = find(expvarCat(ian).dm(:,stidx));
            meandbpower{stset} = computePower(rawpwr(ian).pwr(:,:,sripidx,:), wp, ...
                'dsamp',wp.dsamp, 'run_permutation_test', run_perm);
    end
    
    out(ian).animal = animal;
    out(ian).dm = expvarCat;
    out(ian).expvars = expvars;
    out(ian).wp = wp;
    out(ian).Fp = Fp;
    out(ian).lfptype = lfptype;
    out(ian).meandbpower = meandbpower;
    
    if saveout
        save_data(out, [pconf.andef{2},outdir,'/'], [outdir, '_', Fp.epochEnvironment])
    end
end
end

function [pout] = computePower(pwr,wp,varargin)
fprintf('computing power\n');
tic
% dsamp = 10;
run_permutation_test = 1;
if ~isempty(varargin)
    assign(varargin{:});
end
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
fprintf('baseline %.02f:%.02f sec \n', wp.basewin(1), wp.basewin(2));
% fprintf('prewin %.02f:%.02f sec \n', wp.prewin(1), wp.prewin(2));
% fprintf('postwin %.02f:%.02f sec \n', wp.postwin(1), wp.postwin(2));

% pre and post mean per ntrode per event per freq
% pout.premean = squeeze(nanmean(pwr(:,preind(1):preind(2),:,:),2));
% pout.postmean = squeeze(nanmean(pwr(:,postind(1):postind(2),:,:),2));

pout.bl_pwr_mean = mean(mean(pwr(:,baseind(1):baseind(2),:,:),2),3);
pout.pwrmean = mean(pwr,3); % mean across events
pout.pwr_mean_db = squeeze(10*log10(bsxfun(@rdivide, pout.pwrmean, pout.bl_pwr_mean)));
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
    pout.permt = run_pwr_perm_test(pwr, pout.pwr_mean_db, pout.bl_pwr_mean, wp);
else
    fprintf('not running permtest\n');
end
fprintf('power took %.02f seconds \n', toc)

end

function permt = run_pwr_perm_test(pwr, dbmean, bl_dbmean, wp)
fprintf('running %d permutes \n', wp.n_permutes);
tic
cutpnts = randsample(1:size(dbmean,2),wp.n_permutes,true);
powerperms = cell(1,length(wp.n_permutes));
parfor p = 1:length(cutpnts)
    powerperms{p} = bsxfun(@rdivide, ...
        mean(pwr(:,[cutpnts(p):end 1:cutpnts(p)-1],:,:),3),bl_dbmean);
end
permutes_vals = 10*log10(cat(2,powerperms{:}));
permutes_max = max(max(permutes_vals, [],1),[],3);
permutes_min = min(min(permutes_vals, [],1),[],3);

zmap = squeeze((dbmean-mean(permutes_vals,2)) ./ std(permutes_vals,[],2));
threshmean = squeeze(dbmean);
threshmean(abs(zmap)<norminv(1-wp.voxel_pval))=0;
permt.zmap = zmap;
permt.threshmean = threshmean;
permt.permutes_max = permutes_max;
permt.permutes_min = permutes_min;
fprintf('perm test took %.03f sec \n', toc)
end