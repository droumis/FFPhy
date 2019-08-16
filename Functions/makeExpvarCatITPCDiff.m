
function out = makeExpvarCatITPCDiff(expvarCat, phase, Fp, varargin)
%
% Fp = filter params (see load_filter_params)
% wp = wave params (see getWaveParams)
pconf = paramconfig;
me = animaldef('Demetris');
noiseEvents = struct;
lfptype = 'eeg';
expvars = {{'rewarded', 'unrewarded'},{'outbound', 'inbound'},{'distalWell', 'proximalWell'},...
    {'rewarded_outbound', 'rewarded_inbound'}};
saveout = 1;
run_perm = 1;
outdir = 'expvarCatITPCDiff';
if ~isempty(varargin)
    assign(varargin{:});
end

for ian = 1:length(phase)
    wp = getWaveParams(Fp.waveSet);
    animal = phase(ian).animal;
    fprintf('animal %s \n', animal);
    andef = animaldef(animal);
    tetinfo = loaddatastruct(andef{2}, animal, 'tetinfo');
    den = cellfetch(tetinfo, '');
    ntrodes = unique(den.index(:,3));
    ITPCDiff = cell(1,length(expvars));
    evcatanidx = find(strcmp({expvarCat.animal}, animal));
    userips = ones(length(rawpwr(ian).day),1);
    if ~isempty(noiseEvents)
        % exclude invalid rips
        noiseanidx = find(strcmp({noiseEvents.animal}, animal));
        if ~isempty(noiseEvents(noiseanidx).events)
            invalidrips = ismember([phase(ian).day phase(ian).epoch phase(ian).ripStartTime], ...
                noiseEvents(noiseanidx).events, 'rows');
            userips(invalidrips) = 0;
        end
    end
    
    for s = 1:length(expvars)
        stidx = find(strcmp(expvars{s}{1}, expvarCat(evcatanidx).expvars));
        Aidx = expvarCat(evcatanidx).dm(:,stidx);
        Aidx = find(all([Aidx, userips], 2));
        stidx = find(strcmp(expvars{s}{2}, expvarCat(evcatanidx).expvars));
        Bidx = expvarCat(evcatanidx).dm(:,stidx);
        Bidx = find(all([Bidx, userips], 2));
        ITPCDiff{s} = computeITPCDiff(phase(ian).ph, Aidx, Bidx, wp, 'dsamp', ...
            wp.dsamp, 'run_permutation_test', run_perm);
    end
    out(ian).animal = animal;
    out(ian).dm = expvarCat(evcatanidx);
    out(ian).expvars = expvars;
    out(ian).wp = wp;
    out(ian).Fp = Fp;
    out(ian).lfptype = lfptype;
    out(ian).ITPCDiff = ITPCDiff;
    if saveout
        save_data(out, [pconf.andef{2},outdir,'/'], [outdir, '_', Fp.epochEnvironment])
    end
end
end
function [pout] = computeITPCDiff(ph,Aidx,Bidx,wp,varargin)
fprintf('computing ITPC\n');
tic
% dsamp = 10;
run_permutation_test = 1;
if ~isempty(varargin)
    assign(varargin{:});
end
phA = ph(:,:,Aidx,:);
phB = ph(:,:,Bidx,:);
% DIMS ntrodes x samples x events x freqs
srate = wp.srate/wp.dsamp;
timeWin = wp.win(1):1/srate:wp.win(2);
baseind(1,1) = dsearchn(timeWin',wp.basewin(1));
baseind(1,2) = dsearchn(timeWin',wp.basewin(2));

pout.bl_ITPCA = mean(abs(mean(exp(1i*phA(:,baseind(1):baseind(2),:,:)),3)),2);
pout.ITPCA = abs(mean(exp(1i*phA),3)); % itpc across events
pout.ITPC_dbA = squeeze(10*log10(bsxfun(@rdivide, pout.ITPCA, pout.bl_ITPCA)));

pout.bl_ITPCB = mean(abs(mean(exp(1i*phB(:,baseind(1):baseind(2),:,:)),3)),2);
pout.ITPCB = abs(mean(exp(1i*phB),3)); % itpc across events
pout.ITPC_dbB = squeeze(10*log10(bsxfun(@rdivide, pout.ITPCB, pout.bl_ITPCB)));

pout.ITPC_diff = pout.ITPC_dbA - pout.ITPC_dbB;

pout.dsamp = dsamp;
pout.dims = {'ntrode', 'sample', 'frequency'};

pout.permt = struct;
if run_permutation_test
    pout.permt = run_perm_test(ph, Aidx, Bidx, pout.ITPC_diff, wp);
else
    fprintf('not running permtest\n');
end
fprintf('took %.02f seconds \n', toc)
end


function permt = run_perm_test(ph,Aidx,Bidx,phdiff,wp)
fprintf('running %d permutes \n', wp.n_permutes);
tic
didx = [Aidx; Bidx];
[~, permsidx] = sort(rand(wp.n_permutes,length(didx)),2); % make permuted indices mat without loop
phperms = cell(1,length(wp.n_permutes));
for p = 1:2%wp.n_permutes
    pA = squeeze(abs(mean(exp(1i*ph(:,:,didx(permsidx(p,1:length(Aidx))),:)),3)));
    pB = squeeze(abs(mean(exp(1i*ph(:,:,didx(permsidx(p,length(Aidx)+1:end)),:)),3))); 
    phperms{p} = pA - pB;
end
phpermCat = cat(4,phperms{:});
% % pixel-wise Z-score the data against the perm null distribution for power
phmean_h0 = mean(phpermCat,4);
phstd_h0  = std(phpermCat,[],4);
zmap = squeeze(bsxfun(@rdivide, bsxfun(@minus, phdiff, phmean_h0),...
    phstd_h0));
threshmean = squeeze(phdiff);
threshmean(abs(zmap)<norminv(1-wp.voxel_pval))=0;
permt.zmap = zmap;
permt.threshmean = threshmean;
fprintf('perm test took %.03f sec \n', toc)
end