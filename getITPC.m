

function out = getITPC(expvarCat, , Fp, varargin)
% Fp = filter params (see load_filter_params)

pconf = paramconfig;
me = animaldef('Demetris');
noiseEvents = struct;
lfptype = 'eeg';
expvars = {'onlywdays'}; %,'rewarded', 'unrewarded', 'inbound' , 'outbound', 'proximalWell', ...
%     'distalWell'};
saveout = 1;
run_perm = 1;
outdir = 'expvarCatITPC';
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
    evCatanidx = find(strcmp({expvarCat.animal}, animal));
    ITPC = cell(1,length(expvarCat(evCatanidx).expvars));
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
    
    for stset = 1:length(expvars)
        fprintf('ripstate %s \n', expvars{stset});
        stidx = find(strcmp(expvars{stset}, expvarCat(evCatanidx).expvars));
        sripidx = expvarCat(evCatanidx).dm(:,stidx);
        sripidx(~userips) = 0;
        ITPC{stset} = computeITPC(phase(ian).ph(:,:,find(sripidx),:), wp, ...
            'dsamp',wp.dsamp, 'run_permutation_test', run_perm);
    end
    out(ian).animal = animal;
    out(ian).dm = expvarCat(evCatanidx);
    out(ian).expvars = expvars;
    out(ian).wp = wp;
    out(ian).Fp = Fp;
    out(ian).lfptype = lfptype;
    out(ian).ITPC = ITPC;
    if saveout
        save_data(out, [pconf.andef{2},outdir,'/'], [outdir, '_', Fp.epochEnvironment])
    end
end
end

function [pout] = computeITPC(ph,wp,varargin)
fprintf('computing ITPC\n');
tic
% dsamp = 10;
run_permutation_test = 1;
if ~isempty(varargin)
    assign(varargin{:});
end
% DIMS ntrodes x samples x events x freqs
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
% pout.premean = squeeze(nanmean(ph(:,preind(1):preind(2),:,:),2));
% pout.postmean = squeeze(nanmean(ph(:,postind(1):postind(2),:,:),2));

pout.bl_ITPC = mean(abs(mean(exp(1i*ph(:,baseind(1):baseind(2),:,:)),3)),2);
pout.ITPC = abs(mean(exp(1i*ph),3)); % itpc across events
pout.ITPC_db = squeeze(10*log10(bsxfun(@rdivide, pout.ITPC, pout.bl_ITPC)));
pout.dsamp = dsamp;
pout.dims = {'ntrode', 'sample', 'frequency'};

pout.permt = struct;
if run_permutation_test
    pout.permt = run_perm_test(ph, pout.ITPC, pout.bl_ITPC, wp);
else
    fprintf('not running permtest\n');
end
fprintf('took %.02f seconds \n', toc)

end

function permt = run_perm_test(ph, itpc, bl_itpc, wp)
fprintf('running perm test : %d iterations \n', wp.n_permutes);
tic
cutpnts = randsample(1:size(itpc,2),wp.n_permutes,true);
perms = cell(1,length(wp.n_permutes));
parfor p = 1:length(cutpnts)
    perms{p} = bsxfun(@rdivide, ...
        mean(ph(:,[cutpnts(p):end 1:cutpnts(p)-1],:,:),3),bl_itpc);
end
permutes_vals = cat(2,perms{:});
permutes_max = max(max(permutes_vals, [],1),[],3);
permutes_min = min(min(permutes_vals, [],1),[],3);

zmap = squeeze((itpc-mean(permutes_vals,2)) ./ std(permutes_vals,[],2));
threshmean = squeeze(itpc);
threshmean(abs(zmap)<norminv(1-wp.voxel_pval))=0;
permt.zmap = zmap;
permt.threshmean = threshmean;
permt.permutes_max = permutes_max;
permt.permutes_min = permutes_min;
fprintf('took %.03f sec \n', toc)

end

function save_ITPC(out, Fp, lfptype)
me = paramconfig;
savedir = sprintf('%s/itpc/', me.andef{2});
savestr = sprintf('/itpc_waveSet-%s_%s_%s',Fp.waveSet,lfptype,Fp.epochEnvironment);
save_data(out, savedir, savestr)
end
