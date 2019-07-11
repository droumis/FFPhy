

function out = getITPC(ripstate, Fp, varargin)
% Fp = filter params (see load_filter_params)

uselfptype = 'eeg';
ripstatetypes = {'all'};
saveresult = 1;
run_permutation_test = 0;
me = animaldef('Demetris');
if ~isempty(varargin)
    assign(varargin{:});
end

out = struct;
for ian = 1:length(Fp.animals)
    wp = getWaveParams(Fp.waveSet);
    animal = Fp.animals{ian};
    fprintf('animal %s \n', animal);
    andef = animaldef(animal);
    tetinfo = loaddatastruct(andef{2}, animal, 'tetinfo');
    den = cellfetch(tetinfo, '');
    ntrodes = unique(den.index(:,3));
    
    if isempty(ripstatetypes)
        ripstatetypes = ripstate(ian).statesetsfields;
    end
    
    tic
    ph = load_data(sprintf('%s/analyticSignal/', me{2}), ...
        sprintf('AS_waveSet-%s_%s_phase', wp.waveSet, uselfptype), animal);
    fprintf('%d seconds to load AS phase\n', toc);
    fprintf('lfptype %s \n', uselfptype);
    ITPC = cell(1,length(ripstate(ian).statesetsfields));
    
    % for each state/condition, compute dbpower, run timeshift permtest vs baseline
    for stset = 1:length(ripstatetypes)
        fprintf('ripstate %s \n', ripstatetypes{stset});
        stidx = find(strcmp(ripstatetypes{stset}, ripstate(ian).statesetsfields));
        sripidx = find(ripstate(ian).statesets(:,stidx));
        if strcmp(ripstatetypes{stset}, 'all') % avoid slice copy if 'all'
            ITPC{stset} = computeITPC(ph.ph, wp, ...
                'dsamp',wp.dsamp, 'run_permutation_test', run_permutation_test);
        else
            ITPC{stset} = computeITPC(ph.ph(:,:,sripidx,:), wp, ...
                'dsamp',wp.dsamp, 'run_permutation_test', run_permutation_test);
        end
    end
    
    out(ian).animal = animal;
    out(ian).ripstate = ripstate;
    out(ian).ripstatetypes = ripstatetypes;
    out(ian).wp = wp;
    out(ian).Fp = Fp;
    out(ian).lfptype = uselfptype;
    out(ian).ITPC = ITPC;
    
    if saveresult
        save_ITPC(out(ian), Fp, uselfptype)
    end
end
end

function [pout] = computeITPC(ph,wp,varargin)
fprintf('computing ITPC\n');
tic
dsamp = 10;
run_permutation_test = 1;
if ~isempty(varargin)
    assign(varargin{:});
end
% DIMS ntrodes x samples x events x freqs
srate = wp.srate/wp.dsamp;
timeWin = wp.win(1):1/srate:wp.win(2);
baseind(1,1) = dsearchn(timeWin',wp.basewin(1));
baseind(1,2) = dsearchn(timeWin',wp.basewin(2));
preind(1,1) = dsearchn(timeWin',wp.prewin(1));
preind(1,2) = dsearchn(timeWin',wp.prewin(2));
postind(1,1) = dsearchn(timeWin',wp.postwin(1));
postind(1,2) = dsearchn(timeWin',wp.postwin(2));

fprintf('event time 0\n');
fprintf('win %.02f:%.02f sec\n',wp.win(1), wp.win(2));
fprintf('baseline %.02f:%.02f sec \n', wp.basewin(1), wp.basewin(2));
fprintf('prewin %.02f:%.02f sec \n', wp.prewin(1), wp.prewin(2));
fprintf('postwin %.02f:%.02f sec \n', wp.postwin(1), wp.postwin(2));

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
