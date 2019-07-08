

function getPower(ripstate, Fp, varargin)
% 
% Fp = filter params (see load_filter_params)
% wp = wave params (see getWaveParams)

lfptypes = {'eeg'};
ripstatetypes = {'all'};
savepower = 1;
run_permutation_test = 0;
dsamp = 2;
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
    
    % batch load and vectorize all AS data across ntrodes for this animal
    tic
    as = cell(1,length(lfptypes));
    for itype = 1:length(lfptypes)
        asnt = cell(1,length(ntrodes));
        parfor nti = 1:length(ntrodes)
            % this takes ~100Gb mem, 13-34 min per animal
            nt = ntrodes(nti);
            asnt{nti} = loadAS(Fp.animals{ian}, nt, Fp.waveSet, lfptypes{itype});
        end
        as{itype} = asnt;
        clear asnt
        a = [as{itype}{:}];
        clear as
        b = {a.as};
        clear a
        asdatacat{itype} = cat(4,b{:});
        clear b
    end
    fprintf('batch load AS took %d seconds \n', toc)
    
    if isempty(ripstatetypes)
        ripstatetypes = ripstate(ian).statesetsfields;
    end
    for itype = 1:length(lfptypes)
        pwrout = struct;
        lfptype = lfptypes{itype};
        fprintf('lfptype %s \n', lfptype);
        meandbpower = cell(1,length(ripstate(ian).statesetsfields));

        % for each state/condition, compute dbpower, run timeshift permtest vs baseline
        for stset = 1:length(ripstatetypes)
            fprintf('ripstate %s \n', ripstatetypes{stset});
            stidx = find(strcmp(ripstatetypes{stset}, ripstate(ian).statesetsfields));
            sripidx = find(ripstate(ian).statesets(:,stidx));
            meandbpower{stset} = computePower(asdatacat{itype}(:,sripidx,:,:), wp, ...
                'dsamp',dsamp, 'run_permutation_test', run_permutation_test);
        end

        pwrout.animal = animal;
        pwrout.ripstate = ripstate;
        pwrout.ripstatetypes = ripstatetypes;
        pwrout.wp = wp;
        pwrout.Fp = Fp;
        pwrout.lfptype = lfptype;
        pwrout.meandbpower = meandbpower;

        if savepower
            save_power(pwrout, animal, Fp.waveSet, lfptype)
        end
    end
end
end

function [pout] = computePower(as,wp,varargin)
fprintf('computing power\n');
tic
dsamp = 10;
run_permutation_test = 1;
if ~isempty(varargin)
    assign(varargin{:});
end
% AS dims = samples X events X freqs X ntrodes
pwr = abs(as(1:dsamp:end,:,:,:)).^2; %default dsamp is 2: 1.5kHz to .75kHz sampling
pout.bl_pwr_mean = mean(mean(abs( as(wp.baseind(1):wp.baseind(2),:,:,:)) .^2,2),1);
pout.pwr_mean_db = 10*log10(bsxfun(@rdivide, mean(pwr,2), pout.bl_pwr_mean));
pout.dsamp = dsamp;
% baseline_std_power = std(mean(abs( as(baseind(1):baseind(2),:,:)) .^2,2),[],1);
% pctpower = 100 * bsxfun(@rdivide, bsxfun(@minus, power, basepower), basepower);
% meanZpower = (power-baseline_power) ./ baseline_std_power;
% medianpower = median(abs(as).^2,2);
% medianbasepower = median(median(abs( as(wp.baseind(1):wp.baseind(2),:,:)) .^2,2),1);
% mediandbpower = 10*log10(bsxfun(@rdivide, medianpower, medianbasepower));

pout.permt = struct;
if run_permutation_test
    pout.permt = run_pwr_perm_test(pwr, pout.bl_pwr_mean, pout.pwr_mean_db);
else
    fprintf('not running permtest\n');
end
fprintf('power took %d seconds \n', toc)

end

function permt = run_pwr_perm_test(pwr, dbmean, bl_dbmean, wp)

    fprintf('running %d permutes \n', wp.n_permutes);
    tic
    cutpnts = randsample(1:size(dbmean,1),wp.n_permutes,true);
    powerperms = cell(1,length(wp.n_permutes));
    parfor p = 1:length(cutpnts) % parfor2 on virga, 4 on typhoon
        powerperms{p} = bsxfun(@rdivide, ...
            mean(pwr([cutpnts(p):end 1:cutpnts(p)-1],:,:,:),2),bl_dbmean);
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
    fprintf('perm test took %d seconds \n', toc)

end

function save_power(pwrout, animal, waveSet, lfptype)
me = animaldef('Demetris');
savedir = sprintf('%s/power', me{2});
if ~isdir(savedir)
    mkdir(savedir);
end
savestr = sprintf('%s/power_%s_%s_waveSet-%s_%s.mat', savedir, animal, ...
    'wtrack', waveSet, lfptype);
save(savestr, 'pwrout', '-v7.3');
fprintf('SAVED POWER SIGNAL RESULTS ++++++++++ %s \n',savestr)
end


% [powerout, basedmeanpower, percbasedpowerout,logbasedpowerout, logbasedpoweroutUnspec,percbasedpoweroutUnspec,basedmeanpowerUnspec]  =...
%     computePower(as,ian,int,iEpTetIndsType,baseind, indsSet)
% fprintf('$$$$$$$$$$$$$$$ POWER for animal %s nt %s  $$$$$$$$$$$$$$$$$$$ \n',ian, int);

% % condition-unspecific baseline mean
% %     basedmeanpower = mean(mean(abs(ias).^2,2),1); %mean for each freq across time and across all events for the baseline period
%     basedmeanpowerUnspec = mean(mean(abs(as(baseind(1):baseind(2),:,:)).^2,2),1);

% % for iSet = 1:size(iEpTetIndsType, 2)
%     % currently AS is samples x events x freq
%     % compute power as the magnitude squared from origin in complex space
%     powertmp = cellfun(@(x) mean(abs(as(:,x,:)).^2,2), iEpTetIndsType(:,iSet), 'un', 0);
%     powerout{iSet} = cat(2,powertmp{:}); %col cat the means for each event state type
%
%     %% baseline data, mean
%     % baselineASignal = as(baseind(1):baseind(2),iEpTetIndsType{1},:); %analytic signal from baseline period
%     % condition-specific baseline mean (mean across events and basline samples (times)
%     powerbaselinemeans = cellfun(@(x) mean(mean(abs(as(baseind(1):baseind(2),x,:)).^2,2),1),iEpTetIndsType(:,iSet), 'un', 0);
%     % powerbaselinemeans = cellfun(@(x) mean(mean(abs(baselineASignal(:,x,:)).^2,2),1),iEpTetIndsType, 'un', 0);
%     basepower{iSet} = cat(2,powerbaselinemeans{:});
%
%     %% Normalization
%     % % z-score normalization using only the baseline window to compute mean, std; condition-specific baseline mean
%     % powerbaselinestds = cellfun(@(x) std(reshape(abs(ias(:,x,:)).^2,[],1,numfrex),1),iEpTetIndsType, 'un', 0);
%     % basedstdpower = cat(2,powerbaselinestds{:});
%     % zscoredbasedpowerOut = bsxfun(@rdivide, bsxfun(@minus, powerOut,basedmeanpower), basedstdpower);
%     %
%     % % z-score normalization using the entire window to compute mean, std;  condition-specific baseline mean
%     % % vert cat all the traces into columns of frex, then std within freq of the baseline period
%     % powermeans = cellfun(@(x) mean(mean(abs(as{int}(:,x,:)).^2,2),1),iEpTetIndsType, 'un', 0);
%     % meanpower = cat(2,powermeans{:});
%     % powerstds = cellfun(@(x) std(reshape(abs(as{int}(:,x,:)).^2,[],1,numfrex),1),iEpTetIndsType, 'un', 0);
%     % stdpower = cat(2,powerstds{:});
%     % zscoredpowerOut = bsxfun(@rdivide, bsxfun(@minus, powerOut,meanpower), stdpower);
%
%     % percent change from baseline - normalization using condition-specific baseline mean
%     % (power - baseline mean) /  baseline mean
%     pctpower{iSet} = 100 * bsxfun(@rdivide, bsxfun(@minus, powerout{iSet}, basepower{iSet}), basepower{iSet});
%
%     % decibel - normalization using condition-specific baseline mean
%     % 10 * log10( power / baseline mean )
%     logpower{iSet} = 10*log10(bsxfun(@rdivide, powerout{iSet},basepower{iSet}));
%
% % %     % percent change from baseline - normalization using condition-UNspecific baseline mean
% % %     percbasedpoweroutUnspec{iSet} = 100 * bsxfun(@rdivide, bsxfun(@minus, powerout{iSet},basedmeanpowerUnspec), basedmeanpowerUnspec);
% % %
% % %     % decibel - normalization using condition-UNspecific baseline mean
% % %     logbasedpoweroutUnspec{iSet} = 10*log10(bsxfun(@rdivide, powerout{iSet},basedmeanpowerUnspec));
% % %
% %     %% compute differential tf maps
% % %     if strcmp(indsSet, 'performance') || strcmp(indsSet, 'performanceByDay')
% % %         try
% %             powerout{iSet}(:,6,:) = powerout{iSet}(:,4,:) - powerout{iSet}(:,5,:);
% %             powerout{iSet}(:,7,:) = powerout{iSet}(:,2,:) - powerout{iSet}(:,3,:);
% %             logbasedpowerout{iSet}(:,6,:) = logbasedpowerout{iSet}(:,4,:) - logbasedpowerout{iSet}(:,5,:);
% %             logbasedpowerout{iSet}(:,7,:) = logbasedpowerout{iSet}(:,2,:) - logbasedpowerout{iSet}(:,3,:);
% %             pctbasedpowerout{iSet}(:,6,:) = pctbasedpowerout{iSet}(:,4,:) - pctbasedpowerout{iSet}(:,5,:);
% %             pctbasedpowerout{iSet}(:,7,:) = pctbasedpowerout{iSet}(:,2,:) - pctbasedpowerout{iSet}(:,3,:);
% %
% %             % % i actually don't think you want to use the raw difference. instead use
% %             % the already normalized, then differenced
% %             % tmpinds = [iEpTetIndsType{4}; iEpTetIndsType{5}];
% %             % basedmeanpower(:,6,:) = mean(mean(abs(baselineASignal(:,tmpinds,:)).^2,2),1);
% %             % tmpinds = [iEpTetIndsType{2}; iEpTetIndsType{3}];
% %             % basedmeanpower(:,7,:) = mean(mean(abs(baselineASignal(:,tmpinds,:)).^2,2),1);
% %
% % %             logbasedpoweroutUnspec{iSet}(:,6,:) = logbasedpoweroutUnspec{iSet}(:,4,:) - logbasedpoweroutUnspec{iSet}(:,5,:);
% % %             logbasedpoweroutUnspec{iSet}(:,7,:) = logbasedpoweroutUnspec{iSet}(:,2,:) - logbasedpoweroutUnspec{iSet}(:,3,:);
% % %             percbasedpoweroutUnspec{iSet}(:,6,:) = percbasedpoweroutUnspec{iSet}(:,4,:) - percbasedpoweroutUnspec{iSet}(:,5,:);
% % %             percbasedpoweroutUnspec{iSet}(:,7,:) = percbasedpoweroutUnspec{iSet}(:,2,:) - percbasedpoweroutUnspec{iSet}(:,3,:);
% % %
% % %         catch
% % %             disp(' ..')
% % %         end
% % %     end
%
