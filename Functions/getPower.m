

function getPower(as, ripstate, Fp, varargin)
% LOAD ANALYTIC SIGNAL PER NTRODE AND COMPUTE POWER
% Fp = filter params (see load_filter_params)
% p = wave params (see getWaveParams)
% TODO: only uses lfpstack to get the list of ntrodes.. so i could do without it

% pwrout is a struct array per animal with dbpower data
% pwrout(an).dbpower is a cell array per ntrode of a cell array per statetype

lfptypes = {'eeg'};
savepower = 1;
if ~isempty(varargin)
   assign(varargin{:}); 
end

for ian = 1:length(Fp.animals)
    wp = getWaveParams(Fp.waveSet, []);
    animal = Fp.animals{ian};
    andef = animaldef(animal);
    tetinfo = loaddatastruct(andef{2}, animal, 'tetinfo');
    tmp = cellfetch(tetinfo, '');
    ntrodes = unique(tmp.index(:,3));

    for itype = 1:length(lfptypes)
        pwrout = struct;
        lfptype = lfptypes{itype};
        meanpwrNTs = cell(1,length(ntrodes));
        medianpwrNTs = cell(1,length(ntrodes));
%         meanZpwr = cell(1,length(ntrodes));
        zmappwrNTs = cell(1,length(ntrodes));
        threshmeanpwrNTs = cell(1,length(ntrodes));
        parfor nti = 1:length(ntrodes) % use parfor
            meandbpwrSTs = cell(1,length(ripstate(ian).statesetsfields));
            mediandbpwrSTs = cell(1,length(ripstate(ian).statesetsfields));
            zmapSTs = cell(1,length(ripstate(ian).statesetsfields));
            threshmeanSTs = cell(1,length(ripstate(ian).statesetsfields));
            nt = ntrodes(nti);
%             as = loadAS(animal, nt, Fp.waveSet, lfptype);
            % for each state/condition, compute power and run perm test vs baseline
            for stset = 1:length(ripstate(ian).statesetsfields)
                stidx = find(ripstate(ian).statesets(:,stset));
                [meandbpwrSTs{stset}, mediandbpwrSTs{stset}, zmapSTs{stset}, threshmeanSTs{stset}] = ...
                    computePower(as{itype}{nti}.as(:,stidx,:), wp, 'baselIdx', wp.baseind);
            end
            meanpwrNTs{nti} = meandbpwrSTs;
            medianpwrNTs{nti} = mediandbpwrSTs;
%             meanZpwr{nti} = meanZpower;
            zmappwrNTs{nti} = zmapSTs;
            threshmeanpwrNTs{nti} = threshmeanSTs;
        end
        pwrout.animal = animal;
        pwrout.ripstate = ripstate;
        pwrout.wp = wp;
        pwrout.Fp = Fp;
        pwrout.lfptype = lfptype;
        
        pwrout.meandbpower = meanpwrNTs;
        pwrout.mediandbpower = medianpwrNTs;
        pwrout.zmap = zmappwrNTs;
        pwrout.threshmean = threshmeanpwrNTs;
%         pwrout.meanZpower = meanZpwr;
        
        if savepower
            save_power(pwrout, animal, Fp.waveSet, lfptype)
        end
    end
end
    
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

function [meanpowerdb, mediandbpower, zmap, threshmean] = computePower(as,wp,varargin)

% baseind = [2 3]*round(6001/6);
run_permutation_test = 1;
voxel_pval   = 0.01;
if ~isempty(varargin)
    assign(varargin{:});
end
% AS dims = samples X events X freqs
power = abs(as).^2;
meanpower = mean(power,2);
baseline_power_mean = mean(mean(abs( as(wp.baseind(1):wp.baseind(2),:,:)) .^2,2),1);
% baseline_std_power = std(mean(abs( as(baseind(1):baseind(2),:,:)) .^2,2),[],1);
% pctpower = 100 * bsxfun(@rdivide, bsxfun(@minus, power, basepower), basepower);
meanpowerdb = 10*log10(bsxfun(@rdivide, meanpower, baseline_power_mean));
% meanZpower = (power-baseline_power) ./ baseline_std_power;

if run_permutation_test
    permuted_vals = zeros(numel(wp.timeWin), wp.n_permutes, wp.numfrex);
    %     powercut = power(1000:5000, :, :);
    for permi=1:wp.n_permutes
        fprintf(' %d ', permi);
        cutpoint = randsample(2:length(power(:,1,1))-diff(wp.baseind)-2,1);
        permuted_vals(:,permi,:) = 10*log10(bsxfun(@rdivide, ...
            mean(power([cutpoint:end 1:cutpoint-1],:,:),2),baseline_power_mean) );
    end
    zmap = squeeze((meanpowerdb-mean(permuted_vals,2)) ./ std(permuted_vals,[],2));
    threshmean = squeeze(meanpowerdb);
    threshmean(abs(zmap)<norminv(1-voxel_pval))=0;
end

medianpower = median(abs(as).^2,2);
medianbasepower = median(median(abs( as(wp.baseind(1):wp.baseind(2),:,:)) .^2,2),1);
mediandbpower = 10*log10(bsxfun(@rdivide, medianpower, medianbasepower));

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
