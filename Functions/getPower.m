

function getPower(ripstate, Fp, varargin)
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
    pwrout = struct;
    wp = getWaveParams(Fp.waveSet, []);
    animal = Fp.animals{ian};
    %     ntrodes = lfpstack(ian).ntrodes;
    andef = animaldef(animal);
    tetinfo = loaddatastruct(andef{2}, animal, 'tetinfo');
    tmp = cellfetch(tetinfo, '');
    ntrodes = unique(tmp.index(:,3));
    
    for itype = 1:length(lfptypes)
        lfptype = lfptypes{itype};
        meanpwr = cell(1,length(ntrodes));
        medianpwr = cell(1,length(ntrodes));
        meanZpwr = cell(1,length(ntrodes));
        parfor nti = 1:length(ntrodes) % use parfor
            meandbpwr = cell(1,length(ripstate(ian).statesetsfields));
            mediandbpwr = cell(1,length(ripstate(ian).statesetsfields));
            meanZpower = cell(1,length(ripstate(ian).statesetsfields));
            nt = ntrodes(nti);
            as = loadAS(animal, nt, Fp.waveSet, lfptype);
            
            for stset = 1:length(ripstate(ian).statesetsfields)
                stidx = find(ripstate(ian).statesets(:,stset));
                if ~isempty(wp.baseind)
                    basidx = wp.baseind;
                else
                    basidx = as.analyticsignal.waveparams.baseind;
                end
                [meandbpwr{stset}, mediandbpwr{stset}, meanZpower{stset}] = computePower(...
                    as.analyticsignal.analyticsignal(:,stidx,:), ...
                    'baselIdx', basidx);
                % for each condition, ntrode, run permutation test to get
                % significant blob mask.. so i was only doing the perm test
                % for the diff sets before.. so i need to figure out the
                % right way to do single set perm testing..
                %                 powerpermtest(dbpwr{stset}, stidx
            end
            % now run the set diff for outb - inb, correct - error, etc..
            % with perm testing..
            meanpwr{nti} = meandbpwr;
            medianpwr{nti} = mediandbpwr;
            meanZpwr{nti} = meanZpower;
        end
        pwrout.animal = animal;
        pwrout.ripstate = ripstate;
        pwrout.wp = wp;
        pwrout.Fp = Fp;
        pwrout.lfptype = lfptype;
        
        pwrout.meandbpower = meanpwr;
        pwrout.mediandbpower = medianpwr;
        pwrout.meanZpower = meanZpwr;
        
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

function [meandbpower, mediandbpower, meanZpower] = computePower(as,varargin)

baseind = [2 3]*round(6001/6);
if ~isempty(varargin)
    assign(varargin{:});
end
% AS dims = samples X events X freqs
power = mean(abs(as).^2,2);
baseline_power = mean(mean(abs( as(baseind(1):baseind(2),:,:)) .^2,2),1);
baseline_std_power = std(mean(abs( as(baseind(1):baseind(2),:,:)) .^2,2),[],1);
% pctpower = 100 * bsxfun(@rdivide, bsxfun(@minus, power, basepower), basepower);
meandbpower = 10*log10(bsxfun(@rdivide, power, baseline_power));
meanZpower = (power-baseline_power) ./ baseline_std_power;

medianpower = median(abs(as).^2,2);
medianbasepower = median(median(abs( as(baseind(1):baseind(2),:,:)) .^2,2),1);
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
