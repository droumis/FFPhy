
%{
 select clusters and swr-designmat,
then cacluate swr mod, shuff, mod-shuffrank
%}

function modF = calcSUmod(F, filtF, varargin)
pconf = paramconfig;
% filtF filter to select animal, epochs, ntrode, clusers
respwin = [0 200]; % response period in ms rel to swr on
basewin = [-300 -100]; % baseline period in ms rel to swr on
minNumSwr = 10;
% minNumSwrSpikes = 10;
nshuffs = 1000;
shuffms = 700;
saveResult = 1;
filetail = '';
if ~isempty(varargin)
    assign(varargin{:})
end

modF = struct;
for ian = 1:length(F)
    animal = F(ian).animal;
    modF(ian).data = struct;
    cn = 0;
    % filter for certain clusters
    Fdtc = [];
    for ie = 1:size(filtF(ian).epochs{1},1)
        ieTC = filtF(ian).data{1}{ie};
        d = repmat(filtF(ian).epochs{1}(ie,1), size(ieTC,1),1);
        Fdtc = [Fdtc; d ieTC];
    end
    modFdata = F(ian).data(ismember(cell2mat({F(ian).data.dtc}'), Fdtc, 'rows'));
    
    for ic = 1:length(modFdata)
        % make swr-design mat for this cluster e.g. [allswr, lickbout-swrs]
        dtc = modFdata(ic).dtc;
        swrStart = modFdata(ic).eventtags(:,2); % col 2 is swr starttime
        epIdx = find(filtF(ian).epochs{1}(:,1) == dtc(1));
        DayexIntervals = [];
        for iep = 1:numel(epIdx)
            DayexIntervals = [DayexIntervals; filtF(ian).excludetime{1}{epIdx(iep)}];
        end
        dmat = logical([ones(size(swrStart,1),1) ~isExcluded(swrStart, DayexIntervals)]);
        if sum(dmat(:,2)) < minNumSwr % filter num lickbout-swrs
            fprintf('%d %d %d only %d swrs in lickbouts. skipping\n',dtc, ...
                sum(dmat(:,2)));
            continue
        end
        %% meanmod score for this cluster, per condition
        clear pctRespAboveShuff meanPctSwrResp pctRespAboveShuff
        for id = 1:size(dmat,2)
            ilickFRpsth = modFdata(ic).instantFR(dmat(:,id),:);
%             meanpsthFR{id} = nanmean(ilickFRpsth);
            time = modFdata(ic).time;
            % response
            respIdx = [knnsearch(time', respwin(1)/1000) knnsearch(time', respwin(2)/1000)];
            respmeanperSWR = nanmean(ilickFRpsth(:,respIdx(1):respIdx(2)),2)+1e-6;
            numSpikesResp{id} = sum(sum(modFdata(ic).psth(dmat(:,id),respIdx(1):respIdx(2))));
            % baseline
            baseIdx = [knnsearch(time', basewin(1)/1000) knnsearch(time', basewin(2)/1000)];
            basemeanperSWR = nanmean(ilickFRpsth(:,baseIdx(1):baseIdx(2)),2)+1e-6;
            numSpikesBase{id} = sum(sum(modFdata(ic).psth(dmat(:,id),baseIdx(1):baseIdx(2))));
            % mean pct change from baseline
            pctSwrResp = 100*(respmeanperSWR-basemeanperSWR)./basemeanperSWR;
            meanPctSwrResp{id} = mean(pctSwrResp);
            if numSpikesResp{id} < 100
                fprintf('%d %d %d only %d spikes in response period\n',dtc, ...
                numSpikesResp{id});
            end
            if numSpikesBase{id} < 100
                fprintf('%d %d %d only %d spikes in baseline period\n',dtc, ...
                numSpikesBase{id});
            end
            
            % Shuffle
            numSamples = size(ilickFRpsth,2);
            r = randi([-shuffms shuffms],size(ilickFRpsth,1), nshuffs); % shuff shift offsets
            shuffMeanPctSwrResp{id} = nan(nshuffs,1);
            %%
            tic
            for ish = 1:nshuffs
                % circshift 
%                 ishLickFRpsth = cell2mat(arrayfun(@(x,y) circshift(ilickFRpsth(x,:), ...
%                     r(x,ish),2),1:size(ilickFRpsth,1), 'uni',0)');
                respShmeanperSWR = nan(size(ilickFRpsth,1),1);
                respShmeanperSWR = nanmean(cell2mat(arrayfun(@(x,y) ...
                    ilickFRpsth(x,respIdx(1)+y:respIdx(2)+y), ...
                    1:size(ilickFRpsth,1), r(:,ish)', 'uni',0)'),2);
                baseShmeanperSWR = nan(size(ilickFRpsth,1),1);
                baseShmeanperSWR = nanmean(cell2mat(arrayfun(@(x,y) ...
                    ilickFRpsth(x,(baseIdx(1)+y):(baseIdx(2)+y)), ...
                    1:size(ilickFRpsth,1), r(:,ish)', 'uni',0)'),2);
                % instead of circshifting.. which i think is taking a long
                % time.. randomly select selection indices and just slice 
                % response
%                 respShmeanperSWR = nanmean(ishLickFRpsth(:,respIdx(1):respIdx(2)),2);
%                 % baseline
%                 baseShmeanperSWR = nanmean(ishLickFRpsth(:,baseIdx(1):baseIdx(2)),2);
                % mean pct change from baseline
                pctSwrResp = 100*(respShmeanperSWR-baseShmeanperSWR)./baseShmeanperSWR;
                pctSwrResp(pctSwrResp == inf) = [];
                shuffMeanPctSwrResp{id}(ish,1) = nanmean(pctSwrResp);
            end
            %%
            % real-mod shuff-pct
            pctRespAboveShuff{id} = 100*(1-(sum(shuffMeanPctSwrResp{id} > ...
                meanPctSwrResp{id})./nshuffs));
            fprintf('%s %d %d %d cond%d mod is > %.02f pct of shuffs. took %.03f seconds\n', animal, dtc,id,...
                pctRespAboveShuff{id}, toc)
        end
        %% ouput
        cn = cn +1;
        modF(ian).animal = animal;
        modF(ian).data(cn).dtc = dtc;
        modF(ian).data(cn).time = time;
        modF(ian).data(cn).inputdata = modFdata(ic); % save the parent data
        
        modF(ian).data(cn).dmat = dmat;
        modF(ian).data(cn).dmatIdx = {'all', 'lickbout'};
%         modF(ian).data(cn).meanpsthFR = meanpsthFR;
        modF(ian).data(cn).numSpikesResp = numSpikesResp;
        modF(ian).data(cn).numSpikesBase = numSpikesBase;
        modF(ian).data(cn).meanPctSwrResp = meanPctSwrResp;
        modF(ian).data(cn).shuffMeanPctSwrResp = shuffMeanPctSwrResp;
        modF(ian).data(cn).respAboveShuffPct = pctRespAboveShuff;
    end
end
% ---------------- Save result ---------------------------------------------------
if saveResult
    save_data(modF, 'results', 'sumod', 'filetail', filetail)
end
end
