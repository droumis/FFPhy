

function modF = calcSUmod(dmat, F, varargin)
%{
need to add this to dfa_eventTrigSpiking.. or add that to this and make
this a dfa.. ... should be called with singleDayCellAnal and somehow
receive a dmat.. 

calculate event-triggered modulation of SU spiking
F is a struct array with F.animal and F.data
F.data contains at least fields:
- time, psth, instantFR
%}

pconf = paramconfig;
% filtF filter to select animal, epochs, ntrode, clusers
rwin = [0 200]; % response period in ms rel to swr on
bwin = [-300 -100]; % baseline period in ms rel to swr on
minNumEvents = 10;
minSpikesResp = 100;
minSpikesBase = 100;
% minNumSwrSpikes = 10;
runShuffle = 1;
nshuffs = 200;
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
    modF(ian).animal = animal;
    cn = 1;
    
    for ic = 1:length(F(ian).data)
        % collect event design mat for this cluster (per day)
        idx = F(ian).data.index;
        fprintf('%d %d %d\n', idx);
        dayIdx = dmat(ian).dayeps(:,1) == idx(1);
        dayDM = dmat(ian).dm(dayIdx,:);
        dayET = dmat(ian).evStart(dayIdx,:);
        
        modF(ian).data(cn).index = idx;
        time = F(ian).data(ic).time;
        modF(ian).data(cn).time = time;
        modF(ian).data(cn).inpData = F(ian).data(ic); % save the parent data
        modF(ian).data(cn).dmat = dmat;
        modF(ian).data(cn).dmatIdx = {'all', 'lickbout'};

        %% meanmod score for this cluster, per condition
        % right now i'm just using the instantFR for computing mod, but num
        % spikes for thresholding within win is from the psth
        for iv = 1:numel(dmat(ian).expvars)
            try
                evIFR = F(ian).data(ic).instantFR(dayDM(:,iv),:);
            catch
                fprintf('no valid events\n')
                continue
            end
            % meanpsthFR{id} = nanmean(ilickFRpsth);
            % response mean FR
            rIdx = [knnsearch(time', rwin(1)/1000) knnsearch(time', rwin(2)/1000)];
            evRespM = nanmean(evIFR(:,rIdx(1):rIdx(2)),2)+1e-6;
            numRSpikes{iv} = sum(sum(F(ian).data(ic).psth(dayDM(:,iv),rIdx(1):rIdx(2))));
            fprintf('spikes in response: %d \n', numRSpikes{iv});
            if numRSpikes{iv} < minSpikesResp
                fprintf('skipping\n')
                continue
            end
            
            % baseline mean FR
            bIdx = [knnsearch(time', bwin(1)/1000) knnsearch(time', bwin(2)/1000)];
            evBaseM = nanmean(evIFR(:,bIdx(1):bIdx(2)),2)+1e-6;
            numBSpikes{iv} = sum(sum(F(ian).data(ic).psth(dayDM(:,iv),bIdx(1):bIdx(2))));
            fprintf('spikes in baseline period: %d \n', numBSpikes{iv});
            if numBSpikes{iv} < minSpikesBase
                fprintf('skipping\n')
                continue
            end
            
            % mean pct change from baseline
            mPctChange = nanmean((evRespM-evBaseM)./evBaseM*100);
            modF(ian).data(cn).numRSpikes{iv} = numRSpikes;
            modF(ian).data(cn).numBSpikes{iv} = numBSpikes;
            modF(ian).data(cn).mPctChange{iv} = mPctChange;
            
            %% Shuffle
            numEv = size(evIFR,1);
            r = randi([-shuffms shuffms], numEv, nshuffs); % shuff shift offsets
            mPctChangeSh = nan(nshuffs,1);
            parfor i = 1:nshuffs
                % for each event, select a shuf shifted slice as r and b, then mean
                evRespMSh = nanmean(cell2mat(arrayfun(@(x,y) ...
                    evIFR(x,rIdx(1)+y:rIdx(2)+y), 1:numEv, r(:,i)', 'uni',0)'),2);
                %                     baseShmeanperSWR = nan(size(evIFR,1),1);
                evBaseMSh = nanmean(cell2mat(arrayfun(@(x,y) ...
                    evIFR(x,bIdx(1)+y:bIdx(2)+y), 1:numEv, r(:,i)', 'uni',0)'),2);
                % mean pct change from baseline
                pctSwrResp = (evRespMSh-evBaseMSh)./evBaseMSh*100;
                pctSwrResp(pctSwrResp == inf) = [];
                mPctChangeSh(i) = nanmean(pctSwrResp);
            end
           %%
            % real-mod shuff-pct
            modPctRank = 100*(1-(sum(mPctChangeSh > mPctChange)./nshuffs));
            fprintf('%d %d %d cond%d mod is > %.02f pct of shuffs.\n',idx, iv, modPctRank)
            
            modF(ian).data(cn).mPctChangeSh{iv} = mPctChangeSh;
            modF(ian).data(cn).modPctRank{iv} = modPctRank;
        end
        cn = cn +1;
    end
end
% ---------------- Save result ---------------------------------------------------
if saveResult
    save_data(modF, 'results', 'sumod', 'filetail', filetail)
end
end
