

function out = calcSUmod(F, varargin)
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
% minNumEvents = 10;
minSpikesResp = 10;
minSpikesBase = 10;
% minNumSwrSpikes = 10;
% runShuffle = 1;
nshuffs = 200;
shuffms = 700;
saveResult = 1;
filetail = '';
dmat = [];
dmatIdx = {'all'};
if ~isempty(varargin)
    assign(varargin{:})
end

out = struct;
for a = 1:length(F)
    animal = F(a).animal{3};
    out(a).animal = F(a).animal;
    fprintf('=========== %s ===========\n', animal);
    out(a).output = {};
    out(a).dmatIdx = dmatIdx;
    ic = 0;
    for c = 1:length(F(a).output{1}) % for each cluster
        OP = init_out();
        % collect event design mat for this cluster (per day)
        idata = F(a).output{1}(c);
        idx = idata.index;
        day = idx(1);
        nt = idx(2);
        cl = idx(3);
        eps = idx(4:end);
        fprintf('%d %d %d\n', day, nt, cl);
        dayDM = ones(length(idata.eventTimes),1);
        if ~isempty(dmat)
            % get the dmat for this day's events
            dayIdx = dmat(a).dayeps(:,1) == idx(1);
            dayDM = dmat(a).dm(dayIdx,:);            
        end
        OP.index = idx;
        time = idata.time;
        OP.time = time;
        %% meanmod score for this cluster, per condition
        % right now i'm just using the instantFR for computing mod, but num
        % spikes for thresholding within win is from the psth
        for iv = 1:size(dayDM,2)
            try
                evIFR = idata.psfr(dayDM(:,iv),:);
            catch
                fprintf('no valid events\n')
                continue
            end
            % meanpsthFR{id} = nanmean(ilickFRpsth);
            % response mean FR
            rIdx = [knnsearch(time', rwin(1)/1000) knnsearch(time', rwin(2)/1000)];
            evRespM = nanmean(evIFR(:,rIdx(1):rIdx(2)),2)+1e-6;
            numRSpikes = sum(sum(idata.psth(dayDM(:,iv),rIdx(1):rIdx(2))));
            fprintf('spikes in response: %d \n', numRSpikes);
            if numRSpikes < minSpikesResp
                fprintf('skipping\n')
                continue
            end
            % baseline mean FR
            bIdx = [knnsearch(time', bwin(1)/1000) knnsearch(time', bwin(2)/1000)];
            evBaseM = nanmean(evIFR(:,bIdx(1):bIdx(2)),2)+1e-6;
            numBSpikes = sum(sum(idata.psth(dayDM(:,iv),bIdx(1):bIdx(2))));
            fprintf('spikes in baseline period: %d \n', numBSpikes);
            if numBSpikes < minSpikesBase
                fprintf('skipping\n')
                continue
            end
            % mean pct change from baseline
            mPctChange = nanmean((evRespM-evBaseM)./evBaseM*100);
            
            OP.evMean{iv} = nanmean(evIFR,1); % full mean per DM.
            OP.numRSpikes{iv} = numRSpikes;
            OP.numBSpikes{iv} = numBSpikes;
            OP.mPctChange{iv} = mPctChange;
            
            %% Shuffle
            numEv = size(evIFR,1);
            r = randi([-shuffms shuffms], numEv, nshuffs); % shuff shift offsets
            mPctChangeSh = nan(nshuffs,1);
            for i = 1:nshuffs % can use parfor
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
            fprintf('iv%d mod > %.02f pctShuffs.\n', iv, modPctRank)
            
            OP.mPctChangeSh{iv} = mPctChangeSh;
            OP.modPctRank{iv} = modPctRank;
            OP.area = idata.area;
            OP.subarea = idata.subarea;
        end
        ic = ic +1;
        out(a).output{1}(ic) = OP;
    end
end
end

function op = init_out()
op.dmat = [];
op.dmatIdx = [];
op.index = [];
op.time = [];

op.evMean = [];

op.numRSpikes = [];
op.numBSpikes = [];
op.mPctChange = [];

op.mPctChangeSh = [];
op.modPctRank = [];
op.area = '';
op.subarea = '';

end
