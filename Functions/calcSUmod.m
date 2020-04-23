

function out = calcSUmod(F, dmat, varargin)
% [out] = calcSUmod(F, varargin)
% TODO: rename func to calcTimeModSU
% calculate event-triggered modulation of SU spiking
%
% args:
% - F: struct with F.animal and F.data
%   - F.data: time, psth, instantFR
% - dmat: struct (make with beer)
%       - dayeps: [day ep;...] per ev
%       - dm: binary. (event x set)
%       - expvars: dmat labels as cell array of strings
% varargs:
%{
Notes;
need to add this func to dfa_eventTrigSpiking..
or add that to this and make this a dfa..
%}

pconf = paramconfig;
rwin = [0 .2]; % response period in s rel to swr on
bwin = [-.7 -.2]; % baseline period in s rel to swr on
% minNumEvents = 10;
% minSpikesResp = 10;
% minSpikesBase = 10;
% minNumSwrSpikes = 10;
minSpikes = 50; % TODO.. need to implement
computeShuf = 0;
nshuffs = 1000;
shuffbyMax = .7; % s to rand shuff by
saveResult = 1;
filetail = '';
use_psfr = 1; % use firing rate. otherwise uses instantaneous firing rate (see rat)

try assign(varargin{:}); catch; end

out = struct;
for a = 1:length(F)
    animal = F(a).animal{3};
    out(a).animal = F(a).animal;
    fprintf('=========== %s ===========\n', animal);
    out(a).output = {};
    out(a).dmatIdx = dmat(a).expvars;
    for iv = 1:size(dmat(a).expvars,2)
        out(a).output{iv} = [];
        for c = 1:length(F(a).output{1}) % could parfor this for each cluster
            cF = init_out();
            cF.animal = animal;
            % collect event design mat for this cluster (per day)
            idata = F(a).output{1}(c);
            idx = idata.index;
            day = idx(1);
            nt = idx(2);
            cl = idx(3);
            eps = idx(4:end);
            fprintf('%s %d %d %d\n', animal, day, nt, cl);
            dayDM = ones(length(idata.eventTimes),1);
            if ~isempty(dmat)
                try
                    % get the dmat for this day's events
                    dayIdx = dmat(a).dayeps(:,1) == idx(1);
                    dayDM = logical(dmat(a).dm(dayIdx,:));
                catch
                    error('dmat needs valid dayeps and dm')
                end
            end
            cF.index = idx;
            if use_psfr
                time = idata.wtime;
                psfr = idata.psfr;
            else
                time = idata.time;
                psfr = idata.psifr;
            end
            cF.time = time;
            %% meanmod score for this cluster, per condition
            ivIdx = find(dayDM(:,iv));
            evIFR = psfr(ivIdx,:);
            cF.evMean = nanmean(evIFR,1); % full mean per DM.
            cF.evMeanZ = zscore(cF.evMean);
            
            % response mean FR
            rIdx = [knnsearch(time', rwin(1)) knnsearch(time', rwin(2))];
            evRespM = nanmean(evIFR(:,rIdx(1):rIdx(2)),2);
            
            psthRIdx = [knnsearch(idata.time', rwin(1)) knnsearch(idata.time', rwin(2))];
            numRSpikes = sum(sum(idata.psth(ivIdx,psthRIdx(1):psthRIdx(2))));
            
            fprintf('spikes in response: %d \n', numRSpikes);
            %             if numRSpikes < minSpikesResp
            %                 fprintf('skipping\n')
            %                 continue
            %             end
            % baseline mean FR
            bIdx = [knnsearch(time', bwin(1)) knnsearch(time', bwin(2))];
            evBaseM = nanmean(evIFR(:,bIdx(1):bIdx(2)),2);
            evBaseSTD = nanstd(evIFR(:,bIdx(1):bIdx(2)),[],2);
            
            psthbIdx = [knnsearch(idata.time', bwin(1)) knnsearch(idata.time', bwin(2))];
            numBSpikes = sum(sum(idata.psth(dayDM(:,iv), psthbIdx(1):psthbIdx(2))));
            
            fprintf('spikes in baseline period: %d \n', numBSpikes);
            if (numBSpikes + numRSpikes) < minSpikes
                fprintf('skipping\n')
                continue
            end
            
            % take the diff of base and resp period for each event
            % mean pct change from baseline
            useB = ~(evBaseM == 0);
            mPctChange = nanmean((evRespM(useB)-evBaseM(useB))./evBaseM(useB)) * 100;
            useBs = ~(evBaseSTD == 0);
            mZChange = nanmean((evRespM(useBs)-evBaseM(useBs))./evBaseSTD(useBs));
            
            evIFRz = zscore(evIFR, [], 2);
            mZResp = nanmean(nanmean(evIFRz(:,rIdx(1):rIdx(2)),2));

            %% Shuffle
            if computeShuf
                numEv = size(evIFR,1);
                shuffbyMaxBins = knnsearch(time', shuffbyMax)-knnsearch(time', 0);
                r = randi([-shuffbyMaxBins shuffbyMaxBins], numEv, nshuffs);
                mPctChangeSh = nan(nshuffs,1);
                mZChangeSh = nan(nshuffs,1);
                % use parfor at unit level above to reduce worker data copies
                for i = 1:nshuffs % can use parfor. 
                    % for each event, select a shuf rand period as r and b
                    % time-mean fr in response period
                    evRespMSh = nanmean(cell2mat(arrayfun(@(x,y) evIFR(x,rIdx(1)+y:rIdx(2)+y), ...
                        1:size(r,1), r(:,i)', 'uni',0)'),2);
                    evBaseMSh = nanmean(cell2mat(arrayfun(@(x,y) ...
                        evIFR(x,bIdx(1)+y:bIdx(2)+y), 1:size(r,1), r(:,i)', 'uni',0)'),2);
                    evBaseSTDSh = nanstd(cell2mat(arrayfun(@(x,y) ...
                        evIFR(x,bIdx(1)+y:bIdx(2)+y), 1:size(r,1), r(:,i)', 'uni',0)'),[],2);
                    % mean pct change from baseline
                    useB = ~(evBaseMSh == 0);
                    mPctChangeSh(i) = nanmean((evRespMSh(useB)-evBaseMSh(useB))./evBaseMSh(useB))*100;
                    useBs = ~(evBaseSTDSh == 0);
                    mZChangeSh(i) = nanmean((evRespMSh(useBs)-evBaseMSh(useBs))./evBaseSTDSh(useBs));
                end
                % real-mod shuff-pct
                modPctRank = 100*(1-(sum(mPctChangeSh > mPctChange)./nshuffs));
                modZRank = 100*(1-(sum(mZChangeSh > mZChange)./nshuffs));
                fprintf('iv%d Pctmod > %.02f pctShuffs.\n', iv, modPctRank)
                fprintf('iv%d Zmod > %.02f pctShuffs.\n', iv, modZRank)
            end
            %% unit output
            cF.mPctChange = mPctChange;
            cF.mPctChangeSh = mPctChangeSh;
            cF.modPctRank = modPctRank;
            cF.mZChangeSh = mZChangeSh;
            cF.mZChange = mZChange;
            cF.modZRank = modZRank;
            cF.mZResp = mZResp;
            cF.area = idata.area;
            cF.subarea = idata.subarea;
            cF.cellInfo = idata.cellInfo;
        end
        out(a).output{iv} = [out(a).output{iv}; cF];
    end
end
end

function cF = init_out()
% are empty lists the right init for all these? phasmod uses mostly nan
cF.mZChange = [];
cF.animal = '';
cF.dmat = [];
cF.dmatIdx = [];
cF.index = [];
cF.time = [];

cF.evMean = [];
cF.evMeanZ = [];
% out.numRSpikes = [];
% out.numBSpikes = [];
cF.mPctChange = [];
cF.mPctChangeSh = [];
cF.modPctRank = [];

cF.mZChangeSh = [];
cF.mZChange = [];
cF.modZRank = [];

cF.mZResp = [];

cF.area = '';
cF.subarea = '';
cF.cellInfo = [];
end
