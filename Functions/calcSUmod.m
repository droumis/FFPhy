

function out = calcSUmod(F, varargin)
% [out] = calcSUmod(F, varargin)
% calculate event-triggered modulation of SU spiking
%
%                         Wheedling Wheelbarrow
%                                               _______
%                  ___________________________.'.------`
%                 '---------------------------.'
%                   `.                      .'
%                 .-//`.                  .'
%              .' .//.'/`================'
%             =[=:====:=]=           \\||
%              '. `--' .'             \_|
%                `-  -'
% args:
% - F: struct with F.animal and F.data
%   - F.data: time, psth, instantFR
%
% varargs:
% - dmat: struct (make with beer)
%       - dayeps: [day ep;...] per ev 
%       - dm: binary. (event x set)

%{
Notes:
    - barn:rat:beer:wheelbarrow

need to add just add this func to dfa_eventTrigSpiking.. 
or add that to this and make this a dfa.. 

FFPhy V0.1
@DR
%}

pconf = paramconfig;
% filtF filter to select animal, epochs, ntrode, clusers
rwin = [0 .1]; % response period in s rel to swr on
bwin = [-.7 -.2]; % baseline period in s rel to swr on
% minNumEvents = 10;
% minSpikesResp = 10;
% minSpikesBase = 10;
% minNumSwrSpikes = 10;
minSpikes = 50;
% runShuffle = 1;
nshuffs = 100;
shuffbyMax = .7; % s to rand shuff by
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
        OP.animal = animal;
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
            % get the dmat for this day's events
            dayIdx = dmat(a).dayeps(:,1) == idx(1);
            dayDM = dmat(a).dm(dayIdx,:);            
        end
        OP.index = idx;
        use_psfr = 1;
        if use_psfr
            time = idata.wtime;
            psfr = idata.psfr;
        else
            time = idata.time;
            psfr = idata.psifr;
        end
        OP.time = time;
        %% meanmod score for this cluster, per condition
        % right now i'm just using the instantFR for computing mod, but num
        % spikes for thresholding within win is from the psth
        for iv = 1:size(dayDM,2)
%             try
            ivIdx = find(dayDM(:,iv));
            evIFR = psfr(ivIdx,:);
            OP.evMean{iv} = nanmean(evIFR,1); % full mean per DM.
            OP.evMeanZ{iv} = zscore(OP.evMean{iv}); % full mean per DM.
%                   evIFR = idata.psthW(dayDM(:,iv),:);
%             catch
%                 fprintf('no valid events\n')
%                 continue
%             end
            
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
%             OP.numRSpikes{iv} = numRSpikes;
%             OP.numBSpikes{iv} = numBSpikes;
            %% Shuffle
            numEv = size(evIFR,1);
            binsize = diff(time(1:2));
            shuffbyMaxBins = knnsearch(time', shuffbyMax)-knnsearch(time', 0);
            r = randi([-shuffbyMaxBins shuffbyMaxBins], numEv, nshuffs); 
            mPctChangeSh = nan(nshuffs,1);
            for i = 1:nshuffs % can use parfor
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
           %%
            % real-mod shuff-pct
            modPctRank = 100*(1-(sum(mPctChangeSh > mPctChange)./nshuffs));
            modZRank = 100*(1-(sum(mZChangeSh > mZChange)./nshuffs));
            fprintf('iv%d Pctmod > %.02f pctShuffs.\n', iv, modPctRank)
            fprintf('iv%d Zmod > %.02f pctShuffs.\n', iv, modZRank)
            
            OP.mPctChange{iv} = mPctChange;
            OP.mPctChangeSh{iv} = mPctChangeSh;
            OP.modPctRank{iv} = modPctRank;
            
            OP.mZChangeSh{iv} = mZChangeSh;
            OP.mZChange{iv} = mZChange;
            OP.modZRank{iv} = modZRank;
            
            OP.mZResp{iv} = mZResp;
            
            OP.area = idata.area;
            OP.subarea = idata.subarea;
        end
        ic = ic +1;
        out(a).output{1}(ic) = OP;
    end
end
end

function out = init_out()
out.mZChange = [];
out.animal = '';
out.dmat = [];
out.dmatIdx = [];
out.index = [];
out.time = [];

out.evMean = [];
out.evMeanZ = [];
% out.numRSpikes = [];
% out.numBSpikes = [];
out.mPctChange = [];
out.mPctChangeSh = [];
out.modPctRank = [];

out.mZChangeSh = [];
out.mZChange = [];
out.modZRank = [];

out.mZResp = [];

out.area = '';
out.subarea = '';

end
