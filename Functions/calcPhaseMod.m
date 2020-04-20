
function out = calcPhaseMod(F, dmat, varagin)
% [out] = calcSUmod(F, varargin)
% calculate event based phase modulation of spiking
%
% args:
% - F: struct with F.animal and F.data
%   - F.data: time, psth, instantFR
% - dmat: struct
%       - dayeps: [day ep;...] per ev
%       - dm: binary. (event x set)
%       - expvars: dmat labels as cell array of strings
%{
Saw
- barn:rat:beer:saw
@DKR
%}

pconf = paramconfig;
comuteShuf = 1;
numShuf = 1000;
bin = .001;
minILIthresh = .06; % seconds
maxILIthresh = .250; % seconds
numPhaseBins = 36;
try
    assign(varagin{:})
catch
end

for a = 1:length(F)
    animal = F(a).animal{3};
    out(a).animal = F(a).animal;
    fprintf('=========== %s ===========\n', animal);
    out(a).output = {};
    out(a).outputIdx = dmat(a).expvars;
    for iv = 1:size(dmat(a).expvars,2)
        out(a).output{iv} = [];
        for c = 1:length(F(a).output{1})
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
            
            % for each condition, get events' psth,
            eT = F(a).output{1}(c).eventTimes(dayDM(:,iv));
            time = F(a).output{1}(c).time;
            psth = F(a).output{1}(c).psth(dayDM(:,iv),:);
            % then collect spikes into inter event bins and compute their pct of interval
            ILI = diff(eT);
            cIdx = knnsearch(time', 0);
            ILIidx = knnsearch(time', ILI);
            pSinceLick = [];
            spkOffsetAll = [];
            for e = 1:length(ILI) % per event interval
                if ILI(e) < maxILIthresh && ILI(e) > minILIthresh
                    % get spikes between eventA and eventB
                    spIliIdx = find(psth(e,cIdx:ILIidx(e)));
                    if ~isempty(spIliIdx)
                        % get event interval percent of spikes
                        spkOffset = spIliIdx*bin; % idx distance from center (event), scaled to time
                        pSinceLick{e,1} = [spkOffset ./ ILI(e)]';
                        spkOffsetAll{e,1} = spkOffset;
                    end
                else
                    continue
                end
            end
            spikePctSinceLick = cell2mat(pSinceLick);
            if ~isempty(spikePctSinceLick)
                % convert percent into phase
                spikeIEIphase = 2*pi*spikePctSinceLick;                
                % get mean resultant vec mag and angle
                meanvec = mean(exp(1i*spikeIEIphase));
                meanMRVmag = abs(meanvec);
                vecang = angle(meanvec);
                % get log normalized Raleigh Z as 'phasemod'
                [pval, z] = circ_rtest(spikeIEIphase);
                phasemod = log(z);
                
                cF.pval  = pval;
                cF.ILI = ILI;
                cF.spkOffsetAll = spkOffsetAll;
                cF.spikePctSinceLick = spikePctSinceLick;
                cF.spikeIEIPhase = spikeIEIphase;
                cF.meanMRVmag = meanMRVmag;
                cF.vecang = vecang;
                cF.phasemod = phasemod;
                cF.area = idata.area;
                cF.subarea = idata.subarea;
                
%                 get phaseHist
                binphaseEdges = linspace(0, 2*pi, numPhaseBins);
                cf.spikeIEIPhaseHistProb = histcounts(spikeIEIphase, binphaseEdges, ...
                    'Normalization', 'probability');
%                 
                if comuteShuf
                    % shuffle percent of interval for all spikes
                    r = rand(length(spikePctSinceLick), numShuf);
                    cF.phasemodSh = [];
                    cF.vecangSh = [];
                    cF.meanMRVmagSh = [];
                    for s = 1:numShuf
                        spikeIEIphase = 2*pi*r(:,s);
                        meanvec = mean(exp(1i*spikeIEIphase));
                        cF.meanMRVmagSh = [cF.meanMRVmagSh; abs(meanvec)];
                        cF.vecangSh = [cF.vecangSh; angle(meanvec)];
                        [~, z] = circ_rtest(spikeIEIphase);
                        cF.phasemodSh = [cF.phasemodSh; log(z)];
                    end
                    % real-mod shuff-pct
                    cF.modPctRank = 100*(1-(sum(cF.phasemodSh > cF.phasemod)./numShuf));
                    fprintf('iv%d mod > %.02f pctShufs.\n', iv, cF.modPctRank)
                end
            end
            out(a).output{iv} = [out(a).output{iv}; cF];
        end
    end
end
end

function cf = init_out()
cf.pval = nan;
cf.ILI = nan;
cf.spkOffsetAll = nan;
cf.spikePctSinceLick = nan;
cf.spikeIEIPhase = nan;
cf.meanMRVmag = nan;
cf.vecang = nan;
cf.phasemod = nan;
cf.modPctRank = nan;

cf.area = '';
cf.subarea = '';

cf.spikeIEIPhaseHistProb = [];
cf.meanMRVmagSh = nan;
cf.vecangSh = nan;
cf.phasemodSh = nan;
end