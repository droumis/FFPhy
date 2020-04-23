
function out = calcPhaseMod(F, dmat, varagin)
% [out] = calcSUmod(F, varargin)
% TODO: rename func to calcPhaseModSU
% calculate event based phase modulation of spiking
%
% args:
% - F: struct with F.animal and F.data
%   - F.data: time, psth, instantFR
% - dmat: struct
%       - dayeps: [day ep;...] per ev
%       - dm: binary. (event x set)
%       - expvars: dmat labels as cell array of strings

pconf = paramconfig;
comuteShuf = 1;
numShuf = 1000;
bin = .001;
minILIthresh = .06; % seconds
maxILIthresh = .250; % seconds
numPhaseBins = 36;
minILIspikes = 100; % TODO: implement min ILI spikes cumulative per unit

try assign(varagin{:}); catch; end

for a = 1:length(F)
    animal = F(a).animal{3};
    out(a).animal = F(a).animal;
    fprintf('=========== %s ===========\n', animal);
    out(a).output = {};
    out(a).outputIdx = dmat(a).expvars;
    for iv = 1:size(dmat(a).expvars,2) % for each experimental condition
        out(a).output{iv} = [];
        for c = 1:length(F(a).output{1}) % could parfor this for each unit
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
            % for each condition, get events' psth
            eT = F(a).output{1}(c).eventTimes(dayDM(:,iv));
            time = F(a).output{1}(c).time;
            psth = F(a).output{1}(c).psth(dayDM(:,iv),:);
            % collect and transform spikes into inter event bins and 
            % compute spike bin pct of interval
            ILI = diff(eT); % diff gives duration from each lick to the next
            cIdx = knnsearch(time', 0); % find center (event time index)
            ILIidx = knnsearch(time', ILI); % get index into time for each ILI duration
            pSinceLick = [];
            spkOffsetAll = [];
            for e = 1:length(ILI) % per event interval
                if ILI(e) < maxILIthresh && ILI(e) > minILIthresh % if valid ILI
                    % get spikes within current ILI
                    spIliIdx = find(psth(e,cIdx:ILIidx(e)));
                    if ~isempty(spIliIdx)
                        % get event interval percent of spikes
                        spkOffset = spIliIdx*bin; % idx distance from center (event), scaled to time
                        pSinceLick{e,1} = [spkOffset ./ ILI(e)]'; % spiketime pct of ILI
                        spkOffsetAll{e,1} = spkOffset;
                    end
                else
                    continue
                end
            end
            spikePctSinceLick = cell2mat(pSinceLick);
            if ~isempty(spikePctSinceLick)
                % convert percent into phase. since pct is 0-1, just 2*pi*pct
                spikeIEIphase = 2*pi*spikePctSinceLick;                
                % get mean resultant vec mag and angle *****
                meanvec = mean(exp(1i*spikeIEIphase));
                meanMRVmag = abs(meanvec); % magnitude
                vecang = angle(meanvec); % angle of mean vec
                % get log normalized Raleigh Z as 'phasemod'
                [pval, z] = circ_rtest(spikeIEIphase);
                phasemod = log(z); % log normalize    
                % make phase histogram edges
                binphaseEdges = linspace(0, 2*pi, numPhaseBins);
                % compute phase histogram *****
                spikeIEIPhaseHistProb = histcounts(spikeIEIphase, binphaseEdges, ...
                    'Normalization', 'probability');
                %% Shuffle
                if comuteShuf
                    % rand shuffle percent of interval for all spikes
                    % TODO: shuffle spikeIEIphase instead of spikePctSinceLick
                    r = rand(length(spikePctSinceLick), numShuf);
                    meanMRVmagSh = nan(numShuf);
                    vecangSh = nan(numShuf);
                    phasemodSh = nan(numShuf);
                    for s = 1:numShuf
                        spikeIEIphase = 2*pi*r(:,s);
                        meanvec = mean(exp(1i*spikeIEIphase));
                        meanMRVmagSh(s) = abs(meanvec);
                        vecangSh(s) = angle(meanvec);
                        [~, z] = circ_rtest(spikeIEIphase);
                        phasemodSh(s) = log(z);
                    end
                    % real-phasemod vs shuff-phasemod distribution
                    modPctRank = 100*(1-(sum(phasemodSh > phasemod)./numShuf));
                    fprintf('iv%d mod > %.02f pctShufs.\n', iv, modPctRank)
                    % Shuff output
                    cF.meanMRVmagSh = meanMRVmagSh;
                    cF.vecangSh = vecangSh;
                    cF.phasemodSh = phasemodSh;
                    cF.modPctRank = modPctRank;
                end
                %% unit output
                cF.area = idata.area;
                cF.subarea = idata.subarea;
                cF.ILI = ILI;
                cF.spkOffsetAll = spkOffsetAll;
                cF.spikePctSinceLick = spikePctSinceLick;
                cF.spikeIEIPhase = spikeIEIphase;
                cF.meanMRVmag = meanMRVmag;
                cF.vecang = vecang;
                cF.pval  = pval; % circ_r
                cF.phasemod = phasemod;
                cF.spikeIEIPhaseHistProb = spikeIEIPhaseHistProb;
            end
            out(a).output{iv} = [out(a).output{iv}; cF];
        end
    end
end
end

function cf = init_out()
% are nans the right init for all these? timemod uses mostly []
cf.area = '';
cf.subarea = '';
cf.ILI = [];
cf.spkOffsetAll = [];
cf.spikePctSinceLick = [];
cf.spikeIEIPhase = [];
cf.meanMRVmag = nan;
cf.vecang = nan;
cf.phasemod = nan;
cf.modPctRank = nan;
cf.spikeIEIPhaseHistProb = [];
cF.phasemodSh = [];
cF.vecangSh = [];
cF.meanMRVmagSh = [];
cf.pval = nan;
end