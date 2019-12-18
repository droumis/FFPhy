
function out = calcPhaseMod(F, dmat, varagin)
% [out] = calcSUmod(F, varargin)
% calculate event based phase modulation of SU spiking
%
%                   Sagacious Saw
%                                                __
%                                    _____....--' .'
%                          ___...---'._ o      -`(
%                ___...---'            \   .--.  `\
%      ___...---'                      |   \   \ `|
%     |                                |o o |  |  |
%     |                                 \___'.-`.  '.
%     |                                      |   `---'
%     '^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^'
%
% args:
% - F: struct with F.animal and F.data
%   - F.data: time, psth, instantFR
% - dmat: struct (make with beer)
%       - dayeps: [day ep;...] per ev
%       - dm: binary. (event x set)
%       - expvars: dmat labels as cell array of strings

%{
Notes:
    - barn:rat:beer:saw

@DKR
%}

comuteShuf = 1;
numShuf = 1000;
pconf = paramconfig;
bin = .001;
minILIthresh = .06; % seconds
maxILIthresh = .250; % seconds
try
    assign(varagin{:})
catch
end

for a = 1:length(F)
    animal = F(a).animal{3};
    out(a).animal = F(a).animal;
    fprintf('=========== %s ===========\n', animal);
    out(a).output = {};
    out(a).dmatIdx = dmat(a).expvars;
    ic = 0;
    for c = 1:length(F(a).output{1})
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
            try
                % get the dmat for this day's events
                dayIdx = dmat(a).dayeps(:,1) == idx(1);
                dayDM = logical(dmat(a).dm(dayIdx,:));
            catch
                error('dmat needs valid dayeps and dm')
            end
        end
        OP.index = idx;
        for iv = 1:size(dayDM,2)
            eT = F(a).output{1}(c).eventTimes(dayDM(:,iv));
            time = F(a).output{1}(c).time;
            psth = F(a).output{1}(c).psth(dayDM(:,iv),:);
            ILI = diff(eT);
            cIdx = knnsearch(time', 0);
            ILIidx = knnsearch(time', ILI);
            pSinceLick = [];
            spkOffsetAll = [];
            for e = 1:length(ILI)
                if ILI(e) < maxILIthresh && ILI(e) > minILIthresh
                    spIliIdx = find(psth(e,cIdx:ILIidx(e)));
                    if ~isempty(spIliIdx)
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
                spikeILIphase = 2*pi*spikePctSinceLick;
                meanvec = mean(exp(1i*spikeILIphase));
                meanMRVmag = abs(meanvec);
                vecang = angle(meanvec);
                [pval, z] = circ_rtest(spikeILIphase);
                phasemod = log(z);
                OP.pval{iv} = pval;
                OP.ILI{iv} = ILI;
                OP.spkOffsetAll{iv} = spkOffsetAll;
                OP.spikePctSinceLick{iv} = spikePctSinceLick;
                OP.spikeLickPhase{iv} = spikeILIphase;
                OP.meanMRVmag{iv} = meanMRVmag;
                OP.vecang{iv} = vecang;
                OP.phasemod{iv} = phasemod;
                
                OP.area = idata.area;
                OP.subarea = idata.subarea;
                if comuteShuf
                    r = rand(length(spikePctSinceLick), numShuf);
                    for s = 1:numShuf
                        spikeILIphase = 2*pi*r(:,s);
                        meanvec = mean(exp(1i*spikeILIphase));
                        meanMRVmag = abs(meanvec);
                        vecang = angle(meanvec);
                        [~, z] = circ_rtest(spikeILIphase);
                        phasemodSh{s} = log(z);
                    end
                    OP.meanMRVmagSh{iv} = meanMRVmag;
                    OP.vecangSh{iv} = vecang;
                    OP.phasemodSh{iv} = cell2mat(phasemodSh);

                    % real-mod shuff-pct
                    OP.modPctRank{iv} = 100*(1-(sum(OP.phasemodSh{iv} > OP.phasemod{iv})./numShuf));
                    fprintf('iv%d mod > %.02f pctShufs.\n', iv, modPctRank)
                end
            end
        end
        ic = ic +1;
        out(a).output{1}(ic) = OP;
    end
end
end
function out = init_out()
out.pval = [];
out.ILI = [];
out.spkOffsetAll = [];
out.spikePctSinceLick = [];
out.spikeLickPhase = [];
out.meanMRVmag = [];
out.vecang = [];
out.phasemod = [];
out.modPctRank = [];

out.area = '';
out.subarea = '';

out.meanMRVmagSh = [];
out.vecangSh = [];
out.phasemodSh = [];
end