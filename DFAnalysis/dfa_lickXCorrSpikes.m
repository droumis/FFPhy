

function out = dfa_lickXCorrSpikes(idx, excludeIntervals, varargin)

%{
run with filterframework via singleDayCellAnal
- this should really be done after epochs are combined... 
- gets spikes (burst) x lick (burst) xcorr and phase clustering
- example script: lickXcorrSU_20190916.m

idx: [day ntrode cluster eps:]
excludeIntervals: [start end; ...] timefilter (applied to spiketimes)
varargin required data: {'spikes', 'lick'}

- phase modulation score as used in karalis sirota 2018/19
- uses the circ stats toolbox's circ_rtest to get the pval and z which is (R^2/n)..
log(R^2/n) where R is the normed mean resultant vector and n is num samples
circ_r gets you the mean resultant vector, r
- note, karalis and sirota don't use any point process with less than 200 points in
the times of interest..
$DR19
%}

fprintf('%d %d %d %d %d\n', idx)
reqData = {'lick', 'spikes', 'cellinfo'};
for s = 1:length(reqData)
    if ~any(cellfun(@(x) strcmp(x,reqData{s}), varargin(1:2:end), 'un', 1))
        error('missing required data')
    end
end


bin = .01;
tmax = 1;
excShortBin = bin*2;
excLongBin = .250;
minILIthresh = .06; % seconds
maxILIthresh = .250; % seconds
minBoutLicks = 3;
rmsmincounts = 1; % min bin count within rmstamax. otherwise nan
rmstmax = .25; % seconds
computeShuf = 1;
numShufs = 100;
maxShift = 500; % ms

% eventType = 'lick';
% psthbin = .001;
% TF = 1; % legacy tetfilter
% plotfigs = 1;
% displayplots = 0;
% saveplots = 1;
% savefigas = {'png', 'mfig'};
if ~isempty(varargin)
    assign(varargin{:});
end

% init output
out = init_out();
out.index = idx;
out.bin = bin;
out.tmax = tmax;

day = idx(1);
nt = idx(2);
clust = idx(3);
eps = idx(4:5);
%% spiketimes
spikeTimes = [];
for e = eps
    try
        spikeTimes = [spikeTimes; spikes{day}{e}{nt}{clust}.data(:,1)];
    catch
        continue
    end
end

if isempty(spikeTimes)
    fprintf('spiketimes empty\n');
    return
end

includetimes = ~isExcluded(spikeTimes, excludeIntervals);
spikeTimes = spikeTimes(includetimes);
if isempty(spikeTimes)
    fprintf('spikeimes empty\n');
    return
end

%% Get lick-burst intervals
licksVecs = getLickBout([], animal, [day eps(1); day eps(2)], 'lick', lick, 'maxIntraBurstILI', ...
    maxILIthresh, 'minBoutLicks', minBoutLicks);
burstIntvs = [];
for e = eps
try
    burstIntvs = [burstIntvs; vec2list(licksVecs{day}{e}.lickBout, licksVecs{day}{e}.time)];
catch
    fprintf('error defining lick bouts for %d %d \n', day, ep)
    return
end
end
%% get licks within lick-burst intervals
lbLicks = [];
idlicks = [];
for e = eps
    licks = lick{day}{e}.starttime;
    idlicks = [idlicks; lick{day}{e}.id];
    lbLicks = [lbLicks; licks(logical(isExcluded(licks, burstIntvs)))]; % isIncluded
end
fprintf('%d licks within %d bursts \n',numel(lbLicks), size(burstIntvs,1))
%% time mod
%% run spike x lick xcorr
time = (-tmax+bin/2):bin:(tmax-bin/2);
xc = spikexcorr(spikeTimes, lbLicks, bin, tmax);
normxc = xc.c1vsc2 ./ sqrt(xc.nspikes1 * xc.nspikes2); % normalize xc
nstd = round(excShortBin/bin);
g1 = gaussian(nstd, 2*nstd+1);
smthxc = smoothvect(normxc, g1); % smooth xc
excorr = nanmean(excesscorr(xc.time, xc.c1vsc2, xc.nspikes1, xc.nspikes2, excShortBin, ...
    excLongBin));
% T = abs(diff(swr{day}{ep}{1}.timerange)); %total ep time
% xcrms = xcorrrms(xc.time, xc.c1vsc2, rmstmax, rmsmincounts); % RMS
% p1 = xc.nspikes1/T;
% p2 = xc.nspikes2/T; % fr in Hz
% expProb = p1*p2; % per sec. Expected probability

out.time = time;
out.xc = xc;
out.normxc = normxc;
out.smthxc = smthxc;
out.excorr = excorr;
% out.xcrms = xcrms;
% out.p1 = p1;
% out.p2 = p2;
% out.expProb = expProb;

%% phasemod stuff
%% get spike events within lick-burst intervals
burstSpikeTime = spikeTimes(logical(isExcluded(spikeTimes, burstIntvs)));
burstSpikeTime = sort(burstSpikeTime(~isnan(burstSpikeTime)));
fprintf('spikes within %d lickbouts:  %d (%.02f pct swrs) \n', size(burstIntvs,1), ...
    numel(burstSpikeTime), length(burstSpikeTime)/length(spikeTimes)*100)
if isempty(burstSpikeTime)
    fprintf('no spikes in lick bouts for %d %d %d %d %d skipping\n', idx)
    return
end

%% Get the containing licks for each swr
% histc given licks as edges: get bin idx. lickedges(binidx) is the prior edge of bin 
[~,~,spikeLickBinidx] = histcounts(burstSpikeTime, lbLicks);
% exclude swr's in 0 bin, before first lick edge
spikeLickBinidx = spikeLickBinidx(spikeLickBinidx > 0); 
burstSpikeTime = burstSpikeTime(spikeLickBinidx > 0);
% time since pre-containing lick
spikeTimeSinceLick = burstSpikeTime - lbLicks(spikeLickBinidx);
% get inter lick interval containing each spike
spikeBinILI = lbLicks(spikeLickBinidx+1) - lbLicks(spikeLickBinidx); % +1 is post-containing lick
% exclude swr-containing ili out of lick burst range 60 : 250 ms
spikeValidILI = all([spikeBinILI > minILIthresh spikeBinILI < maxILIthresh],2);
fprintf('lBurst spikes in valid ILI: %d (%.02f pct) \n', sum(spikeValidILI), ...
    sum(spikeValidILI)/length(spikeTimes)*100)
burstSpikeTime = burstSpikeTime(spikeValidILI);
if isempty(burstSpikeTime)
    return
end
spikeBinILI = spikeBinILI(spikeValidILI);
spikeTimeSinceLick = spikeTimeSinceLick(spikeValidILI);
% save spike-containing burst id
[~,~,spikeLickBurstIdx] = histcounts(burstSpikeTime, [burstIntvs(:,1); inf]);
burstContainSpike = burstIntvs(spikeLickBurstIdx,:);

%% phasemod
spikePctSinceLick = spikeTimeSinceLick ./ spikeBinILI;
spikeLickPhase = 2*pi*spikePctSinceLick; % pct of full cycle
meanvec = mean(exp(1i*spikeLickPhase)); % get mean resultant vector
meanMRVmag = abs(meanvec); % vector magnitude
vecang = angle(meanvec); % for plotting mrv
[~, z] = circ_rtest(spikeLickPhase); % z is mean res vec
phasemod = log(z); % log variance normalizes (karalis,sirota)

% out.spikeTimes = spikeTimes;
% out.licks = licks;
% out.idlicks = idlicks;
% out.boutIntvs = burstIntvs;
% out.swrTimeSinceLick = spikeTimeSinceLick;
% out.swrBinILI = spikeBinILI;
% out.burstSWRStart = burstSpikeTime;
out.burstContainSpike = burstContainSpike;
% % out.dayEpoch = repmat([day ep], length(burstSpikeTime),1);
out.spikePctSinceLick = spikePctSinceLick;
out.spikeLickPhase = spikeLickPhase;
out.meanMRVmag = meanMRVmag;
out.vecang = vecang;
out.phasemod = phasemod;


%% Compute shuffled lickDin x swr
if computeShuf
    tic
    for i = 1:numShufs
        % shuffle spike times locally 
        r = randi([-maxShift maxShift],length(spikeTimes(:,1)),1)/1e3;
        spikeTimesShuf = sort(spikeTimes+r);
       %% time mod
       %% run spike x lick xcorr
        xc = spikexcorr(spikeTimesShuf, lbLicks, bin, tmax);
        normxc = xc.c1vsc2 ./ sqrt(xc.nspikes1 * xc.nspikes2); % normalize xc
        try
            smthxc = smoothvect(normxc, g1);
        catch % if normxc is empty
            continue
        end
        excorr = nanmean(excesscorr(xc.time, xc.c1vsc2, xc.nspikes1, xc.nspikes2, excShortBin, ...
            excLongBin));
%         xcrms = xcorrrms(xc.time, xc.c1vsc2, rmstmax, rmsmincounts); % RMS
%         p1 = xc.nspikes1/T;
%         p2 = xc.nspikes2/T; % fr in Hz
%         expProb = p1*p2; % per sec. Expected probability
        
        out.xcShuf{end+1} = xc;
        out.normxcShuf{end+1} = normxc;
        out.smthxcShuf{end+1} = smthxc;
        out.excorrShuf{end+1} = excorr;
%         out.xcrmsShuf{end+1} = xcrms;
%         out.p1Shuf{end+1} = p1;
%         out.p2Shuf{end+1} = p2;
%         out.expProbShuf{end+1} = expProb;

       %% phasemod stuff
       %% get spike events within lick-burst intervals
        burstSpikeTime = spikeTimesShuf(logical(isExcluded(spikeTimesShuf, burstIntvs)));
        burstSpikeTime = sort(burstSpikeTime(~isnan(burstSpikeTime)));
%         fprintf('spikes within %d lickbouts:  %d (%.02f pct swrs) \n', size(burstIntvs,1), ...
%             numel(burstSpikeTime), length(burstSpikeTime)/length(spikeTimesShuf)*100)
        if isempty(burstSpikeTime)
            continue
        end
        
        % histc given licks as edges: get bin idx. lickedges(binidx) is the prior edge of bin 
        [~,~,spikeLickBinidx] = histcounts(burstSpikeTime, lbLicks);
        % exclude spikes in 0 bin, before first lick edge
        spikeLickBinidx = spikeLickBinidx(spikeLickBinidx > 0);
        burstSpikeTime = burstSpikeTime(spikeLickBinidx > 0);
        % time since pre-containing lick
        spikeTimeSinceLick = burstSpikeTime - lbLicks(spikeLickBinidx);
        % get inter lick interval containing each spike
        spikeBinILI = lbLicks(spikeLickBinidx+1) - lbLicks(spikeLickBinidx);
        % exclude spike-containing ili out of lick burst range 60 : 250 ms
        spikeValidILI = all([spikeBinILI > minILIthresh spikeBinILI < maxILIthresh],2);
        burstSpikeTime = burstSpikeTime(spikeValidILI);
        spikeBinILI = spikeBinILI(spikeValidILI);
        spikeTimeSinceLick = spikeTimeSinceLick(spikeValidILI);
%         [~,~,spikeLickBurstIdx] = histcounts(burstSpikeTime, [burstIntvs(:,1); inf]);
%         burstContainSpike = burstIntvs(spikeLickBurstIdx,:);
        if isempty(spikeTimeSinceLick)
            return
        end
       % phasemod
        spikePctSinceLick = spikeTimeSinceLick ./ spikeBinILI;
        spikeLickPhase = 2*pi*(spikePctSinceLick);
        meanvec = mean(exp(1i*spikeLickPhase));
        meanMRVmag = abs(meanvec);
        vecang = angle(meanvec);
        [~, z] = circ_rtest(spikeLickPhase); % pval is the stat rayleigh test. z is mean res vec
        phasemod = log(z); % i think log makes it 'variance normalized' (karalis,sirota)
        
%         out.burstContainSpikeShuf{end+1} = burstContainSpike;
%         out.burstSpikeTimeShuf{end+1} = burstSpikeTime;
%         out.swrTimeSinceLickShuf{end+1} = spikeTimeSinceLick;
%         out.spikePctSinceLickShuf{end+1} = spikePctSinceLick;
        out.spikeLickPhaseShuf{end+1} = spikeLickPhase;
        out.meanMRVmagShuf{end+1} = meanMRVmag;
        out.vecangShuf{end+1} = vecang;
        out.phasemodShuf{end+1} = phasemod;
%         
%         out.swrBinILIShuf{end+1} = spikeBinILI;
%         out.swrInBurstStartShuf{end+1} = burstSpikeTime;
%         out.swrTimeSinceLickShuf{end+1} = spikeTimeSinceLick;
%         out.swrBurstIntervalShuf{end+1} = burstContainSwr;
        
    end
    fprintf('shuffle took %.02f s\n', toc);
end

%%


% 
% 
% %% Get the containing licks for each spike
% [~,~,swrLickBinidx] = histcounts(spikeTimes, lbLicks);
% % exclude swr's in 0 bin. when would this happen? maybe if same time as first lick
% swrLickBinidx = swrLickBinidx(swrLickBinidx > 0);
% swrInBurstStart = swrInBurstStart(swrLickBinidx > 0);

end

function out = init_out()
out.index = [];
out.animal = [];
out.time = [];
% out.licks = [];
% out.idlicks = [];
out.boutIntvs = [];
% out.swrStart = [];
% out.dayEpoch = [];
% xcTime
out.xc = [];
out.normxc = [];
out.smthxc = [];
out.excorr = [];
% out.xcrms = [];
% out.p1 = [];
% out.p2 = [];
% out.expProb = [];
% iliPhase
% out.swrBinILI = [];
% out.swrBurstInterval = [];
% out.swrInBurstStart = [];
% out.swrTimeSinceLick = [];
out.meanMRVmag = [];
% out.swrLickPhase = [];
% out.swrPctSinceLick = [];
out.vecang = [];
out.phasemod = [];

%% shuffle
% out.swrStartShuf = [];
% xcTime
out.xcShuf = {};
out.normxcShuf = {};
out.smthxcShuf = {};
out.excorrShuf = {};
% out.xcrmsShuf = {};
% out.p1Shuf = {};
% out.p2Shuf = {};
% out.expProbShuf = {};
% iliPhase
% out.swrBinILIShuf = {};
% out.swrBurstIntervalShuf = {};
% out.swrInBurstStartShuf = {};
% out.swrTimeSinceLickShuf = {};
out.meanMRVmagShuf = {};
out.spikeLickPhaseShuf = {};
% out.swrPctSinceLickShuf = {};
out.vecangShuf = {};
out.phasemodShuf = {};
end

% %% recover lick times from excludeIntervals ? i should just pass it in no?
% evvar = varargin{find(cellfun(@(x) ~isempty(x), strfind(varargin(1:2:end), ...
%     eventType), 'un', 1))*2-1};
% lick = eval(evvar);
% try
%     lickTimes = lick{day}{epoch}.starttime;
%     id = lick{day}{epoch}.id; % to distinguish events at diff reward wells
% catch
%     lickTimes = lick{day}{epoch}{TF}.starttime;
%     id = lick{day}{epoch}{TF}.id;
% end
% 
% %% timefilter lick events
% includetimes = ~isExcluded(lickTimes, excludeIntervals); % result: 1 include, 0 exclude
% lickTimes = lickTimes(includetimes);
% if isempty(lickTimes)
%     fprintf('evenTimes empty\n');
%     return
% end
% 
% out.eventTimes = lickTimes;

%% move plotting to a script
% if plotfigs
%     %% plot
%     Pp=load_plotting_params({'defaults','dfa_lickXCorrSpikes'});
%     if saveplots && ~displayplots; close all;
%         ifig = figure('Visible','off','units','normalized','position', Pp.position, ...
%             'color','white', 'InvertHardcopy', 'off');
%     else
%         ifig = figure('units','normalized','position',Pp.position, 'color','white', ...
%             'InvertHardcopy', 'off');
%     end
%     
%     nrow = 3;
%     ncol = 1;
%     %% lick psth
%     psthtime = (-tmax-0.5*psthbin):psthbin:(tmax+0.5*psthbin);
%     if ~isempty(eventTimes)
%         psth = nan(length(eventTimes),length(psthtime));
%         for r=1:length(eventTimes)
%             shist = histc(spiketimes, eventTimes(r) + psthtime);
%             psth(r,:) = logical(shist);
%         end
%     end
%     sf1 = subaxis(nrow,ncol,1,'SpacingVert', Pp.SpVt, 'SpacingHoriz', Pp.SpHz, ...
%         'MarginLeft', Pp.MgLt, 'MarginRight', Pp.MgRt, 'MarginTop', Pp.MgTp, ...
%         'MarginBottom', Pp.MgBm); set(gca, 'Tag', 'ripkons');
%     h=zoom; h.Motion = 'horizontal'; h.Enable = 'on'; p=pan; p.Motion='horizontal';
%     [xx, yy] = find(psth');
%     uclr = cell2mat(cellfun(@(x) repmat(x(1),x(2),1),num2cell([id sum(psth,2)],2), ...
%         'un',0));
%     f = scatter(xx/1000-1.001,yy,Pp.psthSize, uclr, 'filled');
%     colormap(cbrewer('qual', 'Accent', 3))
%     set(gca, 'color', 'k')
%     sf1.GridColor = 'w';
%     f.Marker = 'd';
%     f.MarkerEdgeAlpha = 0;
%     f.MarkerFaceAlpha = 0.4;
%     set(gca, 'TickDir','out', 'YGrid', 'off', 'XGrid', 'on', 'TickLength', [0.001 0], ...
%         'Tag', 'ndata', 'Xticklabel',[]);
%     axis tight;
%     ylabel('licknum','FontSize',12,'FontWeight','bold', 'FontName','Arial')
%     axis tight
%     line([0 0],ylim, 'color','red', 'linewidth',2, 'Color', [1 0 0 .5]);
%     %% psth pdf
%     sf2 = subaxis(nrow,ncol,2);
%     h = histogram(xx/1000-1.001, psthtime(1:20:end), 'facecolor', 'k', 'Normalization', ...
%         'pdf'); 
%     h.EdgeColor = 'none';
%     axis tight
%     line([0 0],ylim, 'color','red', 'linewidth',2, 'Color', [1 0 0 .5]);
%     ylabel('pdf','FontSize',12,'FontWeight','bold', 'FontName','Arial')
%     set(gca, 'YGrid', 'off', 'XGrid', 'on','TickDir','out', 'TickLength', [0.001 0]);
%     xticklabels([])
%     hold off;
%     %% lick XCORR spikes
%     sf3 = subaxis(nrow,ncol,3);
%     realxcs = spikexcorr(spiketimes,eventTimes, bin, tmax);
%     plot(realxcs.time, realxcs.c1vsc2', 'k', 'linewidth', 2)%, 'facecolor', 'k')%, '-k', 'filled', 'markersize', 5); 
%     hold on;
%     axis tight;
%     plot(xcs.time, xmean, 'color', 'b', 'linewidth', 2); hold on;
%     fill([xcs.time'; flipud(xcs.time')],[xmean'-xstd';flipud(xmean'+xstd')],'b',...
%         'linestyle','none', 'facealpha', .1);
%     line([0 0],ylim, 'color','red', 'linewidth',2, 'Color', [1 0 0 .5]);
%     hold off
%     ylabel('xcorr')
%     xlabel('lag time s')
%     set(gca, 'YGrid', 'off', 'XGrid', 'on','TickDir','out', 'TickLength', [0.001 0]);
% 
%     %% ----- link x axis -----
%     allAxesInFigure = findall(ifig,'type','axes'); linkaxes(allAxesInFigure, 'x');
%     % ---- super axis -----
%     sprtit = sprintf('%s %d %d %d %d lickXCorrSpikes bin%.0f shuff%.0f',animal,day,epoch,...
%         ntrode, clust, bin*1e3, shuff);
%     sprax = axes('Position',[0 0 1 1],'Visible','off', 'Parent', ifig);
%     iStitle = text(.5, .98, {sprtit}, 'Parent', sprax, 'Units', 'normalized');
%     set(iStitle,'FontWeight','bold','FontName','Arial','horizontalAlignment','center', ...
%         'FontSize',12);
%     h = get(gcf,'Children');
%     set(gcf,'Children',flip(h)); % super axis to bottom. allows for zoom/pan
%     % ---- pause, save figs ----
%     if displayplots; pause; end
%     if saveplots
%         [~, fname,~] = fileparts(mfilename('fullpath'));
%         outdir = sprintf('%s/%s/%s/', pconf.andef{4},fname,animal);
%         save_figure(outdir, sprtit, 'savefigas', savefigas); 
%     end
% end