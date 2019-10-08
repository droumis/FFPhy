
% phase modulation score as used in karalis sirota 2018/19
% use the circ stats toolbox's circ_rtest to get the pval and z which is (R^2/n).. 
% log(R^2/n) where R is the normed mean resultant vector and n is num samples
% circ_r gets you the mean resultant vector, r
% right now i'm basically treating swrs as a cluster.. in that the
% procedure should be very similar to quantify the phase modulation with each cluster

% so i can get 1 phase mod score per time period: epoch? or day? or animal?
% karalis, sirota don't use any point process with less than 200 points in
% the times of interest.. which probably means 1 measure per animal (or
% maybe some days there's enough swr's during the lick bouts?).. the thing
% that sucks is the variability of the phase of each beam break has some
% variability probably..

function out = dfa_lickswrcorr(idx, excludeIntervals, varargin)
pconf = paramconfig;
pausefigs = 0;
savefigs = 1;
savefigas = 'png';
plotfigs = 1;
bin = .01;
tmax = 1;
boutNum = 10;
lickgap = .500;
sw1 = bin*2;
% sw1 = 0.005; 
sw2 = .250;
rmstmax = 0.1;
rmsmincounts = 50;
numshuffs = 1000;
shuffOffset = 200; %ms
returnoutput = 1;
sigperc = .975;
if ~isempty(varargin); assign(varargin{:}); end
day = idx(1,1); epoch = idx(1,2);
[~, fname,~] = fileparts(mfilename('fullpath'));
outdir = sprintf('%s/%s/%s/', pconf.andef{4},fname,animal,animal);

% get the ripplekons events, which could be variably named
evid = find(contains(varargin(1:2:end), 'rippleskon'));
o = [1:2:length(varargin)]+1;
events = varargin{o(evid)};

out = make_blank(idx);

try
    eventTime = [events{day}{epoch}{1}.starttime events{day}{epoch}{1}.endtime];
catch; fprintf('no events detected for day%d ep%d\n',day,epoch); return; end
eventTime = eventTime(~isExcluded(eventTime(:,1),excludeIntervals),:);
swrStart = eventTime(:,1);

% total time
T = diff(events{day}{epoch}{1}.timerange);

% licks, bouts
andef = animaldef(animal);
licks = lick{day}{epoch}.starttime;
idlicks = lick{day}{epoch}.id;
licksVecs = getLickBout(andef{2}, animal, [day epoch]);
try
    boutIntvs = vec2list(licksVecs{day}{epoch}.lickBout, licksVecs{day}{epoch}.time);
catch
    fprintf('error defining lick bouts for %d %d\n', day, epoch)
    return
end
swrinbouts = swrStart(logical(isExcluded(swrStart, boutIntvs)));
swrinbouts = sort(swrinbouts(~isnan(swrinbouts)));
fprintf('%d (%.02f perc) swrs within the %d lickbouts\n', numel(swrinbouts),...
    length(swrinbouts)/length(swrStart)*100, size(boutIntvs,1))
if isempty(swrinbouts)
    fprintf('no swrs in lick bouts for %d %d.. skipping\n', day, epoch)
    return
end
lickboutlicks = licks(logical(isExcluded(licks, boutIntvs)));
fprintf('%d (%.02f perc) licks within the %d lickbouts\n', numel(lickboutlicks),...
    length(lickboutlicks)/length(licks)*100, size(boutIntvs,1))

time = (-tmax+bin/2):bin:(tmax-bin/2);
nstd=round(sw1/(time(2) - time(1)));
g1 = gaussian(nstd, 2*nstd+1);

%% shuffles
meanMRVmagShuff = [];
smthxcShuff = [];
excorrShuff = [];
phasemodShuff = [];
r = randi([-shuffOffset shuffOffset],length(swrinbouts),numshuffs)/1e3;
swrStartShift = bsxfun(@plus,swrinbouts,r); %, make shuffle array
for i = 1:numshuffs
    % (mean/std xcorr norm smooth)
%     swrStartShift = sort(swrinbouts+r(:,i));
    xcSh = spikexcorr(sort(swrStartShift(:,i)), lickboutlicks, bin, tmax);
    normxc = xcSh.c1vsc2 ./ sqrt(xcSh.nspikes1 * xcSh.nspikes2); % normalize
    try 
        smthxcShuff(end+1,:) = smoothvect(normxc, g1);
    catch % if normxc is empty
        continue
    end
    % excorr
    excorrShuff(end+1) = nanmean(excesscorr(xcSh.time, xcSh.c1vsc2, xcSh.nspikes1, ...
        xcSh.nspikes2, sw1, sw2)); %0 lag
    % MRVmag
    [~,~,spbin] = histcounts(swrStartShift(:,i), lickboutlicks);
    % exclude swrs falling outside lickbouts
    goodswr = swrStartShift(spbin>0,i);
    spbin = spbin(spbin>0);
    lickbinstart = lickboutlicks(spbin);
    relswrtime = goodswr - lickbinstart; % swrtime - the nearest preceding lick in bout
    ili = lickboutlicks(spbin+1) - lickbinstart; % inter lick interval
    useILI = all([ili > .05, ili < .25],2); % include know lick bout intervals
%     fprintf('%d good swrs in lick bout intervals.. %d %d\n', sum(useILI), day, epoch)
    if sum(useILI) < 5
        fprintf('only %d swrs in lickbouts during a shuffle, skipping %d %d\n', ...
            sum(useILI), day, epoch)
        return
    end
    swrlickcylcephase = relswrtime(useILI) ./ ili(useILI);
    swrLickPhaseShuff = 2*pi*(swrlickcylcephase);
    meanvec = mean(exp(1i*swrLickPhaseShuff));
    meanMRVmagShuff(end+1) = abs(meanvec);
    vecangShuff = angle(meanvec);
    % phasemod
    [~, z] = circ_rtest(swrLickPhaseShuff); % pval is the stat rayleigh test. z is mean res vec
    phasemodShuff(end+1) = log(z); % i think log makes it 'variance normalized' (karalis,sirota)
end
smthxcShuffstd = std(smthxcShuff); % shuff xcorr std
smthxcShuffmean = mean(smthxcShuff); %shuff xcorr mean

%% real
% xcorr norm smooth
xc = spikexcorr(swrinbouts,lickboutlicks, bin, tmax);
normxc = xc.c1vsc2 ./ sqrt(xc.nspikes1 * xc.nspikes2); % normalize xc
smthxc = smoothvect(normxc, g1); % smooth xc
% excorr
excorr = nanmean(excesscorr(xc.time, xc.c1vsc2, xc.nspikes1, xc.nspikes2, sw1, sw2)); % 0 lag
% xcrms = xcorrrms(xc.time, xc.c1vsc2, rmstmax, rmsmincounts); % RMS
% p1 = xc.nspikes1/T;
% p2 = xc.nspikes2/T; % fr in Hz
% expProb = p1*p2; % per sec. Expected probability
% MRVmag
[~,~,spbin] = histcounts(swrinbouts, licks);
relswrtime = swrinbouts - licks(spbin); % swrtime - the nearest preceding lick in bout
fprintf('%.03f mean %.03f max %.03f min swr time since last lick \n', mean(relswrtime), ...
    max(relswrtime), min(relswrtime))
if any(relswrtime) < 0
    error('swr rel time must be pos')
end
ili = licks(spbin+1) - licks(spbin);
useILI = all([ili > .1, ili < .2],2);
if sum(useILI) < 20
    fprintf('only %d swrs in lickbouts, skipping %d %d\n', sum(useILI), day, epoch)
    return
end
fprintf('%.03f mean %.03f max %.03f min ILI \n', mean(ili), max(ili), min(ili))
swrlickcylcephase = relswrtime(useILI) ./ ili(useILI);
swrLickPhase = 2*pi*(swrlickcylcephase);
meanvec = mean(exp(1i*swrLickPhase));
meanMRVmag = abs(meanvec);
vecang = angle(meanvec);
% phasemod
[pval, z] = circ_rtest(swrLickPhase); % pval is the stat rayleigh test. z is mean res vec
phasemod = log(z); % i think log makes it 'variance normalized' (karalis,sirota)

%% plot
if plotfigs
    Pp=load_plotting_params({'defaults','dfa_lickswrcorr'}, 'pausefigs', pausefigs, ...
        'savefigs', savefigs);
    %% Xcorr norm smooth, shuff mean/std
    sf1 = subaxis(2,2,1,Pp.posparams{:});
    sf1.Tag = 'xcorr';
    % shuffled xcorr mean.std ghost trail
    plot(time, smthxcShuffmean, 'color', [0 0 1 .2], 'linewidth', 1);
    hold on;
    fill([time'; flipud(xc.time')],[smthxcShuffmean'-smthxcShuffstd';flipud(smthxcShuffmean'+smthxcShuffstd')],'b',...
        'linestyle','none', 'facealpha', .s1);
    % xcorr norm
    bar(xc.time, normxc, 'k', 'facealpha', .2, 'edgealpha', 0)
    % xcorr norm smooth
    plot(xc.time, smthxc, 'k')
    line([0 0], ylim, 'color', 'k', 'linestyle', '--', 'linewidth', .5)
    
    ylabel('xcorr');
    xlabel('time from lick s');
%     title('xcorr 
    hold off;
    %% excorr over shuff excorr cdf distr
    % relative swr from last lick
    sf2 = subaxis(2,2,2);
    histogram(excorrShuff, 60,'Normalization','probability','edgealpha', 0, 'facecolor', 'k');
    excsort = sort(excorrShuff);
    idxsig = round(sigperc*length(excorrShuff));
    line([excsort(idxsig) excsort(idxsig)], ylim, 'color', [0 0 0 .8], 'linestyle', '--');
    hold on;
    line([excorr excorr], ylim, 'color', 'r');
    excp = 1-sum(excorr>excorrShuff)/length(excorrShuff);
    title(sprintf('excorr %.03f p%.03f', excorr, excp));
    ylabel('probability')
    xlabel('excess corr')
    axis tight
    hold off;
%     %% mag distr
%     sf3 = subaxis(3,2,3);
%     histogram(meanMRVmagShuff, 60, 'Normalization', 'probability', 'edgealpha', 0, 'facecolor', 'k');
%     hold on;
%     mrvsort = sort(meanMRVmagShuff);
%     idxsig = round(sigperc*length(meanMRVmagShuff));
%     line([mrvsort(idxsig) mrvsort(idxsig)], ylim, 'color', [0 0 0 .8], 'linestyle', '--');
%     hold on
%     line([meanMRVmag meanMRVmag], ylim, 'color', 'r');
%     mrvp = 1-sum(meanMRVmag>meanMRVmagShuff)/length(meanMRVmagShuff);
%     title(sprintf('MRVmag %.03f p%.03f', meanMRVmag, mrvp));
%     ylabel('probability')
%     axis tight
%     hold off
    
    %% phase mod
    sf4 = subaxis(2,2,4);
    histogram(phasemodShuff, 100, 'Normalization', 'probability', 'edgealpha', 0, 'facecolor', 'k');
    hold on;
    mrvsort = sort(phasemodShuff);
    idxsig = round(sigperc*length(phasemodShuff));
    line([mrvsort(idxsig) mrvsort(idxsig)], ylim, 'color', [0 0 0 .8], 'linestyle', '--');
    hold on
    line([phasemod phasemod], ylim, 'color', 'r');
    modp = 1-sum(phasemod>phasemodShuff)/length(phasemodShuff);
    title(sprintf('phasemod %.03f p%.03f Rpval%.03f', phasemod, modp, pval));
    ylabel('probability')
    xlabel('log(Rayleigh Z)')
    axis tight
    hold off
    
    %% polar distr phase clustering, swrLickPhase, meanMRVmag
    sf3 = subaxis(2,2,3);
%     a = polarhistogram(swrLickPhase, 16, 'Normalization', 'pdf', 'edgealpha', 0,...
%         'facealpha', .5);
    polarplot([zeros(size(swrLickPhase,1),1) swrLickPhase]', ...
       repmat([0 1],size(swrLickPhase,1),1)', 'color', [0 0 0 .4], 'linewidth', 4);
   hold on
    polarplot([0; vecang], [0; meanMRVmag], 'color', [1 0 .3], 'linewidth', 4)
    grid off
    rticks([])
    thetaticks([])
    title('swr lick-phase')
    hold off
    axis tight
    %% ---- super axis -----
    sprtit = sprintf('%s %d %d lickswrcorr', animal,day,epoch);
    sprax = axes('Position',[0 0 1 1],'Visible','off', 'Parent', gcf);
    iStitle = text(.5, .98, {sprtit}, 'Parent', sprax, 'Units', 'normalized');
    set(iStitle,'FontWeight','bold','FontName','Arial','horizontalAlignment','center', ...
        'FontSize',12);
    h = get(gcf,'Children');
    set(gcf,'Children',flip(h)); % put super axis at bottom of axis stack. allows for zoom
    % ---- pause, save figs ----
    if pausefigs; pause; end
    if savefigs; save_figure(outdir, sprtit, 'savefigas', savefigas); end
end
if returnoutput
    out.index = idx;
    out.time = time;
    out.swrStart = swrStart;
    out.licks = licks;
    out.idlicks = idlicks;
    out.boutIntvs = boutIntvs;
    out.xc = xc;
    out.normxc = normxc;
    out.smthxc = smthxc;
    out.excorr = excorr;
%     out.xcrms = xcrms;
    out.swrinbouts = swrinbouts;
    out.relswrtime = relswrtime;
    out.meanMRVmag = meanMRVmag;
    out.swrLickPhase = swrLickPhase;
    out.vecang = vecang;
    out.phasemod = phasemod;
    out.meanMRVmagShuff = meanMRVmagShuff;
    out.swrLickPhaseShuff = swrLickPhaseShuff;
    out.excorrShuff = excorrShuff;
    out.phasemodShuff = phasemodShuff;
    out.smthxcShuff = smthxcShuff;
end
end
function out = make_blank(idx)
    out.index = idx;
    out.time = [];
    out.swrStart = [];
    out.licks = [];
    out.idlicks = [];
    out.boutIntvs = [];
    out.xc = [];
    out.normxc = [];
    out.smthxc = [];
    out.excorr = [];
%     out.xcrms = [];
    out.swrinbouts = [];
    out.relswrtime = [];
    out.meanMRVmag = [];
    out.swrLickPhase = [];
    out.vecang = [];
    out.phasemod = [];
    out.meanMRVmagShuff = [];
    out.swrLickPhaseShuff = [];
    out.excorrShuff = [];
    out.phasemodShuff = [];
    out.smthxcShuff = [];
end
