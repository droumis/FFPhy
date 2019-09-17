
%{
- find relations between info channels

need to get the coherence/correlation of the different info sources..

- is my previous hypothesis about theta phase/power peri ripple still
valid?
- if the licking results turn out then i guess i'd need to look at ripples
away from the wells.. which i guess is what i had already done.. but now
the distinction is less clean with the licking thing..
- find examples of ripples away from the well.. what do they show?

observationA: licking appears to pace lfp and swrs
    - how to do spike-field pan spectrum coherence?
    - how to do spike-spike coherence?
    spikes:
        - licking (? i think so.. )
        - swrs 
    field: 
        - theta lfp 
    - is licking coherent with theta lfp? (spike-field coherence)
        - if yes, does the first licks reset the theta lfp or does it occur predictable at a
            certainly cycle of it.. both alternatives of this are interesting,
    - are the swrs coherent with theta lfp? (spike-field coherence)
    - is licking coherent with swrs (spike-spike coherence)

observationB: 


observationC:

%}

create_filter =1;
run_ff = 0;
save_ffdata = 0;
load_ffdata = 0;
processData = 0;
plotfigs = 0;
savefigs = 0;
pausefigs = 1;
%%
pconf = paramconfig;
Fp.animals = {'D10'}; %, 'JZ1', 'D12', 'D13', 'JZ1', 'JZ2', 'JZ3', 'JZ4'};
Fp.filtfunction = 'dfa_plotDataChunks';
Fp.add_params = {'wtrackdays','excludePriorFirstWell','<4cm/s', ...
    'excludeAfterLastWell', 'valid_ntrodes', 'exemplar_wepochs'};
Fp = load_filter_params(Fp, 'add_params', Fp.add_params);
animal = Fp.animals{1};
andef = animaldef(animal);
%%
if create_filter
    F = createfilter('animal', Fp.animals, 'epochs', Fp.epochfilter,  ...
       'excludetime', Fp.timefilter,'iterator',Fp.iterator,'eegtetrodes',Fp.tetfilter);
    F = setfilterfunction(F, Fp.filtfunction, Fp.datatypes, Fp.options{:});
end
day = F.epochs{1}(1,1);
epoch = F.epochs{1}(1,2);
%% load data
dio = loaddatastruct(andef{2}, animal, 'DIO', day);
task = loaddatastruct(andef{2}, animal, 'task', day);
events = loaddatastruct(andef{2}, animal, 'ca1rippleskons', day);

%% get licks
lickDIOIdx = task{day}{epoch}.inputdio;
rewardIdx = task{day}{epoch}.outputdio + length(dio{day}{epoch})/2;

isinput = cellfun(@(x) isequal(x.input,1), dio{day}{epoch}, 'un', 1);
dioID = cellfun(@(x) str2double(regexp(x.original_id,'\d*','Match')), dio{day}{epoch}, 'un', 1);
licks = sortrows(cell2mat(cellfun(@(x) [x.times repmat(str2double(regexp(x.original_id,'\d*','Match')), ...
    length(x.times),1)], dio{day}{epoch}(ismember(dioID(isinput), lickDIOIdx)), 'un', 0)'),1);
scatter(licks(:,1),ones(length(licks(:,1)),1),10,licks(:,2), 'marker', '+')
%% filter swrs
excludeIntervals = F.excludetime{1}{ismember(F.epochs{1},[day epoch], 'rows')};
eventTime = [events{day}{epoch}{1}.starttime events{day}{epoch}{1}.endtime];
eventTime = eventTime(~isExcluded(eventTime(:,1),excludeIntervals),:);
swrStart = eventTime(:,1);
%% get ISI's
subplot(3,2,1)
lickISI = diff(licks(:,1));
lickISI = lickISI(lickISI<1);
bins = (0:.01:1);
h = hist(lickISI, bins);
bar(bins,h); hold on;
set(gca, 'XScale', 'linear', 'YScale','linear')
[pks, locs] = findpeaks(h,'SortStr', 'descend', 'MinPeakProminence',50);
% plot(bins(locs), pks, 'or')
text(bins(locs(1)), pks(1), {sprintf('%.0f Hz',1/bins(locs(1)))})
text(bins(locs(2)), pks(2), {sprintf('%.0f Hz',1/bins(locs(2)))})
hold off
axis tight
xlim([0 .5])
xlabel('ILI s')
ylabel('lick counts')
colormap(hot)
%% get autocorr
subplot(3,2,2)
f = 1000;
bin = .01;
tmax = .5;
% axc = xcorr(lickInd(:,1) - mean(lickInd(:,1)), f, 'coeff');
xcs = spikexcorr(licks(:,1),licks(:,1), bin, tmax);
[~,i] = max(xcs.c1vsc2);
xcs.c1vsc2(i) = 0;
bar(xcs.time, xcs.c1vsc2)
% set(gca, 'yscale', 'log');
% n = length(lickInd(:,1));
% line([-f f], [2/sqrt(n) 2/sqrt(n)], 'color', 'r')
% line([-f f], -[2/sqrt(n) 2/sqrt(n)], 'color', 'r')
axis tight
% xlim([-200 200]);
ylabel('acorr')
xlabel('lag s')
%% xcorr
% subplot(3,2,3)
bin = .01;
starttime = min([swrStart; licks(:,1)]);
endtime = max([swrStart; licks(:,1)]);
time = [starttime:bin:endtime]';
 % hist + matlab's xcorr
% swrhist = hist(swrStart, time);
% lickhist = hist(licks(:,1), time);
% % %% lick-swr xcorr
% 
% tmax = .5;
% f = tmax/bin; % 
% axc = xcorr(swrhist, lickhist, f, 'coeff');
% [~,i] = max(axc);
% axc(i) = 0;
% bar(-f*bin:bin:f*bin, axc)
% xticks([-0.5:.5:.5])
% axis tight
% ylabel('xcorr')
% xlabel('lag s')
% hold off
% lab's spiketimes xcorr
subplot(3,2,3)
xcs = spikexcorr(swrStart,licks(:,1), bin, tmax);
bar(xcs.time, xcs.c1vsc2); axis tight; xlim([-.51 .51]); xticks([-0.5:.5:.5])
ylabel('spikexcorr 10ms bin')
xlabel('lag s')
hold off

%% xcorr
bin = .025;
% starttime = min([swrStart; licks(:,1)]);
% endtime = max([swrStart; licks(:,1)]);
% time = [starttime:bin:endtime]';
%  % hist + matlab's xcorr
% swrhist = hist(swrStart, time);
% lickhist = hist(licks(:,1), time);
% % %% lick-swr xcorr
% subplot(3,2,5)
% tmax = .5;
% f = tmax/bin; % 
% axc = xcorr(swrhist, lickhist, f, 'coeff');
% [~,i] = max(axc);
% axc(i) = 0;
% bar(-f*bin:bin:f*bin, axc)
% xticks([-0.5:.5:.5])
% axis tight
% ylabel('xcorr')
% xlabel('lag s')
% hold off
% lab's spiketimes xcorr\

subplot(3,2,4)
xcs = spikexcorr(swrStart,licks(:,1), bin, tmax);
bar(xcs.time, xcs.c1vsc2); axis tight; xlim([-.51 .51]); xticks([-0.5:.5:.5])
ylabel('spikexcorr 25ms bin')
xlabel('lag s')
hold off

%% %% get random xcorr between swr's and licks

subplot(3,2,5)
bin = .01;
xcs = spikexcorr(swrStartShift,licks(:,1), bin, tmax);
bar(xcs.time, xcs.c1vsc2); axis tight; xlim([-.51 .51]); xticks([-0.5:.5:.5])
title('shuffled')
ylabel('spikexcorr 10ms bin')
xlabel('lag s')
%%

figure('units','normalized', 'position', [.1 .1 .3 .3], 'color', 'white')
bin = .02;
tmax = 1;
realxcs = spikexcorr(swrStart,licks(:,1), bin, tmax);
bar(realxcs.time, realxcs.c1vsc2, 'k'); axis tight; hold on;
% shuffle swrs
x = [];
r = randi([-200 200],length(swrStart),1000)/1e3;
for i = 1:1000
    swrStartShift = sort(swrStart+r(:,i));
    xcs = spikexcorr(swrStartShift,licks(:,1), bin, tmax);
    x(i,:) = xcs.c1vsc2;
end
xstd = std(x);
xmean = mean(x);
plot(xcs.time, xmean, 'color', 'r', 'linewidth', 1);
fill([xcs.time'; flipud(xcs.time')],[xmean'-xstd';flipud(xmean'+xstd')],'b',...
    'linestyle','none', 'facealpha', .2);
ylim([min(realxcs.c1vsc2)-10 max(realxcs.c1vsc2)+10])
hold off
ylabel('xcorr')
xlabel('lag s')
titl = 'lickVSswr bin20ms shuff200msSTD';
title(titl)
save_figure(sprintf('%s/lick_swr_xcorr/', pconf.andef{4}), titl)

%%
bar(xcs.time, xcs.c1vsc2); axis tight; xlim([-.51 .51]); xticks([-0.5:.5:.5])
title('shuffled')
ylabel('spikexcorr 25ms bin')
xlabel('lag s')

%%% swrIdx = knnsearch(time, swrStart);
% swrInd = zeros(length(time),1);
% swrInd(swrIdx) = 1;
% 
% lickIdx = knnsearch(time, licks(:,1));
% lickInd = zeros(length(time),1);
% lickInd(lickIdx) = 1;
% plot(cumsum(lickInd)); hold on;
% plot(cumsum(swrInd))

% imagesc(time, [0 1], [lickInd sumswrInd]');
% colormap(1-gray)
% scatter(licks(:,1),zeros(length(licks(:,1)),1),10,licks(:,2), 'marker', '+');
%% spectrogram?
% movingwin = [.5 .05];
% params.fpass = [0 50];
% params.tapers = [2 3];
% [S, T, F] = mtspecgrampb(lickInd

%% filter frameworkifying

create_filter = 1;
run_ff = 1;
save_ffdata = 0;
load_ffdata = 0;
%% filter params
pconf = paramconfig;
Fp.animals = {'D10','D12','D13', 'JZ1', 'JZ2', 'JZ4'};
Fp.filtfunction = 'dfa_lickswrcorr';
Fp.add_params = {'wtrackdays','excludePriorFirstWell','<4cm/s', ...
    'excludeAfterLastWell'};%, 'exemplar_wepochs'}; %'excludeNoise',
Fp = load_filter_params(Fp, 'add_params', Fp.add_params);
%% FF
if create_filter
    F = createfilter('animal', Fp.animals, 'epochs', Fp.epochfilter,  ...
        'excludetime', Fp.timefilter,'iterator',Fp.iterator,'eegtetrodes',Fp.tetfilter);
    F = setfilterfunction(F, Fp.filtfunction, Fp.datatypes, Fp.options{:}); end
if run_ff; F = runfilter(F); for d = 1:length(F); F(d).datafilter_params = Fp; end; end
if save_ffdata
    save_data(F, Fp.paths.filtOutputDirectory, Fp.paths.filenamesave, 'filetail',...
        sprintf('_%s', Fp.epochEnvironment)); end
if load_ffdata
    F = load_data(Fp.paths.filtOutputDirectory, Fp.paths.filenamesave, ...
        Fp.animals, 'filetail', sprintf('_%s', Fp.epochEnvironment)); end



















