


animals = {'D10'};
filtfunction = 'dfa_riptrigspiking';
env = 'wtrack';
me = animaldef('demetris');

loaddata = 1;

paths = make_paths(filtfunction, env);
%% load data
if loaddata
    data = load_filter_output(paths.filtOutputDirectory, paths.filenamesave, ...
        animals);
end
%% get keys,inds into data for this nt,day
data_keys = cell2mat({data{1}.F.output{1}.index}');
use_nt = 11;
use_day = 3;
nt_data_inds = find(all([data_keys(:,3)==11 ismember(data_keys(:,1), ...
    use_day, 'rows')],2));
nt_day_keys = data_keys(nt_data_inds,:);

%% get the MU/SU keys for this nt,day
ian = 1;
andef = animaldef(animals{ian});
cellinfo = loaddatastruct(andef{2}, andef{3}, 'cellinfo', 3);
mu_keys = evaluatefilter(cellinfo, 'isequal($tags, {''mua''})');
su_keys = evaluatefilter(cellinfo, '~isequal($tags, {''mua''})');
nt_mu_keys = nt_day_keys(find(ismember(nt_day_keys, mu_keys, 'rows')),:);
nt_su_keys = nt_day_keys(find(ismember(nt_day_keys, su_keys, 'rows')),:);
%%
% i previously wrote this to vert stack epochs in a day but i want to do
% the precession R2 per epoch..
psth = combine_riptrigspiking_perntrode(data);
%%
% but first just get ISI distribution.. per swr
nt_psth = psth{1}([psth{1}.ntrode]'==use_nt);
%%
% for row in the mu psth, collect diff between spike times
diffs = {};
for r = 1:length(nt_psth.mucluster{use_day}{1}(:,1))
    diffs{r} = diff(find(nt_psth.mucluster{use_day}{1}(r,:)));
end
diffs = cell2mat(diffs);
%%
histogram(diffs(all([diffs<200 diffs>100],2)),1000)
% hmm opk that's not going to work bc there are so many spikes.. try using
% an FR psth
%%
% histc bins, where center bin is centered time 0
winlen = length(nt_psth.mucluster{use_day}{1}(1,:));
numrips = length(nt_psth.mucluster{use_day}{1}(:,1));
window = [0 winlen];
binsize = 1;
frbinsize = 5;
time = (-window(1)-0.5*binsize):binsize:(window(2)+0.5*binsize);
frtime = (-window(1)-0.5*frbinsize):frbinsize:(window(2)+0.5*frbinsize);
%%
f = [];
for r = 1:numrips
    spikesinwin{r} = find(nt_psth.mucluster{use_day}{1}(r,:));
    %     ifr = instantfr(spikesinwin{r}, time);
    %
    bifr = instantfr(spikesinwin{r}, frtime);
    %     hbifr = histc(spikesinwin{r}, frtime);
    f(r,:) = ksdensity(spikesinwin{r},frtime,'Bandwidth',10,'Function','pdf');
end
%%
clf
imagesc(flipud(f(1:nt_psth.eplengths{use_day}(1),:)))
colormap('parula')
%%
subplot(1,2,1)
spk = nt_psth.mucluster{use_day}{1}(r,:);
imagesc(spk);%),[1], '.k')

subplot(1,2,2)
% for r = 1:numrips
Fs = 1000; %sampling
t = 0:1/Fs:1; %time
% f = 5;
% x = sin(2*pi*t*f);
x = spk; % signal
% nfft = 1024; %
nfft = 2^nextpow2(length(x));
X = fft(x,nfft);
X = X(1:nfft/2);
mx = abs(X);
f = (0:nfft/2-1)*Fs/nfft;
% % figure
% % plot(t,x)
figure(1)
plot(f,mx)
xlim([0 50])
ylim([0,80])
title('power spec')
xlabel('Hz')
ylabel('power')

%%  just use 7 to create a wavelet..
% then convolve the Kernel density estimation with this wavelet and get
% phase at swr onset

% how to create a wavelet?
% /home/droumis/Src/Matlab/filterframework_dr
% itpc.m ( mike cohens original matlab example of itpc)
% he uses exp(21i*pi*freq.*time) .* exp(-time.^2./(2*s^2));
% where s is numwavecycles/(2*pi*freq)
% then does fft of the wavelet and the data (just using matlabs fft)
% then multiples the waveletFFT with dataFFT
% then takes inverse fft

%% HERE WE GO dfs_phasecoherence
% /home/droumis/Src/Matlab/filterframework_dr/DFScripts/dfs_phasecoherence.m
% that makes the AS with the first below, then computes IXPC with the other
% -/home/droumis/Src/Matlab/filterframework_dr/Functions/computeAnalyticSignal.m
% -/home/droumis/Src/Matlab/filterframework_dr/Functions/computeIXPC.m

%% try just usiong regular matlab tools to get phase

Fs = 1000;
t = linspace(0,1,Fs);
x = cos(2*pi*100*t)+0.5*randn(size(t));

fc = 150;
Wn = (2/Fs)*fc;
b = fir1(20,Wn,'low',kaiser(21,3));
fvtool(b,1,'Fs',Fs)
%%
animal = animals{ian};
pxxs = [];
for r = 1:numrips
    spk = nt_psth.mucluster{use_day}{1}(r,:);
    [pxx,f] = periodogram(spk,[],[],Fs);
    pxxs(r,:) = 10*log10(pxx);
end
figure
set(gca,'Color', 'w')
mpxx = mean(pxxs);
plot(f,mpxx)
[pks,locs] = findpeaks(mpxx);
ym = ylim;
axis tight
line([f(locs(1)) f(locs(1))], [-100 100], 'Color', [.9 .9 .9]);
ylim(ym)
xlabel('Hz')
ylabel('dB')
hold on
plot(f(locs(1)),pks(1),'or')
text(f(locs(1)),pks(1),sprintf('  %.2f Hz',f(locs(1))))
xlim([4 40])
hold off
title([{'mean riptrig spike periodogram peak detection'};...
    {sprintf('%s %d %d', animal, use_day, nt_psth.ntrode)}])
%% Now make a wavelet with that peak frequency
% then plot the MU raster, with the KdensFR, (create wavelet with peak FFT),
% and finally plot wavelet-convolved signal.
% the phase will come from the analytic signal, created when creating the ...
% wavefiltered signal

%%dfs_phasecoherence
%{
So i could get this to work on the already created riptrig spiking
instantaneous firing rate (insead of lfp).. i just need to alter
gatherRipSnips so that it takes data field name as input so it's not locked to
looking for the ('data') field in the filter output.

gatherRipSnips takes indices and the filter output of riptriglfp and
returns allNTDataCat which get fed, along with wave config (getWaveParams),
into computeAnalyticSignal which also can save the AS.
-- the wave param set defines which frequencies, etc the wavelets are made
and datafiltered by
%}
andays = find(~cellfun(@isempty,data{1}.F.output)); %get nonempty eps
% basically rotate and vertstack to get ::ntrode x sample x trial
% need to create a version of gatherRipSnips specifcally for riptrigspiking
% because it is saved in Foutput differently than lfp is
% [index, allNTDataCat, dataByDay] = gatherRipSnips_spikes(...
%     data{1}.F.output{1}, 'datafield', 'instantFR');
%% combine the MU clusters of the riptrigspiking data
Fout = combine_mu_clusters(data{ian}.F);
% im putting this inside combine_riptrigspiking_perntrode
%%
out = combine_riptrigspiking_perntrode(data);
%% Now all (MU) datamat
ntrodes = [out{1}.ntrode];
alld = [];
for n = 1:length(out{1})
    try
        alld(n,:,:) = permute(cell2mat({out{1}(n).mudata.instantFR}'), [3 2 1]);
    catch
        fprintf('wut %d\n',n)
        continue
    end
end
%% hmm some D10 ntrodes are missing mu on different days.. (nt 9, 19, 29,30)
% maybe try to do the wavelet stuff for each ntrode individually
waveSet = 'riptrigMU_instantFR';
paths = make_paths(waveSet, env);
ntrodes = [out{1}.ntrode];
alld = [];
for n = 1:length(out{1})
    nt_mu_keys = cell2mat({out{1}(n).mudata.index}');
    nt_instaFR_alldayeps = permute(cell2mat({out{1}(n).mudata.instantFR}'), [3 2 1]);
    eplengths = cellfun(@(x) length(x(:,1)), {out{1}(n).mudata.instantFR}, 'un', 1);
    dayepvec = cell2mat(arrayfun(@(x,y,z) repmat([y z],x,1),eplengths',nt_mu_keys(:,1),nt_mu_keys(:,2),'un',0));
    AS(n).ntrode = ntrodes(n);
    AS(n).dayep_perrip = dayepvec;
    
    postwin = nt_instaFR_alldayeps(:,size(nt_instaFR_alldayeps,2)/2:end,:);
    wp = getWaveParams(waveSet, postwin);
    AS(n).postresult = computeAnalyticSignal(postwin, wp, animals{ian}, ...
        paths.filenamesave, 'saveAnalyticSignal', 0);
    
    prewin = nt_instaFR_alldayeps(:,1:size(nt_instaFR_alldayeps,2)/2,:);
    wp = getWaveParams(waveSet, prewin);
    AS(n).preresult = computeAnalyticSignal(prewin, wp, animals{ian}, ...
        paths.filenamesave, 'saveAnalyticSignal', 0);
    
    fullwin = nt_instaFR_alldayeps; %
    wp = getWaveParams(waveSet, fullwin);
    AS(n).result = computeAnalyticSignal(fullwin, wp, animals{ian}, ...
        paths.filenamesave, 'saveAnalyticSignal', 0);
end
%%
figure
usef = 3;
nt = 11;
day = 3;
ep1 = 2;
ep2 = 4;
subplot(2,1,1)
usedayeprips = ismember(AS(nt).dayep_perrip, [day ep1], 'rows');
imagesc(flipud(squeeze(AS(nt).result.ph(:,usedayeprips,usef))'))

subplot(2,1,2)
usedayeprips = ismember(AS(nt).dayep_perrip, [day ep2], 'rows');
imagesc(flipud(squeeze(AS(nt).result.ph(:,usedayeprips,usef))'))

%% plot all of the epochs
figure
imagesc(flipud(squeeze(AS(nt).result.ph(:,:,3))'))
%%
figure
for f = 1:6
    subplot(1,6,f)
    imagesc(flipud(squeeze(AS(nt).result.ph(:,:,f))'))
end
colormap('jet')

%% now plot an ep of the phase data
usef = 3;
nt = 11;
day = 3;
ep = 2;
usedayeprips = ismember(AS(nt).dayep_perrip, [day ep], 'rows');
subplot(2,2,1)
imagesc(flipud(squeeze(AS(nt).preresult.ph(:,usedayeprips,usef))'))
subplot(2,2,2)
imagesc(flipud(squeeze(AS(nt).postresult.ph(:,usedayeprips,usef))'))
colormap('jet')
subplot(2,2,3)
ph = squeeze(AS(nt).preresult.ph(:,usedayeprips,usef)');
prephripstart = ph(:,end);
scatter(prephripstart,1:length(prephripstart), 'filled')
% lsline
axis tight
subplot(2,2,4)
ph = squeeze(AS(nt).postresult.ph(:,usedayeprips,usef)');
phripstart = ph(:,1);
scatter(1:length(phripstart),abs(phripstart), 'filled')
axis tight
lsline
% [r, p] =corrcoef(1:length(phripstart),abs(phripstart))
%% now get a scatter of those phase values with swr num
ph = squeeze(AS(nt).result.ph(:,usedayeprips,usef)');
phripstart = ph(:,length(ph(1,:))/2);
scatter(length(phripstart):-1:1,phripstart, 'filled')
lsline
%% It kinda looks like the phase that the swr happens on is switching instead
% of gradually changing.. maybe look to see if that is change in task phase
% or correct/incorrect trial?

[eState] = calcEventPerformanceState(F, animals, behavestruct);
%%
% figure
% th = linspace(0,360,50);
% r = 0.005*th/10;
% th_radians = deg2rad(th);
cmap = jet(length(1:length(phripstart)));
clf
p = polar(phripstart',length(phripstart):-1:1, '.');

t = 0:pi/20:2*pi;
plot(t,sin(2*t),'-mo',...
    'LineWidth',2,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',[.49 1 .63],...
    'MarkerSize',10)

%% plot the MEC intra-swr response magnitude over single days
% if the phase shift is there then it should also be apparent in the
% gradual change in the response magnitude
% how to do this? mean instantFR within first 100ms?
% then later i can use the actual swr bounds

n = 11;
day = 3;
ep = 2;

nt_mu_keys = cell2mat({out{1}(n).mudata.index}');
usedayeprips = find(ismember(nt_mu_keys(:,[1 2]), [day ep], 'rows'));
a = out{1}(n).mudata(usedayeprips).instantFR;
figure(4)
imagesc(a)
colormap(jet)
%%
figure(13)
a = cell2mat({out{1}(n).mudata.instantFR}');
b = log(a+1);
c = b(:,3:end-2);
i = imagesc(c);
colormap(winter)
set(gca,'color', 'w')
ylabel('swr#')
xlabel('time')
title('FiringRate riptrigFR wtrack alldays D10 nt11')
% i.xtick(-.5:.001:.5)
colorbar
caxis([0 4]);

%%
figure
imagesc(squeeze(AS(nt).result.ph(:,:,6))')
colormap(hot)
set(gca,'color', 'w')
ylabel('swr#')
xlabel('time')
title('PHASE waveletConvFR wtrack alldays D10 nt11 ')
colorbar
%%
figure
imagesc(log10(abs(squeeze(AS(nt).result.as(:,:,6)))'))
colormap(jet)
set(gca,'color', 'w')
ylabel('swr#')
xlabel('time')
title('POWER waveletConvFR wtrack alldays D10 nt11 ')
colorbar

%% maybe overlay an alpha map of the analytic signal to highlight the
% phase locked periods more

%% do all ntrodes FR for D10, then the other animals
filenamesave = 'riptrigFR_overdays';
epenv = 'wtrack';
pausefigs = 1;
savefigs = 0;
%%
Pp = load_plotting_params(filenamesave);
paths = make_paths(filenamesave, epenv);
for ani = 1:length(animals)
    animal = animals{ani}
    for n = 1:length(ntrodes)
        %% ---- init fig----
        if savefigs && ~pausefigs;
            close all
            ifig = figure('Visible','off','units','normalized','position', ...
                Pp.position);
        else
            ifig = figure('units','normalized','position',Pp.position);
        end
        set(gcf,'color','white')
        %%
        nt_mu_keys = cell2mat({out{1}(n).mudata.index}');
        eplengths = cellfun(@(x) length(x(:,1)), {out{1}(n).mudata.instantFR}, 'un', 1);
        dayepvec = cell2mat(arrayfun(@(x,y,z) repmat([y z],x,1),eplengths',nt_mu_keys(:,1),nt_mu_keys(:,2),'un',0));
        a = cell2mat({out{1}(n).mudata.instantFR}');
        b = log(a+1); % adding one so that there's no negative log firing rates
        c = b(:,3:end-2);
        im = imagesc(c);
        colormap(hot)
        colorbar
        ylabel('swr#')
        xlabel('time')

        %% Super Axes
        sprtitleax = axes('Position',[0 0 1 1],'Visible','off', 'Parent', ifig);
        sprtit = sprintf('%s %s %s nt%d', filenamesave, epenv, animal, ntrodes(n));
        iStitle = text(.5, .98, {sprtit}, 'Parent', sprtitleax, 'Units', ...
            'normalized');
        set(iStitle,'FontWeight','bold','Color','k', 'FontName', 'Arial', ...
            'horizontalAlignment', 'center','FontSize', 14);
        
        %% ---- pause, save figs ----
        if pausefigs
            pause
        end
        if savefigs
            save_figure(paths.figdirectory, paths.filenamesave, sprtit)
            close all
        end
    end
end















