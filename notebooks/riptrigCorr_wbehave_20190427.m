

animal = 'D13';
aninfo = animaldef(animal);
ffdir = aninfo{2};

%%
day = 6;
ep = 4;
tetA = 9;
tetB = 4;
clusta = 9;
clustb = 6;
%% load data
spikes = loaddatastruct(ffdir, animal, 'spikes', day);
rips = loaddatastruct(ffdir, animal, 'ca1rippleskons', day);
cellinfo = loaddatastruct(ffdir, animal, 'cellinfo', day);

%%
clAspike = spikes{day}{ep}{tetA}{clusta};
clAinfo = cellinfo{day}{ep}{tetA}{clusta};
clBspike = spikes{day}{ep}{tetB}{clustb};
clBinfo = cellinfo{day}{ep}{tetB}{clustb};

%%
clf
bins = linspace(0, 1, 200);
subplot(2,2,1)
dfa = diff(clAspike.data);
histogram(dfa,bins, 'FaceColor', 'r', 'FaceAlpha', .2)
title(sprintf('%s day %d ep %d - nt %d  cl %d',animal, day,ep,tet,clusta))
xlim([0 .3])
ylabel('count')
xlabel('ISI (s)')
subplot(2,2,2)
df = diff(clBspike.data);
histogram(df,bins, 'FaceColor', 'b', 'FaceAlpha', .2);
title(sprintf('%s day %d ep %d - nt %d  cl %d',animal, day,ep,tet,clustb))
xlim([0 .3])
ylabel('count')
xlabel('ISI (s)')

%% get rips and rip spiking
rips_start_end = [rips{day}{ep}{1}.starttime rips{day}{ep}{1}.endtime];
ripspikesA = clAspike.data(logical(isExcluded(clAwspike.data, rips_start_end)));
ripspikesB = clBspike.data(logical(isExcluded(clBspike.data, rips_start_end)));

%%

bins = linspace(0, 1, 200);
subplot(2,2,3)
dfa = diff(ripspikesA);
histogram(dfa,bins, 'FaceColor', 'r', 'FaceAlpha', .2)
title(sprintf('spikes in rips %s day %d ep %d - nt %d  cl %d',animal, day,ep,tet,clusta))
xlim([0 .3])
ylabel('count')
xlabel('ISI (s)')
subplot(2,2,4)
df = diff(ripspikesB);
histogram(df,bins, 'FaceColor', 'b', 'FaceAlpha', .2);
title(sprintf('spikes in rips %s day %d ep %d - nt %d  cl %d',animal, day,ep,tet,clustb))
xlim([0 .3])
ylabel('count')
xlabel('ISI (s)')


%% plot time of rip events
nr = 10;
patch([rips_start_end(1:nr,:)'; flipud(rips_start_end(1:nr,:)')], ...
    [zeros(1,nr); zeros(1,nr); ones(1,nr); ones(1,nr)], 'red', 'FaceAlpha', .5, 'EdgeAlpha', 0)
%% spikes in rips
usespikeinds = ~isExcluded(ispike.data, rips_start_end);
spikesinrip = ispike.data(usespikeinds);

%% plot spikes % too much.. will crash
% plot(spikesinrip, ones(length(spikesinrip))*.5, '.')

%% setup pipeline to screen operators across spiketrainas and lfp per event, gather, plot
%{
1. Load the behave state
2. get derivative of est prob correct
3. get time of each state change
4. load position time
5. interp est prob and its der to pos times
6. save per day file
7. load mu and su spikes for two tets
8. filter to just spikes in rip times
9. 0 lag xcorr full spike trains and get per day xcorr score
9b. 0 lag xcorr for each rip 
10. for each rip, get the est prob corr and deriv based on time knn
11. save as result

1. load result for all days 
2. collect across days for each tetAB pair
3. regress and score est corr/deriv with tetAB xcorr per rip
4. save as result with tetAB index

1. load tetAB results 
2. sort by r2

%}


%%
spikesA_incr = zeros(length(postime), 1);
sres = knnsearch(postime, ripspikesA);
spikesA_incr(sres) = 1;
A = spikesA_incr;
% spend = round(length(spikesA_incr)-mod(length(spikesA_incr),2));
% A = sum(reshape(spikesA_incr(1:spend),spend/2,2),2);

spikesB_incr = zeros(length(postime), 1);
sres = knnsearch(postime, spikesB_incr);
spikesB_incr(sres) = 1;
B = spikesB_incr;
% spend = round(length(spikesB_incr)-mod(length(spikesB_incr),2));
% B = sum(reshape(spikesB_incr(1:spend),spend/2,2),2);

[xcAB, lags] = xcorr(A,B,500);
zerolagxc = xcAB(lags==0);

out{day}{ep}.xcorr = xcAB;
out{day}{ep}.lags = lags;
out{day}{ep}.zerolagxc = zerolagxc;
%%
clf
figure(2)
plot(lags,xcAB)

%% get zerolag for each rip
for i = 1:length(rips_start_end)
    spIripA = ripspikesA([ripspikesA>rips_start_end(i,1) & ripspikesA<rips_start_end(i,2)]);
    spIripB = ripspikesB([ripspikesB>rips_start_end(i,1) & ripspikesB<rips_start_end(i,2)]);
    rtrigA = zeros(1000, 1);
    rtrigB = zeros(1000, 1);
    sresA= knnsearch(postime, spIripA);
    sresB = knnsearch(postime, spIripB);
    rtrigA(sres) = 1;
    rtrigB(sres) = 1;
    
    [xcAB, lags] = xcorr(rtrigA,rtrigB,1);
    zerolagxc = xcAB(lags==0);

%     out{day}{ep}.xcorr = xcAB;
%     out{day}{ep}.lags = lags;
    out{day}{ep}(i) = zerolagxc;
    error('this isnt working yet')
end

%%
pos = loaddatastruct(ffdir, animal, 'pos', day);
postime = pos{day}{ep}.data(:,1);
bs = load([ffdir, animal, 'BehaveState.mat']);
gt = bs.BehaveState.statespace.allbound(:,5) == day;

% [~,index_A,~] = intersect(bs.BehaveState.statespace.allbound(:,4:5),[day ep],'rows');
d = bs.BehaveState.statespace.allbound(gt,1:4);
dt = bs.BehaveState.statespace.allepsMat(gt, 3);
estcorr_interp = interp1(dt,d,postime);
%%
% gaussian(sigma,numpoints)
a = abs(diff(estcorr_interp));
sma = conv(a(:,1), gaussian(2,1000), 'same') *300000;

%%
% plot(dt,d(:,1))
clf
subaxis(3,1,1)
plot(postime, estcorr_interp(:,1))
ylabel('est prob corr')
title('est prob correct')

subaxis(3,1,2)
plot(postime(1:end-1), a(:,1)*300000)
ylabel('abs change in est prob corr')
xlabel('time')

subaxis(3,1,3)
plot(a(logical(spikesA_incr)))

%%
% 
% function out = runxcorr(index, excludetimes, spikes, pos)
% 
% day = index(1);
% ep = index(2);
% 
% 
% 
% end
% 
% 
% 
% 
% 
% 
