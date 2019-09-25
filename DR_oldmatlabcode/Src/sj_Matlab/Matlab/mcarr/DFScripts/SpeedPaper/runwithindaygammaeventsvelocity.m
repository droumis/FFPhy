%% Run Filter

%animal selection
animals = {'Corriander','Eight','Five','Miles','Ten'};

% Epoch selection
epochfilterA = [];
epochfilterB = [];
for i = 1:12
    epochfilterA{i} = ['($experimentday == ',num2str(i),') & $dailyexposure == 1 & isequal($description,''TrackA'')'];
    epochfilterB{i} = ['($experimentday == ',num2str(i),') & $dailyexposure == 1 & isequal($description,''TrackB'')'];
end

% Select iterator
iterator = 'singleepochanal';

%% Run Filter for Slow Gamma

% Time selection
timefilter = {{'getgammatimes', '($ngamma >= 1)', [], 'low','tetfilter', 'isequal($area,''CA1'') & $numcells>1','minthresh',2}, ...
      {'getgammatimes', '($ngamma == 0)', [], 'high','tetfilter', 'isequal($area,''CA1'') & $numcells>1','minthresh',2}};

% Create and Run filter velocity
a = createfilter('animal',animals,'epochs',epochfilterA,'excludetime',timefilter, 'iterator',iterator);
a = setfilterfunction(a,'calcvelocitysegment',{'pos'},'smooth',1);
a = runfilter(a);
velA = a;

a = setfilterfunction(a,'getmaxthresh',{'gammal','cellinfo'},'tetlist','isequal($area,''CA1'') & $numcells >=5');
a = runfilter(a);
gamA = a;

b = createfilter('animal',animals,'epochs',epochfilterB,'excludetime',timefilter, 'iterator',iterator);
b = setfilterfunction(b,'calcvelocitysegment',{'pos'},'smooth',1);
b = runfilter(b);
velB = b;

b = setfilterfunction(b,'getmaxthresh',{'gammal','cellinfo'},'tetlist','isequal($area,''CA1'') & $numcells>=5');
b = runfilter(b);
gamB = b;

save('/data13/mcarr/VelocityPaper/SlowGamma/velA.mat','velA')
save('/data13/mcarr/VelocityPaper/SlowGamma/velB.mat','velB')
save('/data13/mcarr/VelocityPaper/SlowGamma/gamA.mat','gamA')
save('/data13/mcarr/VelocityPaper/SlowGamma/gamB.mat','gamB')

%% Run Filter for Fast Gamma

% Time selection
timefilter = {{'getgammatimes', '($ngamma == 0)', [], 'low','tetfilter', 'isequal($area,''CA1'') & $numcells>1','minthresh',2}, ...
      {'getgammatimes', '($ngamma >= 1)', [], 'high','tetfilter', 'isequal($area,''CA1'') & $numcells>1','minthresh',2}};

% Create and Run filter velocity
a = createfilter('animal',animals,'epochs',epochfilterA,'excludetime',timefilter, 'iterator',iterator);
a = setfilterfunction(a,'calcvelocitysegment',{'pos'},'smooth',1);
a = runfilter(a);
velA = a;

a = setfilterfunction(a,'getmaxthresh',{'gammah','cellinfo'},'tetlist','isequal($area,''CA1'') & $numcells>=5');
a = runfilter(a);
gamA = a;

b = createfilter('animal',animals,'epochs',epochfilterB,'excludetime',timefilter, 'iterator',iterator);
b = setfilterfunction(b,'calcvelocitysegment',{'pos'},'smooth',1);
b = runfilter(b);
velB = b;

b = setfilterfunction(b,'getmaxthresh',{'gammah','cellinfo'},'tetlist','isequal($area,''CA1'') & $numcells >=5');
b = runfilter(b);
gamB = b;

save('/data13/mcarr/VelocityPaper/FastGamma/velA.mat','velA')
save('/data13/mcarr/VelocityPaper/FastGamma/velB.mat','velB')
save('/data13/mcarr/VelocityPaper/FastGamma/gamA.mat','gamA')
save('/data13/mcarr/VelocityPaper/FastGamma/gamB.mat','gamB')


%% Load data

% To analyze slow gamma
load '/data13/mcarr/VelocityPaper/SlowGamma/velA.mat'
load '/data13/mcarr/VelocityPaper/SlowGamma/velB.mat'
load '/data13/mcarr/VelocityPaper/SlowGamma/gamA.mat'
load '/data13/mcarr/VelocityPaper/SlowGamma/gamB.mat'

% To analyze fast gamma
load '/data13/mcarr/VelocityPaper/FastGamma/velA.mat'
load '/data13/mcarr/VelocityPaper/FastGamma/velB.mat'
load '/data13/mcarr/VelocityPaper/FastGamma/gamA.mat'
load '/data13/mcarr/VelocityPaper/FastGamma/gamB.mat'

%% PLOT SLOPE FOR EACH ANIMAL
for an = 1:length(velA)
    velocityA = [];
    threshA = [];
    velocityB = [];
    threshB = [];

    count = 0;
    for d = 1:length(velA(an).output)
        if ~isempty(velA(an).output{d}) && ~isempty(velB(an).output{d})
            count = count+1;
            if count == 1
            velocityA = [velocityA; log(velA(an).output{d}.velocity)];
            tmpA = nan(size(gamA(an).output{d}.maxthresh));
            for e = 1:size(gamA(an).output{d}.maxthresh,2)
                tmpA(:,e) = gamA(an).output{d}.maxthresh(:,e).*gamA(an).output{d}.std(e) + gamA(an).output{d}.baseline(e);
            	tmpA(:,e) = (tmpA(:,e)-gamA(an).output{d}.totalbaseline(e))./gamA(an).output{d}.totalstd(e);
            end
            threshA = [threshA; nanmean(tmpA,2)];
            velocityB = [velocityB; log(velB(an).output{d}.velocity)];
            tmpB = nan(size(gamB(an).output{d}.maxthresh));
            for e = 1:size(gamB(an).output{d}.maxthresh,2)
            	tmpB(:,e) = gamB(an).output{d}.maxthresh(:,e).*gamB(an).output{d}.std(e) + gamB(an).output{d}.baseline(e);
                tmpB(:,e) = (tmpB(:,e)-gamB(an).output{d}.totalbaseline(e))./gamB(an).output{d}.totalstd(e);
            end
            threshB = [threshB; nanmean(tmpB,2)];
            end
        end
    end

    
    threshA(isinf(velocityA) | isnan(velocityA) | velocityA<log(3) | velocityA>4) = [];
    velocityA(isinf(velocityA) | isnan(velocityA) | velocityA<log(3) | velocityA>4) = [];
    threshB(isinf(velocityB) | isnan(velocityB) | velocityB<log(3) | velocityB>4) = [];
    velocityB(isinf(velocityB) | isnan(velocityB) | velocityB<log(3) | velocityB>4) = [];
    velocityA(isnan(threshA) | threshA>10 | threshA<1) = []; threshA(isnan(threshA) | threshA>10 | threshA<1) = [];
    velocityB(isnan(threshB) | threshB>10 | threshB<1) = []; threshB(isnan(threshB) | threshB>10 | threshB<1) = [];

    figure(1)
    [Ab Abint] = regress(threshA,[ones(length(velocityA),1) velocityA],0.05);
    [Bb Bbint] = regress(threshB,[ones(length(velocityB),1) velocityB],0.05);
    hold on
    plot([1 2],[Bb(2) Ab(2)],'k', [1 2], [Bb(2) Ab(2)],'ko','MarkerFace','k')
    set(gca,'xtick',[1 2],'xticklabel','Novel|Familiar','xlim',[0.8 2.2])
    set(gca,'yLim',[-0.15 0.25],'ytick',-0.2:0.1:0.2,'FontSize',18)
    ylabel('Slope','FontSize',18)
    title({'Within Day Comparison of Slope' , 'Speed vs. Slow Gamma Amplitude'})
    box off
end

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/SpeedPaper/%d_%d_%d_withindayexample_fasttgammaslope.pdf', m, d, y);
print('-dpdf', savestring)

%% PLOT GROUP DATA, subsampling to control for differences in behavior
bin = log([2 4 8 16]);
tmpcount = nan(length(velA),2);
for an = 1:length(velA)
    exposure = 0;
    for d = 1:length(velA(an).output)
        if ~isempty(velA(an).output{d}) && ~isempty(velB(an).output{d})
            exposure = exposure+1;
            if exposure == 1
                tmpv = log(velA(an).output{d}.velocity);
                tmpA = nan(size(gamA(an).output{d}.maxthresh));
                for e = 1:size(gamA(an).output{d}.maxthresh,2)
                    tmpA(:,e) = gamA(an).output{d}.maxthresh(:,e).*gamA(an).output{d}.std(e) + gamA(an).output{d}.baseline(e);
                    tmpA(:,e) = (tmpA(:,e)-gamA(an).output{d}.totalbaseline(e))./gamA(an).output{d}.totalstd(e);
                end
                tmpt = nanmean(tmpA,2);
                tmpt(isinf(tmpv) | isnan(tmpv) | tmpv<log(1/4) | tmpv>4) = [];
                tmpv(isinf(tmpv) | isnan(tmpv) | tmpv<log(1/4) | tmpv>4) = [];
                tmpv(isnan(tmpt) | tmpt>10 | tmpt < 1) = [];
                tmpt(isnan(tmpt) | tmpt>10 | tmpt < 1) = [];
                subs = lookup(tmpv,bin);
                tmpcount(an,1) = min(hist(subs,length(bin)));
                tmpv = log(velB(an).output{d}.velocity);
                tmpB = nan(size(gamB(an).output{d}.maxthresh));
                for e = 1:size(gamB(an).output{d}.maxthresh,2)
                    tmpB(:,e) = gamB(an).output{d}.maxthresh(:,e).*gamB(an).output{d}.std(e) + gamB(an).output{d}.baseline(e);
                    tmpB(:,e) = (tmpB(:,e)-gamB(an).output{d}.totalbaseline(e))./gamB(an).output{d}.totalstd(e);
                end
                tmpt = nanmean(tmpB,2);                
                tmpt(isinf(tmpv) | isnan(tmpv) | tmpv<log(1/4) | tmpv>4) = [];
                tmpv(isinf(tmpv) | isnan(tmpv) | tmpv<log(1/4) | tmpv>4) = [];
                tmpv(isnan(tmpt) | tmpt>10 | tmpt < 1) = [];
                tmpt(isnan(tmpt) | tmpt>10 | tmpt < 1) = [];
                subs = lookup(tmpv,bin);
                tmpcount(an,2) = min(hist(subs,length(bin)));
            end
        end
    end
end

nboot = 1000;
bootslopeA = nan(nboot,1);
bootslopeB = nan(nboot,1);
bootmeanA = nan(nboot,length(bin));
bootmeanB = nan(nboot,length(bin));
bootdifference = nan(nboot,length(bin));
bootoveralldifference = nan(nboot,1);
count = min(tmpcount,[],2);
totalcount = min(count);
for s = 1:nboot
velocityA = [];
threshA = [];
velocityB = [];
threshB = [];
    for an = 1:length(velA)
        exposure = 0;
        for d = 1:length(velA(an).output)
            if ~isempty(velA(an).output{d}) && ~isempty(velB(an).output{d})
                exposure = exposure+1;
                if exposure == 1
                tmpv = log(velA(an).output{d}.velocity);
                tmpA = nan(size(gamA(an).output{d}.maxthresh));
                for e = 1:size(gamA(an).output{d}.maxthresh,2)
                    tmpA(:,e) = gamA(an).output{d}.maxthresh(:,e).*gamA(an).output{d}.std(e) + gamA(an).output{d}.baseline(e);
                    tmpA(:,e) = (tmpA(:,e)-gamA(an).output{d}.totalbaseline(e))./gamA(an).output{d}.totalstd(e);
                end
                tmpt = nanmean(tmpA,2);
                tmpt(isinf(tmpv) | isnan(tmpv) | tmpv<log(1/4) | tmpv>4) = [];
                tmpv(isinf(tmpv) | isnan(tmpv) | tmpv<log(1/4) | tmpv>4) = [];
                tmpv(isnan(tmpt) | tmpt>10 | tmpt < 1) = [];
                tmpt(isnan(tmpt) | tmpt>10 | tmpt < 1) = [];
                subs = lookup(tmpv,bin);
                boot = ceil(length(tmpv)*rand(2*length(tmpv),1));
                for q = 1:max(subs)
                    tmpboot = find(subs(boot)==q);
                    tmpboot = boot(tmpboot);
                    tmp = tmpv(tmpboot);
                    tmp = tmp(1:totalcount);
                    velocityA = [velocityA; tmp];
                    tmp = tmpt(tmpboot);
                    tmp = tmp(1:totalcount);
                    threshA = [threshA; tmp];
                end

                tmpv = log(velB(an).output{d}.velocity);
                tmpB = nan(size(gamB(an).output{d}.maxthresh));
                for e = 1:size(gamB(an).output{d}.maxthresh,2)
                    tmpB(:,e) = gamB(an).output{d}.maxthresh(:,e).*gamB(an).output{d}.std(e) + gamB(an).output{d}.baseline(e);
                    tmpB(:,e) = (tmpB(:,e)-gamB(an).output{d}.totalbaseline(e))./gamB(an).output{d}.totalstd(e);
                end
                tmpt = nanmean(tmpB,2);
                tmpt(isinf(tmpv) | isnan(tmpv) | tmpv<log(1/4) | tmpv>4) = [];
                tmpv(isinf(tmpv) | isnan(tmpv) | tmpv<log(1/4) | tmpv>4) = [];
                tmpv(isnan(tmpt) | tmpt>10 | tmpt < 1) = [];
                tmpt(isnan(tmpt) | tmpt>10 | tmpt < 1) = [];

                subs = lookup(tmpv,bin);
                boot = ceil(length(tmpv)*rand(2*length(tmpv),1));
                for q = 1:max(subs)
                    tmpboot = find(subs(boot)==q);
                    tmpboot = boot(tmpboot);
                    tmp = tmpv(tmpboot);
                    tmp = tmp(1:totalcount);
                    velocityB = [velocityB; tmp];
                    tmp = tmpt(tmpboot);
                    tmp = tmp(1:totalcount);
                    threshB = [threshB; tmp];
                end
                clear tmpA tmpB
                end
            end
        end
    end
    subs = lookup(velocityA,bin);
	bootmeanA(s,:) = accumarray(subs,threshA,[max(subs) 1],@(x) nanmean(x),NaN);
    bootmeanB(s,:) = accumarray(subs,threshB,[max(subs) 1],@(x) nanmean(x),NaN);
    tmp = [];
    for e = 1:max(subs)
        bootdifference(s,e) = nanmean(threshB(subs==e)-threshA(subs==e));
        tmp = [tmp; threshB(subs==e)-threshA(subs==e)];
    end
    bootoveralldifference(s) = nanmean(tmp);
end

% Plot mean
figure
hold on
seLB = mean(bootmeanB)-prctile(bootmeanB,2.5);
seUB = -mean(bootmeanB)+prctile(bootmeanB,97.5);
seLA = mean(bootmeanA)-prctile(bootmeanA,2.5);
seUA = -mean(bootmeanA)+prctile(bootmeanA,97.5);
errorbar2(bin,mean(bootmeanB),[seLB' seUB']',0.025,'c')
errorbar2(bin,mean(bootmeanA),[seLA' seUA']',0.025,'b')
legend('Novel','Familiar')
plot(bin,mean(bootmeanB),'c',bin,mean(bootmeanB),'co','MarkerFace','c')
plot(bin,mean(bootmeanA),'b',bin, mean(bootmeanA),'bo','MarkerFace','b')
set(gca,'xtick',bin,'xticklabel',exp(bin),'FontSize',18)
set(gca,'ylim',[1.9 2.4],'ytick',2:0.1:5)
ylabel('Mean Slow Gamma Power')
xlabel('Speed cm/sec')

% Plot difference
A = mean(bootdifference);
seL = A-prctile(bootdifference,2.5);
seU = -A+prctile(bootdifference,97.5);
figure
hold on
errorbar2(bin,A,[seL' seU']',0.001,'k')
plot(bin,A,'ko','MarkerFace','k')
set(gca,'xtick',bin,'xtickLabel',exp(bin))
set(gca,'yLim',[-0.3 .3],'ytick',-0.3:0.1:1.5,'FontSize',18)
xlabel('Speed (cm/sec)','FontSize',18)
ylabel('Delta Slow Gamma Power: Novel - Familiar','FontSize',18)
box off

%--------------------------------------------------------------------------
% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/SpeedPaper/%d_%d_%d_withinday_slowgammapower.pdf', m, d, y);
print('-dpdf', savestring)

%% PLOT WITHIN DAY COMPARISON FOR EACH ANIMAL
bin = log([4 8 16]);
tmpcount = nan(length(velA),2);
for an = 1:length(velA)
    exposure = 0;
    for d = 1:length(velA(an).output)
        if ~isempty(velA(an).output{d}) && ~isempty(velB(an).output{d})
            exposure = exposure+1;
             if exposure == 1
                tmpv = log(velA(an).output{d}.velocity);
                tmpA = nan(size(gamA(an).output{d}.maxthresh));
                for e = 1:size(gamA(an).output{d}.maxthresh,2)
                    tmpA(:,e) = gamA(an).output{d}.maxthresh(:,e).*gamA(an).output{d}.std(e) + gamA(an).output{d}.baseline(e);
                    tmpA(:,e) = (tmpA(:,e)-gamA(an).output{d}.totalbaseline(e))./gamA(an).output{d}.totalstd(e);
                end
                tmpt = nanmean(tmpA,2);
                tmpt(isinf(tmpv) | isnan(tmpv) | tmpv<log(3) | tmpv>4) = [];
                tmpv(isinf(tmpv) | isnan(tmpv) | tmpv<log(3) | tmpv>4) = [];
                tmpv(isnan(tmpt) | tmpt>10 | tmpt < 1) = [];
                tmpt(isnan(tmpt) | tmpt>10 | tmpt < 1) = [];
                subs = lookup(tmpv,bin);
                boot = ceil(length(tmpv)*rand(length(tmpv),1));
                tmpcount(an,1) = min(hist(subs,length(bin)));
                
                tmpv = log(velB(an).output{d}.velocity);
                tmpB = nan(size(gamB(an).output{d}.maxthresh));
                for e = 1:size(gamB(an).output{d}.maxthresh,2)
                    tmpB(:,e) = gamB(an).output{d}.maxthresh(:,e).*gamB(an).output{d}.std(e) + gamB(an).output{d}.baseline(e);
                    tmpB(:,e) = (tmpB(:,e)-gamB(an).output{d}.totalbaseline(e))./gamB(an).output{d}.totalstd(e);
                end
                tmpt = nanmean(tmpB,2);
                tmpt(isinf(tmpv) | isnan(tmpv) | tmpv<log(3) | tmpv>4) = [];
                tmpv(isinf(tmpv) | isnan(tmpv) | tmpv<log(3) | tmpv>4) = [];
                tmpv(isnan(tmpt) | tmpt>10 | tmpt<1) = [];
                tmpt(isnan(tmpt) | tmpt>10 | tmpt<1) = [];
                subs = lookup(tmpv,bin);
                tmpcount(an,2) = min(hist(subs,length(bin)));
            end
        end
    end
end
nboot = 500;
for an = 1:length(velA)
    bootslopeA = nan(nboot,1);
    bootslopeB = nan(nboot,1);
    bootmeanA = nan(nboot,length(bin));
    bootmeanB = nan(nboot,length(bin));
    bootdifference = nan(nboot,length(bin));
    bootoveralldifference = nan(nboot,1);
    count = min(tmpcount,[],2);
    totalcount = count(an);
    for s = 1:nboot
        velocityA = [];
        threshA = [];
        velocityB = [];
        threshB = [];
        exposure = 0;
        for d = 1:length(velA(an).output)
            if ~isempty(velA(an).output{d}) && ~isempty(velB(an).output{d})
                exposure = exposure+1;
                if exposure == 1
                tmpv = log(velA(an).output{d}.velocity);
                tmpA = nan(size(gamA(an).output{d}.maxthresh));
                for e = 1:size(gamA(an).output{d}.maxthresh,2)
                    tmpA(:,e) = gamA(an).output{d}.maxthresh(:,e).*gamA(an).output{d}.std(e) + gamA(an).output{d}.baseline(e);
                    tmpA(:,e) = (tmpA(:,e)-gamA(an).output{d}.totalbaseline(e))./gamA(an).output{d}.totalstd(e);
                end
                tmpt = nanmean(tmpA,2);
                tmpt(isinf(tmpv) | isnan(tmpv) | tmpv<log(3) | tmpv>4) = [];
                tmpv(isinf(tmpv) | isnan(tmpv) | tmpv<log(3) | tmpv>4) = [];
                tmpv(isnan(tmpt) | tmpt>10 | tmpt<1) = [];
                tmpt(isnan(tmpt) | tmpt>10 | tmpt<1) = [];
                subs = lookup(tmpv,bin);
                boot = ceil(length(tmpv)*rand(3*length(tmpv),1));
                for q = 1:max(subs)
                    tmpboot = find(subs(boot)==q);
                    tmpboot = boot(tmpboot);
                    tmp = tmpv(tmpboot);
                    tmp = tmp(1:totalcount);
                    velocityA = [velocityA; tmp];
                    tmp = tmpt(tmpboot);
                    tmp = tmp(1:totalcount);
                    threshA = [threshA; tmp];
                end
                
                tmpv = log(velB(an).output{d}.velocity);
                tmpB = nan(size(gamB(an).output{d}.maxthresh));
                for e = 1:size(gamB(an).output{d}.maxthresh,2)
                    tmpB(:,e) = gamB(an).output{d}.maxthresh(:,e).*gamB(an).output{d}.std(e) + gamB(an).output{d}.baseline(e);
                    tmpB(:,e) = (tmpB(:,e)-gamB(an).output{d}.totalbaseline(e))./gamB(an).output{d}.totalstd(e);
                end
                tmpt = nanmean(tmpB,2);
                tmpt(isinf(tmpv) | isnan(tmpv) | tmpv<log(3) | tmpv>4) = [];
                tmpv(isinf(tmpv) | isnan(tmpv) | tmpv<log(3) | tmpv>4) = [];
                tmpv(isnan(tmpt) | tmpt>10 | tmpt<1) = [];
                tmpt(isnan(tmpt) | tmpt>10 | tmpt<1) = [];
                subs = lookup(tmpv,bin);
                boot = ceil(length(tmpv)*rand(3*length(tmpv),1));
                for q = 1:max(subs)
                    tmpboot = find(subs(boot)==q);
                    tmpboot = boot(tmpboot);
                    tmp = tmpv(tmpboot);
                    tmp = tmp(1:totalcount);
                    velocityB = [velocityB; tmp];
                    tmp = tmpt(tmpboot);
                    tmp = tmp(1:totalcount);
                    threshB = [threshB; tmp];
                end
                clear tmpA tmpB
                end
            end
        end
        subs = lookup(velocityA,bin);
        bootmeanA(s,:) = accumarray(subs,threshA,[max(subs) 1],@(x) nanmean(x),NaN);
        bootmeanB(s,:) = accumarray(subs,threshB,[max(subs) 1],@(x) nanmean(x),NaN);

%       b = regress(threshA,[ones(length(velocityA),1) velocityA],0.05);
%       bootslopeA(s) = b(2);
%       [b bint] = regress(threshB,[ones(length(velocityB),1) velocityB],0.05);
%       bootslopeB(s) = b(2);
        tmp = [];
        for e = 1:max(subs)
            bootdifference(s,e) = nanmean(threshB(subs==e)-threshA(subs==e));
            tmp = [tmp; threshB(subs==e)-threshA(subs==e)];
        end
        bootoveralldifference(s) = nanmean(tmp);
    end
    
%   % Plot slopes
%   figure(1)
%   hold on
%   seLB = mean(bootslopeB)-prctile(bootslopeB,2.5);
%   seUB = -mean(bootslopeB)+prctile(bootslopeB,97.5);
%   seLA = mean(bootslopeA)-prctile(bootslopeA,2.5);
%   seUA = -mean(bootslopeA)+prctile(bootslopeA,97.5);
%   errorbar2([1 2],mean([bootslopeB bootslopeA]),[seLB seUB; seLA seUA],0.001,'k')
%   plot([1 2],mean([bootslopeB bootslopeA]),'k',[1 2],mean([bootslopeB bootslopeA]),'ko','MarkerFace','k')
%   set(gca,'xtick',[1 2],'xticklabel','Novel|Familiar','FontSize',18)
%   ranksum(bootslopeB,bootslopeA);
%   set(gca,'xlim',[0.5 2.5],'ylim',[-0.2 -0.05],'ytick',-0.2:0.05:-0.05)
%   ylabel('Slope')
 
    % Plot mean
    figure(1+an)
    hold on
    seLB = mean(bootmeanB)-prctile(bootmeanB,2.5);
    seUB = -mean(bootmeanB)+prctile(bootmeanB,97.5);
    seLA = mean(bootmeanA)-prctile(bootmeanA,2.5);
    seUA = -mean(bootmeanA)+prctile(bootmeanA,97.5);
    errorbar2(bin,mean(bootmeanB),[seLB' seUB']',0.025,'m')
    errorbar2(bin,mean(bootmeanA),[seLA' seUA']',0.025,'r')
    legend('Novel','Familiar')
    plot(bin,mean(bootmeanB),'m',bin,mean(bootmeanB),'mo','MarkerFace','m')
    plot(bin,mean(bootmeanA),'r',bin, mean(bootmeanA),'ro','MarkerFace','r')
    set(gca,'xtick',bin,'xticklabel',exp(bin),'FontSize',18)
    set(gca,'ylim',[1.9 2.4],'ytick',2:0.1:5)
    ylabel('Mean Fast Gamma Power')
    xlabel('Speed cm/sec')
    
    %Plot differences
    figure(1)
    A = mean(bootdifference);
    seL = A-prctile(bootdifference,2.5);
    seU = -A+prctile(bootdifference,97.5);
    hold on
    errorbar2(bin,A,[seL' seU']',0.001,'k')
    plot(bin,A,'ko','MarkerFace','k')
    set(gca,'xtick',bin,'xtickLabel',exp(bin))
    set(gca,'yLim',[0 1.5],'ytick',[0:0.5:1.5],'FontSize',20)
    xlabel('Speed (cm/sec)','FontSize',22)
    ylabel({'Mean Difference in Novel Environment Power'; 'Compared to Familiar Environment'},'FontSize',22)
    box off
     
     prctile(bootoveralldifference,2.5)
end

