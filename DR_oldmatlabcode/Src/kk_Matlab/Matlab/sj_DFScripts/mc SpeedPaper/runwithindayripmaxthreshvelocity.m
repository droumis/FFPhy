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

% Time selection
timefilter = {{'getriptimes', '($nripples >= 1)', [], 'tetfilter', 'isequal($area,''CA1'') & $numcells>2'}};

% Create and Run filter velocity
a = createfilter('animal',animals,'epochs',epochfilterA,'excludetime',timefilter, 'iterator',iterator);
a = setfilterfunction(a,'calcvelocitysegment',{'pos'},'smooth',1);
a = runfilter(a);
velA = a;

a = setfilterfunction(a,'getmaxthresh',{'ripples','cellinfo'},'tetlist','isequal($area,''CA1'') & $numcells>2');
a = runfilter(a);
ripA = a;

b = createfilter('animal',animals,'epochs',epochfilterB,'excludetime',timefilter, 'iterator',iterator);
b = setfilterfunction(b,'calcvelocitysegment',{'pos'},'smooth',1);
b = runfilter(b);
velB = b;

b = setfilterfunction(b,'getmaxthresh',{'ripples','cellinfo'},'tetlist','isequal($area,''CA1'') & $numcells>2');
b = runfilter(b);
ripB = b;

save('/data13/mcarr/VelocityPaper/Ripples/velA.mat','velA')
save('/data13/mcarr/VelocityPaper/Ripples/velB.mat','velB')
save('/data13/mcarr/VelocityPaper/Ripples/ripA.mat','ripA')
save('/data13/mcarr/VelocityPaper/Ripples/ripB.mat','ripB')

%% Load data

load '/data13/mcarr/VelocityPaper/Ripples/velA.mat'
load '/data13/mcarr/VelocityPaper/Ripples/velB.mat'
load '/data13/mcarr/VelocityPaper/Ripples/ripA.mat'
load '/data13/mcarr/VelocityPaper/Ripples/ripB.mat'

%% COMPARE NOVEL AND FAMILIAR WITHIN DAY
bin = log([1/4 1 4 16]);
velocityA = [];
threshA = [];
velocityB = [];
threshB = [];

for an = 2:length(velA)
    count = 0;
    for d = 1:length(velA(an).output)
        if ~isempty(velA(an).output{d}) && ~isempty(velB(an).output{d})
            count = count+1;
            if count == 1
            velocityA = [velocityA; log(velA(an).output{d}.velocity)];
            tmpA = nan(size(ripA(an).output{d}.maxthresh));
            for e = 1:size(ripA(an).output{d}.maxthresh,2)
                tmpA(:,e) = ripA(an).output{d}.maxthresh(:,e).*ripA(an).output{d}.std(e) + ripA(an).output{d}.baseline(e);
            	tmpA(:,e) = (tmpA(:,e)-ripA(an).output{d}.totalbaseline(e))./ripA(an).output{d}.totalstd(e);
            end
            if an == 2
                tmpA = tmpA(:,1:6);
            end
            threshA = [threshA; nanmean(tmpA,2)];
            velocityB = [velocityB; log(velB(an).output{d}.velocity)];
            tmpB = nan(size(ripB(an).output{d}.maxthresh));
            for e = 1:size(ripB(an).output{d}.maxthresh,2)
            	tmpB(:,e) = ripB(an).output{d}.maxthresh(:,e).*ripB(an).output{d}.std(e) + ripB(an).output{d}.baseline(e);
                tmpB(:,e) = (tmpB(:,e)-ripB(an).output{d}.totalbaseline(e))./ripB(an).output{d}.totalstd(e);
            end
            threshB = [threshB; nanmean(tmpB,2)];
            end
        end
    end
end

invalid = isinf(velocityA) | isnan(velocityA) | velocityA<log(1/8) | isnan(threshA) | threshA>15;
threshA(invalid) = [];
threshA = threshA - 3;
velocityA(invalid) = [];
subsF = lookup(velocityA,bin);

invalid = isinf(velocityB) | isnan(velocityB) | velocityB<log(1/8) | isnan(threshB) | threshB>15;
threshB(invalid) = [];
threshB = threshB - 3;
velocityB(invalid) = [];
subsN = lookup(velocityB,bin);

RN = accumarray(subsN,threshB,[length(bin) 1],@(x) nanmean(x),NaN);
RF = accumarray(subsF,threshA,[length(bin) 1],@(x) nanmean(x),NaN);
errN = accumarray(subsN,threshB,[length(bin) 1],@(x) nanstd(x)./sqrt(sum(~isnan(x))-1),NaN);
errF = accumarray(subsF,threshA,[length(bin) 1],@(x) nanstd(x)./sqrt(sum(~isnan(x))-1),NaN);

figure
barerror([1 3 5 7],RN,errN,0.5,0.5,'r');
barerror([2 4 6 8],RF,errF,0.5,0.5);
set(gca,'xlim',[0 9],'ylim',[-0.5 1.5],'xtick',1.5:2:7.5,'xticklabel',[1/4 1 4 16])
xlabel('Speed')
ylabel('Normalized ripple power')

pvalue = nan(length(bin),1);
for i = 1:length(bin)
    pvalue(i,1) = ranksum(threshA(subsF==i),threshB(subsN==i));
end
%--------------------------------------------------------------------------
% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/SpeedPaper/%d_%d_%d_withinday_ripplepowergroup.pdf', m, d, y);
print('-dpdf', savestring)

%% PLOT SLOPE FOR EACH ANIMAL
for an = 2:length(velA)
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
            tmpA = nan(size(ripA(an).output{d}.maxthresh));
            for e = 1:size(ripA(an).output{d}.maxthresh,2)
                tmpA(:,e) = ripA(an).output{d}.maxthresh(:,e).*ripA(an).output{d}.std(e) + ripA(an).output{d}.baseline(e);
            	tmpA(:,e) = (tmpA(:,e)-ripA(an).output{d}.totalbaseline(e))./ripA(an).output{d}.totalstd(e);
            end
            if an == 2
                tmpA = tmpA(:,1:6);
            end
            threshA = [threshA; nanmean(tmpA,2)];
            velocityB = [velocityB; log(velB(an).output{d}.velocity)];
            tmpB = nan(size(ripB(an).output{d}.maxthresh));
            for e = 1:size(ripB(an).output{d}.maxthresh,2)
            	tmpB(:,e) = ripB(an).output{d}.maxthresh(:,e).*ripB(an).output{d}.std(e) + ripB(an).output{d}.baseline(e);
                tmpB(:,e) = (tmpB(:,e)-ripB(an).output{d}.totalbaseline(e))./ripB(an).output{d}.totalstd(e);
            end
            threshB = [threshB; nanmean(tmpB,2)];
            end
        end
    end

    invalid = isinf(velocityA) | isnan(velocityA) | velocityA<log(1/8) | isnan(threshA) | threshA>15;
    velocityA(invalid) = [];
    threshA(invalid) = [];
    
    invalid = isinf(velocityB) | isnan(velocityB) | velocityB<log(1/8) | isnan(threshB) | threshB>15;
    velocityB(invalid) = [];
    threshB(invalid) = [];
    
    figure(4)
    [Ab Abint] = regress(threshA,[ones(length(velocityA),1) velocityA],0.05);
    [Bb Bbint] = regress(threshB,[ones(length(velocityB),1) velocityB],0.05);
    hold on
    plot([1 2],[Bb(2) Ab(2)],'k', [1 2], [Bb(2) Ab(2)],'ko','MarkerFace','k')
    set(gca,'xtick',[1 2],'xticklabel','Novel|Familiar','xlim',[0.8 2.2],'FontSize',18)
    set(gca,'yLim',[-0.3 0.0],'ytick',-0.3:0.1:0.1,'FontSize',18)
    ylabel('Slope','FontSize',18)
    title({'Within Day Comparison of Slope' , 'Speed vs. Ripple Amplitude'})
    box off

end
%Significant at p < 0.05 by regression

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/SpeedPaper/%d_%d_%d_withinday_rippleslope_peranimal.pdf', m, d, y);
print('-dpdf', savestring)

%% PLOT GROUP DATA, subsampling to control for differences in behavior
bin = log([1/4 1 4 16]);

nboot = 1000;
bootslopeA = nan(nboot,1);
bootslopeB = nan(nboot,1);
bootmeanA = nan(nboot,length(bin));
bootmeanB = nan(nboot,length(bin));
bootdifference = nan(nboot,length(bin));
totalcount = 10;
for s = 1:nboot
velocityA = [];
threshA = [];
velocityB = [];
threshB = [];
    for an = 2:length(velA)
        exposure = 0;
        for d = 1:length(velA(an).output)
            if ~isempty(velA(an).output{d}) && ~isempty(velB(an).output{d})
                exposure = exposure+1;
                if exposure == 1
                tmpv = log(velA(an).output{d}.velocity);
                tmpA = nan(size(ripA(an).output{d}.maxthresh));
                for e = 1:size(ripA(an).output{d}.maxthresh,2)
                    tmpA(:,e) = ripA(an).output{d}.maxthresh(:,e).*ripA(an).output{d}.std(e) + ripA(an).output{d}.baseline(e);
                    tmpA(:,e) = (tmpA(:,e)-ripA(an).output{d}.totalbaseline(e))./ripA(an).output{d}.totalstd(e);
                end
                tmpt = nanmean(tmpA,2);
                tmpt(isinf(tmpv) | isnan(tmpv) | tmpv<log(1/8)) = [];
                tmpv(isinf(tmpv) | isnan(tmpv) | tmpv<log(1/8)) = [];
                tmpv(isnan(tmpt)) = []; tmpt(isnan(tmpt)) = [];
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
                tmpB = nan(size(ripB(an).output{d}.maxthresh));
                for e = 1:size(ripB(an).output{d}.maxthresh,2)
                    tmpB(:,e) = ripB(an).output{d}.maxthresh(:,e).*ripB(an).output{d}.std(e) + ripB(an).output{d}.baseline(e);
                    tmpB(:,e) = (tmpB(:,e)-ripB(an).output{d}.totalbaseline(e))./ripB(an).output{d}.totalstd(e);
                end
                tmpt = nanmean(tmpB,2);
                tmpt(isinf(tmpv) | isnan(tmpv) | tmpv<log(1/8)) = [];
                tmpv(isinf(tmpv) | isnan(tmpv) | tmpv<log(1/8)) = [];
                tmpv(isnan(tmpt)) = []; tmpt(isnan(tmpt)) = [];
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
    end
    subs = lookup(velocityA,bin);
    b = regress(threshA,[ones(length(velocityA),1) velocityA],0.05);
    bootslopeA(s) = b(2);
    [b bint] = regress(threshB,[ones(length(velocityB),1) velocityB],0.05);
    bootslopeB(s) = b(2);
end

% Plot slope
figure
seLB = mean(bootslopeB)-prctile(bootslopeB,2.5);
seUB = -mean(bootslopeB)+prctile(bootslopeB,97.5);
seLA = mean(bootslopeA)-prctile(bootslopeA,2.5);
seUA = -mean(bootslopeA)+prctile(bootslopeA,97.5);
errorbar2([1 2],mean([bootslopeB bootslopeA]),[seLB seUB; seLA seUA],0.1,'k')
hold on
bar([1 2],mean([bootslopeB bootslopeA]))
set(gca,'xtick',[1 2],'xticklabel','Novel|Familiar','FontSize',18)
set(gca,'xlim',[0.5 2.5],'ylim',[-0.45 0.05],'ytick',-0.5:0.1:0.05)
ylabel('Slope')

% By bootstrap statistics, p < 0.05
%--------------------------------------------------------------------------
% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/SpeedPaper/%d_%d_%d_withinday_rippleslope.pdf', m, d, y);
print('-dpdf', savestring)
