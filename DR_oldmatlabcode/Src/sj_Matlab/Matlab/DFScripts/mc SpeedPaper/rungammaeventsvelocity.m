%% Run Filter

% Track A
% Animal Selection
animals = {'Conley','Corriander','Dudley','Eight','Five','Miles','Ten'};

% Epoch selection
epochfilter = [];
for i = 1:14
    epochfilter{i} = ['($exposure == ',num2str(i),') & $dailytracksexperienced == 1 & isequal($description,''TrackA'')'];
end

% Select iterator
iterator = 'singleepochanal';

%--------------------------------------------------------------------------
% Slow Gamma
% Time selection
timefilter = {{'getgammatimes', '($ngamma >= 1)', [], 'low','tetfilter', 'isequal($area,''CA1'') & $numcells>1','minthresh',2}, ...
      {'getgammatimes', '($ngamma == 0)', [], 'high','tetfilter', 'isequal($area,''CA1'') & $numcells>1','minthresh',2}};

% Create and Run filter velocity
l = createfilter('animal',animals,'epochs',epochfilter,'excludetime',timefilter, 'iterator',iterator);
l = setfilterfunction(l,'calcvelocitysegment',{'pos'},'smooth',1);
l = runfilter(l);

velocity = l;
l = setfilterfunction(l,'getmaxthresh',{'gammal','cellinfo'},'tetlist','isequal($area,''CA1'') & $numcells>=5');
l = runfilter(l);
gammaxthresh = l;
clear l

save('/data13/mcarr/VelocityPaper/SlowGamma/velocity.mat','velocity');
save('/data13/mcarr/VelocityPaper/SlowGamma/gammaxthresh.mat','gammaxthresh');

%--------------------------------------------------------------------------
% Fast Gamma
% Time selection
timefilter = {{'getgammatimes', '($ngamma >= 1)', [], 'high','tetfilter', 'isequal($area,''CA1'') & $numcells>1','minthresh',2}, ...
     {'getgammatimes', '($ngamma == 0)', [], 'low','tetfilter', 'isequal($area,''CA1'') & $numcells>1','minthresh',2}};
 
% Create and Run filter velocity
h = createfilter('animal',animals,'epochs',epochfilter,'excludetime',timefilter, 'iterator',iterator);
h = setfilterfunction(h,'calcvelocitysegment',{'pos'},'smooth',1);
h = runfilter(h);

velocity = h;
h = setfilterfunction(h,'getmaxthresh',{'gammah','cellinfo'},'tetlist','isequal($area,''CA1'') & $numcells>=5');
h = runfilter(h);
gammaxthresh = h;
clear h

save('/data13/mcarr/VelocityPaper/FastGamma/velocity.mat','velocity');
save('/data13/mcarr/VelocityPaper/FastGamma/gammaxthresh.mat','gammaxthresh');

%% Load filters

% To analyze slow gamma
load '/data13/mcarr/VelocityPaper/SlowGamma/velocity.mat'
vel = velocity;
clear velocity

load '/data13/mcarr/VelocityPaper/SlowGamma/gammaxthresh.mat'
g = gammaxthresh;
clear gammaxthresh

%--------------------------------------------------------------------------
% To analyze high gamma
load '/data13/mcarr/VelocityPaper/FastGamma/velocity.mat'
vel = velocity;
clear velocity

load '/data13/mcarr/VelocityPaper/FastGamma/gammaxthresh.mat'
g = gammaxthresh;
clear gammaxthresh

%% Plot Group speed vs. threshold for exposure 1

bin = log([4 8 16 32]);
velocity = [];
thresh = [];
% Combine data across all animals for day 1
for an = 1:length(vel)
    if ~isempty(vel(an).output)
        for d = 1
            if ~isempty(vel(an).output{d})
                velocity = [velocity; log(vel(an).output{d}.velocity)];
                tmp = nan(size(g(an).output{d}.maxthresh));
                for e = 1:size(g(an).output{d}.maxthresh,2)
                    tmp(:,e) = (g(an).output{d}.maxthresh(:,e)*g(an).output{d}.std(e)) + g(an).output{d}.baseline(e);
                    tmp(:,e) = (tmp(:,e) - g(an).output{d}.totalbaseline(e))./g(an).output{d}.totalstd(e);
                end
                thresh = [thresh; nanmean(tmp,2)];
                % Get rid of NaN & Inf, and invalid velocity
                 thresh(isinf(velocity) | isnan(velocity) | velocity<log(3) | velocity>4) = [];
                 velocity(isinf(velocity) | isnan(velocity) | velocity<log(3) | velocity>4) = [];
                 velocity(isnan(thresh) | thresh<1 | thresh>9) = [];
                 thresh(isnan(thresh) | thresh<1 | thresh>9) = [];
            end
        end
    end
end

% Compute means with confidence bounds and best fit line
subs = lookup(velocity,bin);
A = accumarray(subs,thresh,[length(bin) 1],@(x) nanmean(x),NaN);
Lstd = accumarray(subs,thresh,[length(bin) 1],@(x) nanstd(x),NaN);
Lsqrt = accumarray(subs,thresh,[length(bin) 1],@(x) sqrt(sum(~isnan(x))),NaN);
seL = -Lstd./Lsqrt;
seU = +Lstd./Lsqrt;
[b,bint,r,rint,stats] = regress(thresh,[ones(length(velocity),1) velocity],0.05);

x = [bin(1) bin(end)];
LL = min([bint(2,1)*x(1)+bint(1,1), bint(2,1)*x(1)+bint(1,2), ...
    bint(2,2)*x(1)+bint(1,1), bint(2,2)*x(1)+bint(1,2)]);
LU = max([bint(2,1)*x(1)+bint(1,1), bint(2,1)*x(1)+bint(1,2), ...
    bint(2,2)*x(1)+bint(1,1), bint(2,2)*x(1)+bint(1,2)]);
RU = max([bint(2,2)*x(2)+bint(1,1), bint(2,2)*x(2)+bint(1,2), ...
    bint(2,1)*x(2)+bint(1,1), bint(2,1)*x(2)+bint(1,2)]);
RL = min([bint(2,2)*x(2)+bint(1,1), bint(2,2)*x(2)+bint(1,2), ...
    bint(2,1)*x(2)+bint(1,1), bint(2,1)*x(2)+bint(1,2)]);

% Plot regression with 95% conf limits and means with conf bounds

% For slow gamma
figure
fill([x(1) x(1) x(2) x(2)], [LL LU RU RL],'b','FaceAlpha',1,'EdgeColor','none')
hold on
plot(bin,bin*b(2)+b(1),'b')    
plot(bin,A,'bo','MarkerFace','b')
errorbar2(bin,A,[seL seU]',0.001,'b','plottype','semilogx')
set(gca,'xtick',bin,'xtickLabel',exp(bin),'xLim',[bin(1)-0.1 bin(end)+0.1])
xlabel('Speed (cm/sec)','FontSize',18)
ylabel('Slow Gamma Normalized Power','FontSize',18)
set(gca,'FontSize',18,'ylim',[1.875 2.55],'ytick',[1.5:0.25:3])
box off

% For fast gamma
figure
fill([x(1) x(1) x(2) x(2)], [LL LU RU RL],'r','FaceAlpha',1,'EdgeColor','none')
hold on
plot(bin,bin*b(2)+b(1),'r')    
plot(bin,A,'ro','MarkerFace','r')
errorbar2(bin,A,[seL seU]',0.001,'r','plottype','semilogx')
set(gca,'xtick',bin,'xtickLabel',exp(bin),'xLim',[bin(1)-0.1 bin(end)+0.1])
xlabel('Speed (cm/sec)','FontSize',18)
ylabel('Fast Gamma Normalized Power','FontSize',18)
set(gca,'FontSize',18,'ylim',[1.875 2.55],'ytick',[1.5:0.25:3])
box off

% Save figure: Slow Gamma
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/SpeedPaper/%d_%d_%d_velocity_slowgammathresh_exposure1.pdf', m, d, y);
print('-dpdf', savestring)
% b = [2.5047 -0.1405] bint = [2.4614 2.5480; -0.1611 -0.1199];
% stats:
% R = 0.0095 F = 179.1502 p = 0

% Save figure: Fast Gamma
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/SpeedPaper/%d_%d_%d_velocity_fastgammathresh_exposure1.pdf', m, d, y);
print('-dpdf', savestring)
% b = [2.0179 0.1217] bint = [1.9828 2.0512; 0.1064 0.1369];
% stats:
% R = 0.0089 F = 244.3877 p = 0

%% LOOK AT CHANGE IN SLOPE OVER DAYS FOR GROUP DATA

bin = log([4 8 16 32]);

velocity = cell(14,1);
thresh = cell(14,1);
for an =1:length(g)
    for d = 1:length(g(an).output)
        if ~isempty(g(an).output{d})
            velocity{d} =[velocity{d}; log(vel(an).output{d}.velocity)];
            tmp = nan(size(g(an).output{d}.maxthresh));
            for e = 1:size(g(an).output{d}.maxthresh,2)
            	tmp(:,e) = (g(an).output{d}.maxthresh(:,e)*g(an).output{d}.std(e)) + g(an).output{d}.baseline(e);
                tmp(:,e) = (tmp(:,e) - g(an).output{d}.totalbaseline(e))./g(an).output{d}.totalstd(e);
            end
            thresh{d} = [thresh{d}; nanmean(tmp,2)];
            thresh{d}(isinf(velocity{d}) | isnan(velocity{d}) | velocity{d}<log(3)) = [];
            velocity{d}(isinf(velocity{d}) | isnan(velocity{d}) | velocity{d}<log(3)) = [];
            velocity{d}(isnan(thresh{d}) | thresh{d}<1 | thresh{d}>9) = [];
            thresh{d}(isnan(thresh{d}) | thresh{d}<1 | thresh{d}>9) = [];
        end
    end
end
nboot = 1000;
Bboot = nan(nboot,length(velocity));
for d = 1:length(velocity)
    if ~isempty(velocity{d})
        for s = 1:nboot
            boot = ceil(length(velocity{d})*rand(length(velocity{d}),1));
            xboot = velocity{d}(boot);
            yboot = thresh{d}(boot);
            b = regress(yboot,[ones(length(xboot),1) xboot],0.05);
            Bboot(s,d) = b(2);
        end
    end
end

subs = repmat([1:length(velocity)],nboot,1);
subs = reshape(subs,nboot*length(velocity),1);
boot = reshape(Bboot,nboot*length(velocity),1);
novelboot = boot(subs<=10);
novelsubs = subs(subs<=10);
A = accumarray(novelsubs,novelboot,[max(novelsubs) 1],@(x) nanmean(x),NaN);
seL = A-accumarray(novelsubs,novelboot,[max(novelsubs) 1],@(x) prctile(x,[2.5]));
seU = A-accumarray(novelsubs,novelboot,[max(novelsubs) 1],@(x) prctile(x,[97.5]));
%--------------------------------------------------------------------------
% Plot slope over days for Slow Gamma
figure
hold on
plot(1:max(novelsubs),A,'bo','MarkerFace','b')
plot(1:max(novelsubs),A,'b')
errorbar2(1:max(novelsubs),A,[seL seU]',0.001,'b')
xlabel('Exposure to Environment','FontSize',18)
ylabel('Slope of Speed vs. Normalized Slow Gamma Power','FontSize',18)
set(gca,'xtick',1:max(novelsubs),'xlim',[0 max(novelsubs)+1])
set(gca,'ylim',[-0.2 0],'ytick',[-0.2:0.05:0],'FontSize',18)
box off

% Plot slope over days for Fast Gamma
figure
hold on
plot(1:max(novelsubs),A,'ro','MarkerFace','r')
plot(1:max(novelsubs),A,'r')
errorbar2(1:max(novelsubs),A,[seL seU]',0.001,'r')
xlabel('Exposure to Environment','FontSize',18)
ylabel('Slope of Speed vs. Normalized Fast Gamma Power','FontSize',18)
set(gca,'xtick',1:max(novelsubs),'xlim',[0 max(novelsubs)+1])
set(gca,'ylim',[-0.05 0.3],'ytick',[-0.2:0.05:0],'FontSize',18)
box off
%--------------------------------------------------------------------------
% Save figure
% For slow gamma
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/SpeedPaper/%d_%d_%d_slopes_velvsslowgamma_exposures1_10.pdf', m, d, y);
print('-dpdf', savestring)

% For fast gamma
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/SpeedPaper/%d_%d_%d_slopes_velvsfastgamma_exposures1_10.pdf', m, d, y);
print('-dpdf', savestring)


%% LOOK AT CHANGE IN SLOPE OVER DAYS FOR GROUP DATA, CONTROL FOR BEHAVIOR

bin = log([4 8 16]);

count = [];
for an = 1:length(vel)
    for d = 1:length(vel(an).output)
        if ~isempty(vel(an).output{d})
        	tmpv = log(vel(an).output{d}.velocity);
            tmp = nan(size(g(an).output{d}.maxthresh));
            for e = 1:size(g(an).output{d}.maxthresh,2)
            	tmp(:,e) = (g(an).output{d}.maxthresh(:,e)*g(an).output{d}.std(e)) + g(an).output{d}.baseline(e);
                tmp(:,e) = (tmp(:,e) - g(an).output{d}.totalbaseline(e))./g(an).output{d}.totalstd(e);
            end
            tmpt = nanmean(tmp,2);
            tmpt(isinf(tmpv) | isnan(tmpv)) = []; tmpv(isinf(tmpv) | isnan(tmpv)) = [];
            tmpt(tmpv<log(3)) = []; tmpv(tmpv<log(3)) = [];
            tmpv(tmpt<1 | tmpt>10) = []; tmpt(tmpt<1 | tmpt>10) = [];
            subs = lookup(tmpv,bin);
        	count= [count; min(hist(subs,length(bin)))];
        end
    end
end
count = min(count(count>0));

nboot = 1000;
qnboot = 1000;
Bboot = nan(nboot,14);
for s = 1:nboot
    velocity = cell(14,1);
    thresh = cell(14,1);
    for an =1:length(g)
        for d = 1:length(g(an).output)
            if ~isempty(g(an).output{d})
                tmpv = log(vel(an).output{d}.velocity);
                tmp = nan(size(g(an).output{d}.maxthresh));
                for e = 1:size(g(an).output{d}.maxthresh,2)
                    tmp(:,e) = (g(an).output{d}.maxthresh(:,e)*g(an).output{d}.std(e)) + g(an).output{d}.baseline(e);
                    tmp(:,e) = (tmp(:,e) - g(an).output{d}.totalbaseline(e))./g(an).output{d}.totalstd(e);
                end
                tmpt = nanmean(tmp,2);
                tmpt(isinf(tmpv) | isnan(tmpv)) = []; tmpv(isinf(tmpv) | isnan(tmpv)) = [];
                tmpt(tmpv<log(3)) = []; tmpv(tmpv<log(3)) = [];
                tmpv(tmpt<1 | tmpt>10) = []; tmpt(tmpt<1 | tmpt>10) = [];
                
                if length(tmpv)>count
                    subs = lookup(tmpv,bin);
                    boot = ceil(length(tmpv)*rand(2*length(tmpv),1));
                    for q = 1:max(subs)
                        tmpboot = find(subs(boot)==q);
                        tmpboot = boot(tmpboot);
                        tmp = tmpv(tmpboot);
                        tmp = tmp(1:count);
                        velocity{d} = [velocity{d}; tmp];
                        tmp = tmpt(tmpboot);
                        tmp = tmp(1:count);
                        thresh{d} = [thresh{d}; tmp];
                    end
                end
            end
        end
    end
    for d = 1:length(velocity)
        if ~isempty(velocity{d})
            b = regress(thresh{d},[ones(length(velocity{d}),1) velocity{d}],0.05);
            Bboot(s,d) = b(2);
        end
    end
end

subs = repmat([1:length(velocity)],nboot,1);
subs = reshape(subs,nboot*length(velocity),1);
boot = reshape(Bboot,nboot*length(velocity),1);
novelboot = boot(subs<=10);
novelsubs = subs(subs<=10);
[b bint r rint stats] = regress(novelboot,[ones(length(novelsubs),1) novelsubs],0.05);
permr = nan(qnboot,1);
for q = 1:qnboot
    perm = novelsubs(randperm(length(novelsubs)));
    [tmpb bint r rint tmpstats] = regress(novelboot,[ones(length(novelsubs),1) perm],0.05);
    permr(q) = tmpb(2);
end
pvalue = sum(b(2)>permr)./qnboot;
A = accumarray(novelsubs,novelboot,[max(novelsubs) 1],@(x) nanmean(x),NaN);
seL = A-accumarray(novelsubs,novelboot,[max(novelsubs) 1],@(x) prctile(x,[2.5]));
seU = -A+accumarray(novelsubs,novelboot,[max(novelsubs) 1],@(x) prctile(x,[97.5]));

%--------------------------------------------------------------------------
% Plot slope over days
% For slow gamma
figure
hold on
plot(1:max(novelsubs),A,'bo','MarkerFace','b')
plot(1:max(novelsubs),A,'b')
errorbar2(1:max(novelsubs),A,[seL seU]',0.25,'b')
xlabel('Exposure to Environment','FontSize',18)
ylabel('Slope of Speed vs. Normalized Slow Gamma Power','FontSize',18)
set(gca,'xtick',1:max(novelsubs),'xlim',[0 max(novelsubs)+1])
set(gca,'ylim',[-0.3 0.1],'ytick',[-0.5:0.1:0.3],'FontSize',18)
box off

% For fast gamma
figure
hold on
plot(1:max(novelsubs),A,'ro','MarkerFace','r')
plot(1:max(novelsubs),A,'r')
errorbar2(1:max(novelsubs),A,[seL seU]',0.25,'r')
xlabel('Exposure to Environment','FontSize',18)
ylabel('Slope of Speed vs. Normalized Fast Gamma Power','FontSize',18)
set(gca,'xtick',1:max(novelsubs),'xlim',[0 max(novelsubs)+1])
set(gca,'ylim',[-0.1 0.3],'ytick',[-0.5:0.1:0.3],'FontSize',18)
box off

% Save figure
% For slow gamma
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/SpeedPaper/%d_%d_%d_slopes_velvsslowgamma_exposures1_10.pdf', m, d, y);
print('-dpdf', savestring)

% For fast gamma
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/SpeedPaper/%d_%d_%d_slopes_velvsfastgamma_exposures1_10.pdf', m, d, y);
print('-dpdf', savestring)

%% LOOK AT CHANGE IN SLOPE OVER DAYS FOR INDIVIDUAL ANIMALS

bin = log([4 8 16 32]);

for an = 1:length(g)
velocity = cell(14,1);
thresh = cell(14,1);
    for d = 1:length(g(an).output)
        if ~isempty(g(an).output{d})
            velocity{d} =[velocity{d}; log(vel(an).output{d}.velocity)];
            tmp = nan(size(g(an).output{d}.maxthresh));
            for e = 1:size(g(an).output{d}.maxthresh,2)
            	tmp(:,e) = (g(an).output{d}.maxthresh(:,e)*g(an).output{d}.std(e)) + g(an).output{d}.baseline(e);
                tmp(:,e) = (tmp(:,e) - g(an).output{d}.totalbaseline(e))./g(an).output{d}.totalstd(e);
            end
            thresh{d} = [thresh{d}; nanmean(tmp,2)];
            thresh{d}(isinf(velocity{d}) | isnan(velocity{d})) = []; velocity{d}(isinf(velocity{d}) | isnan(velocity{d})) = [];
            thresh{d}(velocity{d}<log(3)) = []; velocity{d}(velocity{d}<log(3)) = [];
            velocity{d}(thresh{d}<1 | thresh{d}>9) = []; thresh{d}(thresh{d}<1 | thresh{d}>9) = [];
        end
    end
    
    nboot = 10;
    Bboot = nan(nboot,length(velocity));
    for d = 1:length(velocity)
        if ~isempty(velocity{d})
            for s = 1:nboot
                boot = ceil(length(velocity{d})*rand(length(velocity{d}),1));
                xboot = velocity{d}(boot);
                yboot = thresh{d}(boot);
                b = regress(yboot,[ones(length(xboot),1) xboot],0.05);
                Bboot(s,d) = b(2);
            end
        end
    end

    subs = repmat([1:length(velocity)],nboot,1);
    subs = reshape(subs,nboot*length(velocity),1);
    boot = reshape(Bboot,nboot*length(velocity),1);
    novelboot = boot(subs<=7);
    novelsubs = subs(subs<=7);
    A = accumarray(novelsubs,novelboot,[max(novelsubs) 1],@(x) nanmean(x),NaN);
    seL = A-accumarray(novelsubs,novelboot,[max(novelsubs) 1],@(x) prctile(x,[2.5]));
    seU = A-accumarray(novelsubs,novelboot,[max(novelsubs) 1],@(x) prctile(x,[97.5]));

    %--------------------------------------------------------------------------
    % Plot slope over days
    figure
    hold on
    plot(1:max(novelsubs),A,'ko','MarkerFace','k')
    plot(1:max(novelsubs),A,'k')
    errorbar2(1:max(novelsubs),A,[seL seU]',0.001,'k')
    xlabel('Exposure to Track B','FontSize',24)
    ylabel('Slope of Speed vs. Normalized Gamma Power','FontSize',24)
    set(gca,'xtick',1:max(novelsubs),'xlim',[0 max(novelsubs)+1])
    set(gca,'ylim',[-0.2 0],'ytick',[-0.2:0.1:0],'FontSize',20)
    box off
end
%--------------------------------------------------------------------------
% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/SpeedPaper/%d_%d_%d_slopes_velvsthresh_exposures1_10.pdf', m, d, y);
print('-dpdf', savestring)

%% LOOK AT CHANGE IN SLOPE OVER DAYS FOR INDIVIDUAL ANIMAL, CONTROL FOR BEHAVIOR

bin = log([4 8 16]);

count = [];
for an = 1:length(vel)
    for d = 1:length(vel(an).output)
    	if ~isempty(vel(an).output{d})
        	tmpv = log(vel(an).output{d}.velocity);
            tmp = nan(size(g(an).output{d}.maxthresh));
            for e = 1:size(g(an).output{d}.maxthresh,2)
            	tmp(:,e) = (g(an).output{d}.maxthresh(:,e)*g(an).output{d}.std(e)) + g(an).output{d}.baseline(e);
                tmp(:,e) = (tmp(:,e) - g(an).output{d}.totalbaseline(e))./g(an).output{d}.totalstd(e);
            end
            tmpt = nanmean(tmp,2);
            tmpt(isinf(tmpv) | isnan(tmpv)) = []; tmpv(isinf(tmpv) | isnan(tmpv)) = [];
            tmpt(tmpv<log(3)) = []; tmpv(tmpv<log(3)) = [];
            tmpv(tmpt<1 | tmpt>9) = []; tmpt(tmpt<1 | tmpt>9) = [];
            subs = lookup(tmpv,bin);
        	count= [count; min(hist(subs,length(bin)))];
    	end
    end
    
    count = min(count(count>0));

    nboot = 100;
    Bboot = nan(nboot,14);
    for s = 1:nboot
        velocity = cell(14,1);
        thresh = cell(14,1);
        for d = 1:length(g(an).output)
            if ~isempty(g(an).output{d})
                tmpv = log(vel(an).output{d}.velocity);
                tmp = nan(size(g(an).output{d}.maxthresh));
                for e = 1:size(g(an).output{d}.maxthresh,2)
                    tmp(:,e) = (g(an).output{d}.maxthresh(:,e)*g(an).output{d}.std(e)) + g(an).output{d}.baseline(e);
                    tmp(:,e) = (tmp(:,e) - g(an).output{d}.totalbaseline(e))./g(an).output{d}.totalstd(e);
                end
                tmpt = nanmean(tmp,2);
                tmpt(isinf(tmpv) | isnan(tmpv)) = []; tmpv(isinf(tmpv) | isnan(tmpv)) = [];
                tmpt(tmpv<log(3)) = []; tmpv(tmpv<log(3)) = [];
                tmpv(tmpt<1 | tmpt>9) = []; tmpt(tmpt<1 | tmpt>9) = [];
                
                if length(tmpv)>count
                    subs = lookup(tmpv,bin);
                    boot = ceil(length(tmpv)*rand(2*length(tmpv),1));
                    for q = 1:max(subs)
                        tmpboot = find(subs(boot)==q);
                        tmpboot = boot(tmpboot);
                        tmp = tmpv(tmpboot);
                        tmp = tmp(1:totalcount);
                        velocity{d} = [velocity{d}; tmp];
                        tmp = tmpt(tmpboot);
                        tmp = tmp(1:totalcount);
                        thresh{d} = [thresh{d}; tmp];
                    end
                end
            	b = regress(thresh{d},[ones(length(velocity{d}),1) velocity{d}],0.05);
                Bboot(s,d) = b(2);
            end
        end
    end
    
    subs = repmat([1:length(velocity)],nboot,1);
    subs = reshape(subs,nboot*length(velocity),1);
    boot = reshape(Bboot,nboot*length(velocity),1);
    novelboot = boot(subs<=7);
    novelsubs = subs(subs<=7);
    A = accumarray(novelsubs,novelboot,[max(novelsubs) 1],@(x) nanmean(x),NaN);
    seL = A-accumarray(novelsubs,novelboot,[max(novelsubs) 1],@(x) prctile(x,[2.5]));
    seU = A-accumarray(novelsubs,novelboot,[max(novelsubs) 1],@(x) prctile(x,[97.5]));
    [b,bint,r,rint,stats] = regress(novelboot,[ones(length(novelsubs),1) novelsubs],0.05);

    x = [1 max(novelsubs)];
    LL = min([bint(2,1)*x(1)+bint(1,1), bint(2,1)*x(1)+bint(1,2), ...
        bint(2,2)*x(1)+bint(1,1), bint(2,2)*x(1)+bint(1,2)]);
    LU = max([bint(2,1)*x(1)+bint(1,1), bint(2,1)*x(1)+bint(1,2), ...
        bint(2,2)*x(1)+bint(1,1), bint(2,2)*x(1)+bint(1,2)]);
    RU = max([bint(2,2)*x(2)+bint(1,1), bint(2,2)*x(2)+bint(1,2), ...
        bint(2,1)*x(2)+bint(1,1), bint(2,1)*x(2)+bint(1,2)]);
    RL = min([bint(2,2)*x(2)+bint(1,1), bint(2,2)*x(2)+bint(1,2), ...
        bint(2,1)*x(2)+bint(1,1), bint(2,1)*x(2)+bint(1,2)]);

    %--------------------------------------------------------------------------
    % Plot regression with 95% conf limits and means with conf bounds
    figure
    fill([x(1) x(1) x(2) x(2)], [LL LU RU RL],'b','FaceAlpha',1,'EdgeColor','none')
    hold on
    plot(1:max(novelsubs),[1:max(novelsubs)]*b(2)+b(1),'b')    
    plot(1:max(novelsubs),A,'ko','MarkerFace','k')
    errorbar2(1:max(novelsubs),A,[seL seU]',0.001,'k')
    xlabel('Exposure to Track A','FontSize',24)
    ylabel('Slope of Speed vs. Normalized Slow Gamma Power','FontSize',24)
    set(gca,'xtick',1:max(novelsubs),'xlim',[0 max(novelsubs)+1])
    %set(gca,'ytick',-0.25:0.25:0.25,'ylim',[-0.25 0.25],'FontSize',20)
    box off
end
