%% Run Filter

% CA1
%animal selection
animals = {'Bond','Conley','Corriander','Dudley','Eight','Five','Frank','Miles','Ten'};

% Epoch selection
epochfilter = [];
for i = 1:14
    epochfilter{i} = ['($exposure == ',num2str(i),') & $dailytracksexperienced == 1 & isequal($description,''TrackA'')'];
end

% Select iterator
iterator = 'singleepochanal';

% Time selection
timefilter = {{'getriptimes', '($nripples >= 1)', [], 'tetfilter', 'isequal($area,''CA1'')'}};

% Create and Run filter velocity
g = createfilter('animal',animals,'epochs',epochfilter,'excludetime',timefilter, 'iterator',iterator);
g = setfilterfunction(g,'calcvelocitysegment',{'pos'},'smooth',1);
g = runfilter(g);
velocity = g;

g = setfilterfunction(g,'getmaxthresh',{'ripples','cellinfo'},'tetlist','isequal($area,''CA1'')');
g = runfilter(g);
ripmaxthresh = g;

save('/data13/mcarr/VelocityPaper/Ripples/velocity.mat','velocity')
save('/data13/mcarr/VelocityPaper/Ripples/ripmaxthresh.mat','ripmaxthresh')

%% Load Filters

%CA1
load '/data13/mcarr/VelocityPaper/Ripples/velocity.mat'
vel = velocity;
clear velocity
load '/data13/mcarr/VelocityPaper/Ripples/ripmaxthresh.mat'
rip = ripmaxthresh;
clear ripmaxthresh

%% Plot example of velocity and threshold relationship for exposures 1
bin = log([1/4 1/2 1 2 4 8 16 32]);

velocity = [];
thresh = [];
for an = 1:length(vel)
    if ~isempty(vel(an).output)
        for d = 1
            if ~isempty(vel(an).output{d})
                velocity = [velocity; log(vel(an).output{d}.velocity)];
                thresh = [thresh; nanmean(rip(an).output{d}.maxthresh,2)];
            end
        end
    end
end
thresh(isinf(velocity) | isnan(velocity) | velocity < log(1/8)) = [];
velocity(isinf(velocity) | isnan(velocity) | velocity < log(1/8)) = [];

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

figure
fill([x(1) x(1) x(2) x(2)], [LL LU RU RL],'k','FaceAlpha',1,'EdgeColor','none')
hold on
plot(bin,bin*b(2)+b(1),'k')    
plot(bin,A,'ko','MarkerFace','k')
errorbar2(bin,A,[seL seU]',0.001,'k','plottype','semilogx')
set(gca,'xtick',bin,'xtickLabel',exp(bin))
set(gca,'yLim',[3.2 5.2],'ytick',[3:0.5:5.5],'FontSize',20)
xlabel('Speed (cm/sec)','FontSize',22)
ylabel('Normalized Ripple Power','FontSize',22)
box off

% b: [4.5374 -0.3306] bint: [4.4691 4.6058; -0.3670 -0.2942]
% stats:
% R = 0.104 F = 317.0303 p = 0
%--------------------------------------------------------------------------
% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/SpeedPaper/%d_%d_%d_velocity_ripple_exposure1.pdf', m, d, y);
print('-dpdf', savestring)

%% Look at the change in slopes over days
bin = log([1/4 1/2 1 2 4 8 16 32]);

velocity = cell(23,1);
thresh = cell(23,1);
for an =1:length(rip)
    for d = 1:length(rip(an).output)
        if ~isempty(rip(an).output{d})
            velocity{d} =[velocity{d}; log(vel(an).output{d}.velocity)];
            thresh{d} = [thresh{d}; nanmean(rip(an).output{d}.maxthresh,2)];
        end
    end
end
nboot = 1000;
Bboot = nan(nboot,length(velocity));
B0boot = nan(nboot,length(velocity));
for d = 1:length(velocity)
    if ~isempty(velocity{d})
            thresh{d}(isinf(velocity{d}) | isnan(velocity{d}) | velocity{d}<log(1/8)) = [];
            velocity{d}(isinf(velocity{d}) | isnan(velocity{d}) | velocity{d}<log(1/8)) = [];
        for s = 1:nboot
            boot = ceil(length(velocity{d})*rand(length(velocity{d}),1));
            xboot = velocity{d}(boot);
            yboot = thresh{d}(boot);
            b = regress(yboot,[ones(length(xboot),1) xboot],0.05);
            Bboot(s,d) = b(2);
            B0boot(s,d) = b(1);
        end
    end
end

subs = repmat([1:length(velocity)],nboot,1);
subs = reshape(subs,nboot*length(velocity),1);
boot = reshape(Bboot,nboot*length(velocity),1);
novelboot = boot(subs<=10);
novelsubs = subs(subs<=10);
b = regress(novelboot,[ones(length(novelsubs),1) novelsubs],0.05);
qnboot = nboot;
permr = nan(qnboot,1);
for q = 1:qnboot
    perm = novelsubs(randperm(length(novelsubs)));
    tmpb = regress(novelboot,[ones(length(novelsubs),1) perm],0.05);
    permr(q) = tmpb(2);
end
pvalue = sum(b(2)>permr)./qnboot;
A = accumarray(novelsubs,novelboot,[max(novelsubs) 1],@(x) nanmean(x),NaN);
seL = A-accumarray(novelsubs,novelboot,[max(novelsubs) 1],@(x) prctile(x,[2.5]));
seU = A-accumarray(novelsubs,novelboot,[max(novelsubs) 1],@(x) prctile(x,[97.5]));

%--------------------------------------------------------------------------
% Plot slope over days
figure
hold on
plot(1:max(novelsubs),A,'ko','MarkerFace','k')
plot(1:max(novelsubs),A,'k')
errorbar2(1:max(novelsubs),A,[seL seU]',0.25,'k')
xlabel('Exposure to Track B','FontSize',24)
ylabel('Slope of Speed vs. Normalized Ripple Power','FontSize',24)
set(gca,'xtick',1:max(novelsubs),'xlim',[0 max(novelsubs)+1])
set(gca,'ylim',[-0.4 0],'ytick',[-0.4:0.1:0],'FontSize',20)
box off

%--------------------------------------------------------------------------
% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/SpeedPaper/%d_%d_%d_slopes_velvsripthresh_controlspeed_exposures1_10.pdf', m, d, y);
print('-dpdf', savestring)

%% Combine data across animals, control for unequal behavior
bin = log([1/4 1 4 16]);

count = [];
for an = 1:length(vel)
    for d = 1:length(vel(an).output)
    	if ~isempty(vel(an).output{d})
        	tmpv = log(vel(an).output{d}.velocity);
            tmpt = nanmean(rip(an).output{d}.maxthresh,2);
            tmpt(isinf(tmpv) | isnan(tmpv)) = [];
            tmpv(isinf(tmpv) | isnan(tmpv)) = [];
            subs = lookup(tmpv,bin);
        	count= [count; min(hist(subs,length(bin)))];
    	end
    end
end
count = min(count(count>10));

nboot = 1000;
Bboot = nan(nboot,14);
for s = 1:nboot
    velocity = cell(14,1);
    thresh = cell(14,1);
    for an =1:length(rip)
        for d = 1:length(rip(an).output)
            if ~isempty(rip(an).output{d})
                tmpv = log(vel(an).output{d}.velocity);
                tmpt = nanmean(rip(an).output{d}.maxthresh,2);
                tmpt(isinf(tmpv) | isnan(tmpv)) = [];
                tmpv(isinf(tmpv) | isnan(tmpv)) = [];
                if length(tmpv)>count
                    subs = lookup(tmpv,bin);
                    boot = ceil(length(tmpv)*rand(2*length(tmpv),1));
                    for q = 1:max(subs)
                        tmp = tmpv(boot);
                        tmp = tmp(subs == q);
                        tmp = tmp(1:count);
                        velocity{d} = [velocity{d}; tmp];
                        tmp = tmpt(boot);
                        tmp = tmp(subs == q);
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
[b,bint,r,rint,stats] = regress(novelboot,[ones(length(novelsubs),1) novelsubs],0.05);
A = accumarray(novelsubs,novelboot,[max(novelsubs) 1],@(x) nanmean(x),NaN);
seL = A-accumarray(novelsubs,novelboot,[max(novelsubs) 1],@(x) prctile(x,[2.5]));
seU = A-accumarray(novelsubs,novelboot,[max(novelsubs) 1],@(x) prctile(x,[97.5]));

%--------------------------------------------------------------------------
% Plot slope over days
figure
hold on
plot(1:max(novelsubs),A,'ko','MarkerFace','k')
errorbar2(1:max(novelsubs),A,[seL seU]',0.001,'k')
plot(1:max(novelsubs),[1:max(novelsubs)]*b(2)+b(1),'k')    
xlabel('Exposure to Track','FontSize',18)
ylabel('Slope of Speed vs. Normalized Ripple Power','FontSize',18)
set(gca,'xtick',1:max(novelsubs),'xlim',[0 max(novelsubs)+1])
set(gca,'ylim',[-0.7 0.1],'ytick',[-0.7:0.1:0],'FontSize',18)
box off
x = [1 max(novelsubs)];
LL = min(bint(2,1)*x(1)+bint(1,1), bint(2,1)*x(1)+bint(1,2));
LU = max(bint(2,2)*x(1)+bint(1,1), bint(2,2)*x(1)+bint(1,2));
RU = max(bint(2,2)*x(2)+bint(1,1), bint(2,2)*x(2)+bint(1,2));
RL = min(bint(2,1)*x(2)+bint(1,1), bint(2,1)*x(2)+bint(1,2));
fill([x(1) x(1) x(2) x(2)], [LL LU RU RL],'k','FaceAlpha',1,'EdgeColor','none')

% b: [-0.3092 0.0268] bint: [-0.3105 -0.3080; 0.0266 0.0270]
% stats:
% R = 0.6388 F = 68942 p = 0
% numpoints:100000
%% Look at the change in slopes over days for each animal
bin = log([1/4 1/2 1 2 4 8 16 32]);
pvalue = nan(length(rip),1);
for an =1:length(rip)
    
    velocity = cell(23,1);
    thresh = cell(23,1);
    for d = 1:length(rip(an).output)
        if ~isempty(rip(an).output{d})
            velocity{d} =[velocity{d}; log(vel(an).output{d}.velocity)];
            thresh{d} = [thresh{d}; nanmean(rip(an).output{d}.maxthresh,2)];
        end
    end

    nboot = 1000;
    Bboot = nan(nboot,length(velocity));
    B0boot = nan(nboot,length(velocity));
    for d = 1:length(velocity)
        if ~isempty(velocity{d})
                thresh{d}(isinf(velocity{d})) = [];
                velocity{d}(isinf(velocity{d})) = [];
            for s = 1:nboot
                boot = ceil(length(velocity{d})*rand(length(velocity{d}),1));
                xboot = velocity{d}(boot);
                yboot = thresh{d}(boot);
                b = regress(yboot,[ones(length(xboot),1) xboot],0.05);
                Bboot(s,d) = b(2);
                B0boot(s,d) = b(1);
            end
        end
    end

    subs = repmat([1:length(velocity)],nboot,1);
    subs = reshape(subs,nboot*length(velocity),1);
    boot = reshape(Bboot,nboot*length(velocity),1);
    novelboot = boot(subs<=10);
    novelsubs = subs(subs<=10);
    if sum(isnan(novelboot)) ~= length(novelboot)
        b = regress(novelboot,[ones(length(novelsubs),1) novelsubs],0.05);
        qnboot = nboot;
        permr = nan(qnboot,1);
        for q = 1:qnboot
            perm = novelsubs(randperm(length(novelsubs)));
            tmpb = regress(novelboot,[ones(length(novelsubs),1) perm],0.05);
            permr(q) = tmpb(2);
        end
        pvalue(an) = 1-sum(b(2)>permr)./qnboot;
    end
end
