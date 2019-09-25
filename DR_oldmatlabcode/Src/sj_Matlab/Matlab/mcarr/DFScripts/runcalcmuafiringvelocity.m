%% Run filter and look at relationship between speed and firing rate

animals = {'Conley','Corriander','Dudley','Eight','Five','Miles','Ten'};

% Epoch selection
epochfilter = [];
for i = 1:10
    epochfilter{i} = ['($exposure == ',num2str(i),') & $dailytracksexperienced == 1 & isequal($description,''TrackA'')'];
end

% Cell selection
ca1cellfilter = '(isequal($area, ''CA3'') & ($meanrate < 7) & ($peakrate<5))';
 
%Select iterator
iterator = 'multicellanal';

% Create and Run Filter

f = createfilter('animal',animals,'epochs',epochfilter,'cells',ca1cellfilter,'iterator', iterator);
f = setfilterfunction(f, 'calcspeedfiringrate', {'spikes','pos'});
f = runfilter(f);
firingCA3 = f;

save('/data13/mcarr/VelocityPaper/firingCA3.mat','firingCA3')


%% LOAD DATA
load('/data13/mcarr/VelocityPaper/firingCA1.mat','firingCA1')
load('/data13/mcarr/VelocityPaper/firingCA3.mat','firingCA3')

%% PLOT EXPOSURE ONE
bin = log([1/4 1/2 1 2 4 8 16 32]);
vel1 = []; vel3 = [];
mua1 = []; mua3 = [];
for an = 1:length(firingCA1)
    if ~isempty(firingCA1(an).output)
        for d = 1
            vel1 = [vel1; log(firingCA1(an).output{d}.speed)];
            mua1 = [mua1; firingCA1(an).output{d}.firing./firingCA1(an).output{d}.numcells];
        end
    end
    if ~isempty(firingCA3(an).output)
        for d = 1
            vel3 = [vel3; log(firingCA3(an).output{d}.speed)];
            mua3 = [mua3; firingCA3(an).output{d}.firing./firingCA3(an).output{d}.numcells];
        end
    end
end

invalid = (vel1 < log(1/8)) | isinf(vel1);
mua1(invalid) = []; vel1(invalid) = [];
invalid = (vel3 < log(1/8)) | isinf(vel3);
mua3(invalid) = []; vel3(invalid) = [];

subs = lookup(vel1,bin);
A1 = accumarray(subs,mua1,[length(bin) 1],@(x) nanmean(x),NaN);
Lstd = accumarray(subs,mua1,[length(bin) 1],@(x) nanstd(x),NaN);
Lsqrt = accumarray(subs,mua1,[length(bin) 1],@(x) sqrt(sum(~isnan(x))),NaN);
seL1 = -Lstd./Lsqrt;
seU1 = +Lstd./Lsqrt;

subs = lookup(vel3,bin);
A3 = accumarray(subs,mua3,[length(bin) 1],@(x) nanmean(x),NaN);
Lstd = accumarray(subs,mua3,[length(bin) 1],@(x) nanstd(x),NaN);
Lsqrt = accumarray(subs,mua3,[length(bin) 1],@(x) sqrt(sum(~isnan(x))),NaN);
seL3 = -Lstd./Lsqrt;
seU3 = +Lstd./Lsqrt;

[L1 U1 slope1 intercept1] = regress_fill(vel1,mua1,bin);
[L3 U3 slope3 intercept3] = regress_fill(vel3,mua3,bin);

figure
fill([bin fliplr(bin)], [L1 fliplr(U1)],'r','FaceAlpha',1,'EdgeColor','none')
hold on
fill([bin fliplr(bin)], [L3 fliplr(U3)],'k','FaceAlpha',1,'EdgeColor','none')
legend('CA1','CA3','Location','NorthWest')

plot(bin,bin*slope1+intercept1,'r')    
plot(bin,A1,'ro','MarkerFace','r')
plot(bin,bin*slope3+intercept3,'k')    
plot(bin,A3,'ko','MarkerFace','k')

errorbar2(bin,A1,[seL1 seU1]',0.001,'r','plottype','semilogx')
errorbar2(bin,A3,[seL3 seU3]',0.001,'k','plottype','semilogx')

set(gca,'xtick',bin,'xtickLabel',exp(bin),'yLim',[0 1])
set(gca,'FontSize',18)
xlabel('Speed (cm/sec)','FontSize',18)
ylabel('Multiunit Firing Rate','FontSize',18)
box off

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/SpeedPaper/%d_%d_%d_muafiring_day1.pdf', m, d, y);
print('-dpdf', savestring)

%% PLOT SLOPE ACROSS EXPOSURES

%bin = log([1/4 1/2 1 2 4 8 16 32]);
bin = log([1/4 1 4 16]);
vel1 = cell(10,1); vel3 = cell(10,1);
mua1 = cell(10,1); mua3 = cell(10,1);
for an = 1:length(firingCA1)
    if ~isempty(firingCA1(an).output)
        for d = 1:10
            if length(firingCA1(an).output)>=d
                if firingCA1(an).output{d}.numcells > 3
                    vel1{d} = [vel1{d}; log(firingCA1(an).output{d}.speed)];
                    mua1{d} = [mua1{d}; firingCA1(an).output{d}.firing./firingCA1(an).output{d}.numcells];
                    invalid = (vel1{d} < log(1/8)) | isinf(vel1{d});
                    mua1{d}(invalid) = []; vel1{d}(invalid) = [];
                end
            end
        end
    end
    if ~isempty(firingCA3(an).output)
        for d = 1:10
            if length(firingCA3(an).output)>=d
                if ~isempty(firingCA3(an).output{d})
                    if firingCA3(an).output{d}.numcells > 3
                        vel3{d} = [vel3{d}; log(firingCA3(an).output{d}.speed)];
                        mua3{d} = [mua3{d}; firingCA3(an).output{d}.firing./firingCA3(an).output{d}.numcells];
                        invalid = (vel3{d} < log(1/8)) | isinf(vel3{d});
                        mua3{d}(invalid) = []; vel3{d}(invalid) = [];
                    end
                end
            end
        end
    end
end

nboot = 500;
Bboot1 = nan(nboot,10);
Bboot3 = nan(nboot,10);
count = 50;


for d = 1:10
   
    subs1 = lookup(vel1{d},bin);
    subs3 = lookup(vel3{d},bin);

    for s = 1:nboot
        speed1_boot = []; speed3_boot = []; mua1_boot = []; mua3_boot = [];
        boot1 = ceil(length(subs1)*rand(2*length(subs1),1));
        boot3 = ceil(length(subs3)*rand(2*length(subs3),1));
        for q = 1:max([subs1; subs3])
            tmpboot = boot1(subs1(boot1)==q);
            tmp = vel1{d}(tmpboot);
            if ~isempty(tmp)
                tmp = tmp(1:count);
            end
            speed1_boot = [speed1_boot; tmp];
            tmp = mua1{d}(tmpboot);
            if ~isempty(tmp)
                tmp = tmp(1:count);
            end
            mua1_boot = [mua1_boot; tmp];
        
            tmpboot = boot3(subs3(boot3)==q);
            tmp = vel3{d}(tmpboot);
            if ~isempty(tmp)
                tmp = tmp(1:count);
            end
            speed3_boot = [speed3_boot; tmp];
            tmp = mua3{d}(tmpboot);
            if ~isempty(tmp)
                tmp = tmp(1:count);
            end
            mua3_boot = [mua3_boot; tmp];
        end
        
        if ~isempty(mua1_boot)
            b = regress(mua1_boot,[ones(length(speed1_boot),1) speed1_boot]);
            Bboot1(s,d) = b(2);
        end
        if ~isempty(mua3_boot)
            b = regress(mua3_boot,[ones(length(speed3_boot),1) speed3_boot]);
            Bboot3(s,d) = b(2);
        end
    end
end

figure
plot(Bboot1','r.')
hold on
plot(Bboot3','k.')


% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/SpeedPaper/%d_%d_%d_muafiring_slope_acrossday.pdf', m, d, y);
print('-dpdf', savestring)

