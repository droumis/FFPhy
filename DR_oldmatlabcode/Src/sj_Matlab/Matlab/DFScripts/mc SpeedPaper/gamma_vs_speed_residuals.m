%% Determine the relationship between gamma and residuals as well as speed 
%  and residual correlation.

load('/data13/mcarr/VelocityPaper/CA1residuals.mat')
load('/data13/mcarr/VelocityPaper/power.mat')

speed = []; low = []; high = [];
for an = 1:length(f)
    for d = 1:3%length(f(an).output)
        lowgam = []; highgam = [];
        time = power(an).output{d}.time;
        cell_time = f(an).output{d}.time + median(diff(f(an).output{d}.time))/2;
        
        index = lookup(time,cell_time);
        
        for tet = 1:length(power(an).output{d})
            lowgam = [lowgam power(an).output{d}(tet).lowgamma_power];
            highgam = [highgam power(an).output{d}(tet).highgamma_power];
        end
        lowgam = (mean(lowgam,2) - mean(mean(lowgam,2)))./std(mean(lowgam,2));
        highgam = (mean(highgam,2) - mean(mean(highgam,2)))./std(mean(highgam,2));
        
        lowind = lowgam > highgam & lowgam>1;
        highind = highgam > lowgam & highgam>1;
        if size(f(an).output{d}.resid,1) >=2
            pairind = nchoosek(1:size(f(an).output{d}.resid,1),2);
            for cell = 1:size(pairind,1)
                resid1 = f(an).output{d}.resid(pairind(cell,1),index);
                resid2 = f(an).output{d}.resid(pairind(cell,2),index);
                nonzero = (resid1 & resid2)';
                if sum(~isnan(resid1(nonzero & lowind)) & ~isnan(resid2(nonzero & lowind)))*0.5 > 10
                    tmp = corr(resid1(nonzero & lowind)',resid2(nonzero & lowind)','type','Spearman','rows','complete');
                    low = [low tmp];
                end
                if sum(~isnan(resid1(nonzero & highind)) & ~isnan(resid2(nonzero & highind)))*0.5 > 10
                    tmp = corr(resid1(nonzero & highind)',resid2(nonzero & highind)','type','Spearman','rows','complete');
                    high = [high tmp];
                end
                speed = [speed; f(an).output{d}.cell_pair_info(cell,8:11)];
            end
        end
    end
end

%% Plot CDFs
color = ['b' 'g' 'y' 'r'];

figure
for s = 1:4
    [e x] = ecdf(speed(:,s));
    stairs(x,e,color(s));
    hold on
end
[e x] = ecdf(low);
stairs(x,e,'c')
[e x] = ecdf(high);
stairs(x,e,'m')

% Use KSTEST2 to determine significance
pvalue_low = nan(4,2);
pvalue_high = nan(4,2);
for s = 1:4
    [h pvalue_low(s,1) pvalue_low(s,2)] = kstest2(low,speed(:,s));
    [h pvalue_high(s,1) pvalue_high(s,2)] = kstest2(high,speed(:,s));
end

% Plot Coordination index over days

ci1 = abs(speed(:,1));
ci2 = abs(speed(:,2));
ci3 = abs(speed(:,3));
ci4 = abs(speed(:,4));
l = abs(low)';
h = abs(high)';

l_se = prctile(l,[25 75]);
h_se = prctile(h,[25 75]);

ci1_se = prctile(ci1,[25 75]);
ci2_se = prctile(ci2,[25 75]);
ci3_se = prctile(ci3,[25 75]);
ci4_se = prctile(ci4,[25 75]);

% Plot mean and standard error for all cell pairs
figure
bar(1,nanmedian(ci1),color(1))
hold on
bar(2,nanmedian(ci2),color(2))
bar(3,nanmedian(ci3),color(3))
bar(4,nanmedian(ci4),color(4))
bar(5,nanmedian(l),'b')
bar(6,nanmedian(h),'r')
errorbar2([1 2 3 4 5 6],[nanmedian([ci1 ci2 ci3 ci4]) nanmedian(l) nanmedian(h)],...
    [ci1_se(1) ci2_se(1) ci3_se(1) ci4_se(1) l_se(1) h_se(1)],[ci1_se(2) ci2_se(2) ci3_se(2) ci4_se(2) l_se(2) h_se(2)],'k')
set(gca,'xtick',[1 2 3 4 5 6])
ylabel('Coordination Index','FontSize',24)
xlabel('Speed (cm/sec)','FontSize',24)
set(gca,'xtick',[1 2 3 4 5 6],'xticklabel',[{'1/4'} {'1'} {'4'} {'16'} {'Slow Gamma'} {'Fast Gamma'}],'FontSize',18)
%Kruskalwallis anova to test for significance
group = [ones(size(ci1)); 2*ones(size(ci2)); 3*ones(size(ci3)); 4*ones(size(ci4)); 5*ones(size(l)); 6*ones(size(h))];
x = [ci1; ci2; ci3; ci4; l; h];

[p,table,stats] = kruskalwallis(x,group);
[c] = multcompare(stats,'alpha',0.0001,'ctype','bonferroni');

%% RUN WITHIN DAY COMPARISONS

binsize = 0.5;
min_spikes = 100;
minrate = 0;
speed_bin = [1/4 1 4 16];

%animal selection
animals = {'Eight','Five','Miles','Ten'};

% Epoch selection
epochfilterA = []; epochfilterB = [];
epochfilterA{1} = ['($experimentday == 4) & $dailyexposure == 1 & isequal($description,''TrackA'')'];
epochfilterB{1} = ['($experimentday == 4) & $dailyexposure == 1 & isequal($description,''TrackB'')'];

% Time Filter
timefilter = {{'getriptimes', '($nripples ==0)', [], 'cellfilter', 'isequal($area, ''CA1'')'}};

% Cell Filter
cellfilter = 'isequal($area, ''CA1'') & $meanrate < 7';

iterator = 'multicellanal';

A = createfilter('animal',animals,'epochs',epochfilterA,'cells',cellfilter,'excludetimefilter', timefilter, 'iterator', iterator);
A = setfilterfunction(A, 'calcresiduals', {'spikes','linpos','pos'},'bin',binsize,'minspikes',min_spikes,'minrate',minrate,'speed_bin',speed_bin);
A = runfilter(A);

B = createfilter('animal',animals,'epochs',epochfilterB,'cells',cellfilter,'excludetimefilter', timefilter, 'iterator', iterator);
B = setfilterfunction(B, 'calcresiduals', {'spikes','linpos','pos'},'bin',binsize,'minspikes',min_spikes,'minrate',minrate,'speed_bin',speed_bin);
B = runfilter(B);

% LOAD FILTERS
load '/data13/mcarr/VelocityPaper/withinday_eegA.mat'
load '/data13/mcarr/VelocityPaper/withinday_eegB.mat'

eegA = withinday_eegA; clear withinday_eegA
eegB = withinday_eegB; clear withinday_eegB


speedF = []; lowF = []; highF = []; speedN = []; lowN = []; highN = [];
f = A;
power = eegA;

for an = 1:length(f)
    for d = 1:length(f(an).output)
        lowgam = []; highgam = [];
        time = power(an).output{d}.time;
        cell_time = f(an).output{d}.time + median(diff(f(an).output{d}.time))/2;
        
        index = lookup(time,cell_time);
        
        for tet = 1:length(power(an).output{d})
            lowgam = [lowgam power(an).output{d}(tet).lowgamma_power];
            highgam = [highgam power(an).output{d}(tet).highgamma_power];
        end
        lowgam = (mean(lowgam,2) - mean(mean(lowgam,2)))./std(mean(lowgam,2));
        highgam = (mean(highgam,2) - mean(mean(highgam,2)))./std(mean(highgam,2));
        
        lowind = lowgam > highgam & lowgam>1;
        highind = highgam > lowgam & highgam>1;
        if size(f(an).output{d}.resid,1) >=2
            pairind = nchoosek(1:size(f(an).output{d}.resid,1),2);
            for cell = 1:size(pairind,1)
                resid1 = f(an).output{d}.resid(pairind(cell,1),index);
                resid2 = f(an).output{d}.resid(pairind(cell,2),index);
                nonzero = (resid1 & resid2)';
                if sum(~isnan(resid1(nonzero & lowind)) & ~isnan(resid2(nonzero & lowind)))*0.5 > 10
                    tmp = corr(resid1(nonzero & lowind)',resid2(nonzero & lowind)','type','Spearman','rows','complete');
                    lowF = [lowF tmp];
                end
                if sum(~isnan(resid1(nonzero & highind)) & ~isnan(resid2(nonzero & highind)))*0.5 > 10
                    tmp = corr(resid1(nonzero & highind)',resid2(nonzero & highind)','type','Spearman','rows','complete');
                    highF = [highF tmp];
                end
                speedF = [speedF; f(an).output{d}.cell_pair_info(cell,8:11)];
            end
        end
    end
end

f = B;
power = eegB;

for an = 1:length(f)
    for d = 1:length(f(an).output)
        lowgam = []; highgam = [];
        time = power(an).output{d}.time;
        cell_time = f(an).output{d}.time + median(diff(f(an).output{d}.time))/2;
        
        index = lookup(time,cell_time);
        
        for tet = 1:length(power(an).output{d})
            lowgam = [lowgam power(an).output{d}(tet).lowgamma_power];
            highgam = [highgam power(an).output{d}(tet).highgamma_power];
        end
        lowgam = (mean(lowgam,2) - mean(mean(lowgam,2)))./std(mean(lowgam,2));
        highgam = (mean(highgam,2) - mean(mean(highgam,2)))./std(mean(highgam,2));
        
        lowind = lowgam > highgam & lowgam>1;
        highind = highgam > lowgam & highgam>1;
        if size(f(an).output{d}.resid,1) >=2
            pairind = nchoosek(1:size(f(an).output{d}.resid,1),2);
            for cell = 1:size(pairind,1)
                resid1 = f(an).output{d}.resid(pairind(cell,1),index);
                resid2 = f(an).output{d}.resid(pairind(cell,2),index);
                nonzero = (resid1 & resid2)';
                if sum(~isnan(resid1(nonzero & lowind)) & ~isnan(resid2(nonzero & lowind)))*0.5 > 10
                    tmp = corr(resid1(nonzero & lowind)',resid2(nonzero & lowind)','type','Spearman','rows','complete');
                    lowN = [lowN tmp];
                end
                if sum(~isnan(resid1(nonzero & highind)) & ~isnan(resid2(nonzero & highind)))*0.5 > 10
                    tmp = corr(resid1(nonzero & highind)',resid2(nonzero & highind)','type','Spearman','rows','complete');
                    highN = [highN tmp];
                end
                speedN = [speedN; f(an).output{d}.cell_pair_info(cell,8:11)];
            end
        end
    end
end

corrN = [nanmean(abs(speedN)) nanmean(abs(lowN)) nanmean(abs(highN))];
corrF = [nanmean(abs(speedF)) nanmean(abs(lowF)) nanmean(abs(highF))];
errN = [nanstd(abs(speedN))./sqrt(sum(~isnan(speedN))-1)...
    nanstd(abs(lowN))./sqrt(sum(~isnan(lowN))-1)...
    nanstd(abs(highN))./sqrt(sum(~isnan(highN))-1)];
errF = [nanstd(abs(speedF))./sqrt(sum(~isnan(speedF))-1)...
    nanstd(abs(lowF))./sqrt(sum(~isnan(lowF))-1)...
    nanstd(abs(highF))./sqrt(sum(~isnan(highF))-1)];


figure
hold on
bar([1 3 5 7 10 12],corrN,0.4,'r')
bar([2 4 6 8 11 13],corrF,0.4)
legend('Novel','Familiar')
errorbar2([1 3 5 7 10 12],corrN,errN,0.5,'k')
errorbar2([2 4 6 8 11 13],corrF,errF,0.5,'k')
set(gca,'xtick',[1.5 3.5 5.5 7.5 10.5],'xticklabel',{'1/4' '1' '4' '16' 'Gamma'})
xlabel('Speed')
ylabel('Coordination Index')

%Just Plot Gamam
figure
hold on
bar([1 3],corrN(:,5:6),0.8,'r')
bar([2 4],corrF(:,5:6),0.8)
legend('Novel','Familiar')
errorbar2([1 3],corrN(:,5:6),errN(:,5:6),0.5,'k')
errorbar2([2 4],corrF(:,5:6),errF(:,5:6),0.5,'k')
set(gca,'xtick',[1.5 3.5],'xticklabel',{'Slow Gamma' 'Fast Gamma'},'ylim', [0 0.4])
xlabel('Speed')
ylabel('Coordination Index')

%Save Figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/SpeedPaper/%d_%d_%d_residualCorrelation_gamma.pdf', m, d, y);
print('-dpdf', savestring)
