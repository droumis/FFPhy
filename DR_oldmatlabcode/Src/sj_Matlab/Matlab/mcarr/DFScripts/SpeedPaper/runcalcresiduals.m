%Animal Selection
animals = {'Conley','Corriander','Eight','Five','Miles','Ten'};

% Epoch selection
epochfilter = [];

for i = 1:10
    epochfilter{i} = ['($exposure == ',num2str(i),') & isequal($description,''TrackA'')'];
end

% Time Filter
timefilter = {{'getriptimes', '($nripples == 0)', [], 'cellfilter', 'isequal($area, ''CA1'')'}};

% Cell Filter
cellfilter = 'isequal($area, ''CA1'') & $meanrate < 7';

iterator = 'multicellanal';

f = createfilter('animal',animals,'epochs',epochfilter,'cells',cellfilter,'excludetimefilter', timefilter, 'iterator', iterator);

binsize = 0.125;
min_spikes = 100;
minrate = 0;
speed_bin = [1/4 1 4 16];
f = setfilterfunction(f, 'calcresiduals', {'spikes','linpos','pos'},'bin',binsize,'minspikes',min_spikes,'minrate',minrate,'speed_bin',speed_bin);

f = runfilter(f);

%% Save Filter
save('/data13/mcarr/VelocityPaper/CA1residuals.mat','f')
save('/data13/mcarr/VelocityPaper/CA1residuals_minrate1.mat','f')

%% Load Filter
load('/data13/mcarr/VelocityPaper/CA1residuals.mat')
load('/data13/mcarr/VelocityPaper/CA1residuals_minrate1.mat')

%% Plot speed correlations
% For all times and all animals
color = ['b' 'g' 'y' 'r'];
speedcorr = [];
for an = 1:length(f)
    for d = 1:4%length(f(an).output)
        for e = 1:length(f(an).output{d})
            speedcorr = stack(speedcorr,f(an).output{d}(e).cell_pair_info(:,8:11));
        end
    end
end

figure
for s = 1:4
    [e x] = ecdf(speedcorr(:,s));
    stairs(x,e,color(s));
    hold on
end
set(gca,'xtick',[-1 -.5 0 0.5 1],'ylim',[0 1])
ylabel('Cumulative Probability','FontSize',24)
xlabel('Correlation of Residuals','FontSize',24)
legend('Speed = 1/4 cm/sec','Speed = 1 cm/sec','Speed = 4 cm/sec','Speed = 16 cm/sec')

% Use KSTEST2 to determine significance
pairind = nchoosek(1:4,2);
pvalue = nan(size(pairind,1),2);
for s = 1:size(pairind,1)
    [h pvalue(s,1) pvalue(s,2)] = kstest2(speedcorr(:,pairind(s,1)),speedcorr(:,pairind(s,2)));
end

%Pvalue for all pairs = 1E-6

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/SpeedPaper/%d_%d_%d_noveldaysresidual_cdf.pdf', m, d, y);
print('-dpdf', savestring)

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/SpeedPaper/%d_%d_%d_noveldaysresidual_cdf_binsize125.pdf', m, d, y);
print('-dpdf', savestring)


%% Plot coordination index over days

color = ['b' 'g' 'y' 'r'];
speedcorr = [];
day = [];
for an = 1:length(f)
    for d = 1:length(f(an).output)
        for e = 1:length(f(an).output{d})
            speedcorr = stack(speedcorr,abs(f(an).output{d}(e).cell_pair_info(:,8:11)));
            day = stack(day,d.*ones(size(f(an).output{d}(e).cell_pair_info(:,8:11),1),1));
        end 
    end
end
A1 = accumarray(day,speedcorr(:,1),[max(day) 1],@(x) nanmean(x));
E1 = accumarray(day,speedcorr(:,1),[max(day) 1],@(x) nanstd(x)./sqrt(sum(~isnan(x))-1));
A2 = accumarray(day,speedcorr(:,2),[max(day) 1],@(x) nanmean(x));
E2 = accumarray(day,speedcorr(:,2),[max(day) 1],@(x) nanstd(x)./sqrt(sum(~isnan(x))-1));
A3 = accumarray(day,speedcorr(:,3),[max(day) 1],@(x) nanmean(x));
E3 = accumarray(day,speedcorr(:,3),[max(day) 1],@(x) nanstd(x)./sqrt(sum(~isnan(x))-1));
A4 = accumarray(day,speedcorr(:,4),[max(day) 1],@(x) nanmean(x));
E4 = accumarray(day,speedcorr(:,4),[max(day) 1],@(x) nanstd(x)./sqrt(sum(~isnan(x))-1));

figure
plot(1:max(day),A1,color(1),1:max(day),A2,color(2),1:max(day),A3,color(3),1:max(day),A4,color(4))
legend('Speed = 1/4 cm/sec','Speed = 1 cm/sec','Speed = 4 cm/sec','Speed = 16 cm/sec')
hold on
fill([1:max(day) fliplr(1:max(day))],[A1+E1; flipud(A1 - E1)],color(1),'EdgeColor',color(1))
fill([1:max(day) fliplr(1:max(day))],[A2+E2; flipud(A2 - E2)],color(2),'EdgeColor',color(2))
fill([1:max(day) fliplr(1:max(day))],[A3+E3; flipud(A3 - E3)],color(3),'EdgeColor',color(3))
fill([1:max(day) fliplr(1:max(day))],[A4+E4; flipud(A4 - E4)],color(4),'EdgeColor',color(4))
set(gca,'xtick',1:10,'ylim',[0 0.4])
ylabel('Correlation of Residuals')
xlabel('Exposure')
set(gca,'xlim',[0 11])

% Does the residual correlation increase with learning
pvalue = nan(4,2);
for s = 1:4
    [pvalue(s,1) pvalue(s,2)]=corr(day(~isnan(speedcorr(:,s))),speedcorr(~isnan(speedcorr(:,s)),s),'type','Spearman','rows','all');
end

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/SpeedPaper/ResidualCorrelation_overdays.pdf', m, d, y);
print('-dpdf', savestring)

%% Sort by peak rate and plot
y = cell(10,1);
for an = 1:length(f)
    for d = 1:4%length(f(an).output)
        for e = 1:length(f(an).output{d})
            y{d} = stack(y{d},f(an).output{d}(e).cell_pair_info);
        end 
    end
end

ratebin = [5 10 15 20 25 30];

rate1=[]; rate2=[]; rate=[]; speed1=[]; speed2=[]; speed3=[]; speed4=[]; all=[];
for d = 1:4%10
    % Assign peak rate of each cell to a rate bin;
    rate1 = stack(rate1,lookup(y{d}(:,3),ratebin));
    rate2 = stack(rate2,lookup(y{d}(:,4),ratebin));
    rate = stack(rate, mean(y{d}(:,[3 4]),2));
    speed1 = stack(speed1,y{d}(:,8));
    speed2 = stack(speed2,y{d}(:,9));
    speed3 = stack(speed3,y{d}(:,10));
    speed4 = stack(speed4,y{d}(:,11));
    all = stack(all,y{d}(:,6));
end

% "Coordination index"
ci1 = abs(speed1);
ci2 = abs(speed2);
ci3 = abs(speed3);
ci4 = abs(speed4);

ci1_se = nanstd(ci1)./sqrt(sum(~isnan(ci1))-1);
ci2_se = nanstd(ci2)./sqrt(sum(~isnan(ci2))-1);
ci3_se = nanstd(ci3)./sqrt(sum(~isnan(ci3))-1);
ci4_se = nanstd(ci4)./sqrt(sum(~isnan(ci4))-1);

% Plot mean and standard error for all cell pairs
figure
bar(1,nanmean(ci1),color(1))
hold on
bar(2,nanmean(ci2),color(2))
bar(3,nanmean(ci3),color(3))
bar(4,nanmean(ci4),color(4))
errorbar2([1 2 3 4],nanmean([ci1 ci2 ci3 ci4]),[ci1_se ci2_se ci3_se ci4_se],0.1,'k')
set(gca,'xtick',[1 2 3 4])
ylabel('Coordination Index')

%Anova to test for significance
group = [ones(size(ci1)); 2*ones(size(ci2)); 3*ones(size(ci3)); 4*ones(size(ci4))];
x = [ci1; ci2; ci3; ci4];

[p,table,stats] = kruskalwallis(x,group);
[c] = multcompare(stats,'alpha',1E-10,'ctype','bonferroni');
% all of the groups are different from one another at p<1E-15

%Save Figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/SpeedPaper/ResidualCorrelation_coordinationindex_novelexposure.pdf', m, d, y);
print('-dpdf', savestring)

%Plot mean and standard error as a function of average peak rate for
%fastest and slowest speed times
color = ['r' 'y' 'g' 'b'];
rind = lookup(rate,ratebin);
figure
hold on
m = nan(length(ratebin),4);
e = nan(size(m));
for s = 1:length(ratebin)
    m(s,:) = [nanmean(ci1(rind==s)) nanmean(ci2(rind==s)) nanmean(ci3(rind==s)) nanmean(ci4(rind==s))];
    e(s,:) = [nanstd(ci1(rind==s))./sqrt(sum(~isnan(ci1(rind==s)))-1) nanstd(ci2(rind==s))./sqrt(sum(~isnan(ci2(rind==s)))-1)...
        nanstd(ci3(rind==s))./sqrt(sum(~isnan(ci3(rind==s)))-1) nanstd(ci4(rind==s))./sqrt(sum(~isnan(ci4(rind==s)))-1)];
end
plot(ratebin,m(:,1),color(1),ratebin,m(:,2),color(2),ratebin,m(:,3),color(3),ratebin,m(:,4),color(4))
legend('Speed = 1/4 cm/sec','Speed = 1 cm/sec','Speed = 4 cm/sec','Speed = 16 cm/sec')
legend('Location','SouthWest')
fill([ratebin fliplr(ratebin)],[m(:,1)+e(:,1); flipud(m(:,1) - e(:,1))],color(1),'EdgeColor',color(1))
fill([ratebin fliplr(ratebin)],[m(:,2)+e(:,2); flipud(m(:,2) - e(:,2))],color(2),'EdgeColor',color(2))
fill([ratebin fliplr(ratebin)],[m(:,3)+e(:,3); flipud(m(:,3) - e(:,3))],color(3),'EdgeColor',color(3))
fill([ratebin fliplr(ratebin)],[m(:,4)+e(:,4); flipud(m(:,4) - e(:,4))],color(4),'EdgeColor',color(4))

set(gca,'xtick',ratebin,'xticklabel',ratebin,'xlim',[0 max(ratebin)+5],'ylim',[0 0.4])
ylabel('Coordination index')
xlabel('Cell pair peak rate')

%Anova to test for significance
group_speed = [ones(size(ci1)); 2*ones(size(ci2)); 3*ones(size(ci3)); 4*ones(size(ci4))];
group_rate = [rind; rind; rind; rind];
x = [ci1; ci2; ci3; ci4];

[p,table,stats] = anovan(x,{group_speed group_rate},'model','full');
[c] = multcompare(stats,'alpha',0.01,'ctype','bonferroni','dim',1);

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/SpeedPaper/ResidualCorrelation_acrossrates.pdf', m, d, y);
print('-dpdf', savestring)

%% Speed correlation within cell pairs
y = cell(10,1);
x = cell(10,1);
for an = 1:length(f)
    for d = 1:length(f(an).output)
        for e = 1:length(f(an).output{d})
            y{d} = stack(y{d},f(an).output{d}(e).cell_pair_info);
        end 
    end
end

% Look at how individual cell pairs covary with speed
rho = cell(10,1);
for d = 1:10
    rho{d} = corr([1/4 1 4 16]',y{d}(:,8:11)','type','Spearman','rows','all');   
end
tmprho = [];
for d = 1:10
    tmprho = stack(tmprho,rho{d});
end
figure
bar(-1:0.05:1,hist(reshape(tmprho,size(tmprho,1)*size(tmprho,2),1),-1:0.05:1)./length(reshape(tmprho,size(tmprho,1)*size(tmprho,2),1)))

%% Look at the relationship between speed and residuals
resid = [];
for an = 1:length(f)
    for d = 1:length(f(an).output)
        for e = 1:length(f(an).output{d})
            for c = 1:size(f(an).output{d}(e).resid,1)
                valid = f(an).output{d}(e).resid(c,:) ~= 0;
                if sum(valid) >0                    
                    [tmp tmp_p] = corr((f(an).output{d}(e).speed(valid)),f(an).output{d}(e).resid(c,valid)','type','Spearman');
                    resid = [resid; tmp' tmp_p'];
                end
            end
        end
    end
end
bar(-0.5:0.05:0.5,hist(resid(:,1),[-0.5:0.05:0.5])./length(resid))
set(gca,'xlim',[-0.5 0.5])
xlabel('Spearman correlation between speed and residuals')
ylabel('Proportion of cells')

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/SpeedPaper/%d_%d_%d_covariation_residuals_speed.pdf', m, d, y);
print('-dpdf', savestring)

%% Look at the relationship between overlap and residual correlation
overlap = [];
rho = [];
pvalue = [];
for an = 1:length(f)
    for d = 1:length(f(an).output)
        for e = 1:length(f(an).output{d})
            if ~isempty(f(an).output{d}(e).cell_pair_info)
                overlap = stack(overlap,f(an).output{d}(e).cell_pair_info(:,5));
                rho = stack(rho,f(an).output{d}(e).cell_pair_info(:,6));
                pvalue = stack(pvalue,f(an).output{d}(e).cell_pair_info(:,7));
            end
        end
    end
end
corr(overlap(pvalue<0.05),rho(pvalue<0.05),'type','Spearman','rows','complete')
corr(overlap,rho,'type','Spearman','rows','complete')

%% LOOK AT FIRST AND LAST DAY
y = cell(2,1);
for an = 1:length(f)
    for d = [1 length(f(an).output)]
        for e = 1:length(f(an).output{d})
            if d ==1
                ind = 1;
            else
                ind = 2;
            end
            y{ind} = stack(y{ind},f(an).output{d}(e).cell_pair_info);
        end 
    end
end

novel = abs(y{1}(:,8:11));
familiar = abs(y{2}(:,8:11));

corrN = nanmean(novel);
corrF = nanmean(familiar);
errN = nanstd(novel)./sqrt(sum(~isnan(novel))+1);
errF = nanstd(familiar)./sqrt(sum(~isnan(familiar))+1);


figure
hold on
bar([1 3 5 7],corrN,0.4,'r')
bar([2 4 6 8],corrF,0.4)
legend('Novel','Familiar')
errorbar2([1 3 5 7],corrN,errN,0.5,'k')
errorbar2([2 4 6 8],corrF,errF,0.5,'k')
set(gca,'xtick',[1.5 3.5 5.5 7.5],'xticklabel',{'1/4' '1' '4' '16'})
xlabel('Speed')
ylabel('Coordination Index')


% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/SpeedPaper/%d_%d_%d_accrossdayNF_residcoordinationindex.pdf', m, d, y);
print('-dpdf', savestring)

%% RUN WITHIN DAY COMPARISONS
binsize = 0.5;
min_spikes = 100;
minrate = 0;
speed_bin = [1/4 1 4 16];

%animal selection
animals = {'Corriander','Eight','Five','Miles','Ten'};

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

%% Plot speed correlations for novel and familiar track
% For all times and all animals
speedcorrA = []; speedcorrB = [];
for an = 1:length(A)
    for d = 1:length(A(an).output)
        for e = 1:length(A(an).output{d})
            speedcorrA = stack(speedcorrA,A(an).output{d}(e).cell_pair_info(:,8:11));
            speedcorrB = stack(speedcorrB,B(an).output{d}(e).cell_pair_info(:,8:11));
        end
    end
end

corrB = nanmean(abs(speedcorrB));
corrA = nanmean(abs(speedcorrA));
%errB = prctile(abs(speedcorrB),[25 75]) - [corrB; corrB];
%errA = prctile(abs(speedcorrA),[25 75]) - [corrA; corrA];
errB = nanstd(abs(speedcorrB))./sqrt(sum(~isnan(speedcorrB))-1);
errA = nanstd(abs(speedcorrA))./sqrt(sum(~isnan(speedcorrA))-1);


figure
hold on
bar([1 3 5 7],corrB,0.4,'r')
bar([2 4 6 8],corrA,0.4)
legend('Novel','Familiar')
errorbar2([1 3 5 7],corrB,errB,0.5,'k')
errorbar2([2 4 6 8],corrA,errA,0.5,'k')
set(gca,'xtick',[1.5 3.5 5.5 7.5],'xticklabel',{'1/4' '1' '4' '16'})
xlabel('Speed')
ylabel('Coordination Index')

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/SpeedPaper/%d_%d_%d_withindayNF_residcoordinationindex.pdf', m, d, y);
print('-dpdf', savestring)
