%% RIPPLE PAPER FIGURES: POWER

%Animals
animals = {'Bond','Frank','Ten'};

%Define epochs
epochfilter = [];
epochfilter{1} = 'isequal($type, ''run'')';
epochfilter{2} = 'isequal($type, ''sleep'')';
%Define iterator
iterator = 'epochbehaveanal';

%Define time filter
timefilter = {{'get2dstate', '($velocity < 4)'}};

f = createfilter('animal',animals,'epochs',epochfilter,'excludetime',timefilter,'iterator', iterator);
f = setfilterfunction(f, 'calchighgammapower', {'rips'});
f = runfilter(f);

% save('/data21/mcarr/RipplePaper/highgammapower.mat','f')

%% PLOT HIGH GAMMA VS. LOW GAMMA
time = f(1).output{1}(1).time;
ca1 =[]; ca3 = []; ca1g = []; ca3g = [];
for an = 1:length(f)
    for d = 1%:length(f(an).output)
        for e = 1:length(f(an).output{d})
            ca1g = [ca1g f(an).output{d}(e).ca1_lowgamma];
            ca1 = [ca1 f(an).output{d}(e).ca1_highgamma];
            ca3g = [ca3g f(an).output{d}(e).ca3_lowgamma];
            ca3 = [ca3 f(an).output{d}(e).ca3_highgamma];
        end
    end
end


mean_1 = mean(ca1,2); se_1 = std(ca1,[],2)./sqrt(size(ca1,2)-1);
mean_1g = mean(ca1g,2); se_1g = std(ca1g,[],2)./sqrt(size(ca1g,2)-1);
mean_3 = mean(ca3,2); se_3 = std(ca3,[],2)./sqrt(size(ca3,2)-1);
mean_3g = mean(ca3g,2); se_3g = std(ca3g,[],2)./sqrt(size(ca3g,2)-1);

x = -0.4:0.1:0.4;

figure
hold on
bar(time,mean_1g,'b')
bar(time,mean_1,'c')
errorbar2(time,mean_1,se_1,'k')
errorbar2(time,mean_1g,se_1g,'k')
legend([{'CA1 low gamma'},{'CA1 high gamma'}],'location','NorthWest')
set(gca,'xtick',time,'xticklabel',x,'ylim',[-0.5 1.5],'xlim',[-0.5 0.5])
ylabel('Normalized gamma power')
xlabel('Time since ripple detection (s)')
box off

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_S2.pdf', m, d, y);
print('-dpdf', savestring)


figure
hold on
bar(time,mean_3g,'r')
bar(time,mean_3,'m')
errorbar2(time,mean_3g,se_3g,'k')
errorbar2(time,mean_3,se_3,'k')
legend([{'CA3 low gamma'},{'CA3 high gamma'}],'location','NorthWest')
set(gca,'xtick',time,'xticklabel',x,'ylim',[-0.5 1.5],'xlim',[-0.5 0.5])
ylabel('Normalized gamma power')
xlabel('Time since ripple detection (s)')
box off

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_S2_right.pdf', m, d, y);
print('-dpdf', savestring)

%Test for significance.
group_time = [ones(size(ca1)) ones(size(ca1g))]; group_gamma = [ones(size(ca1)) 2*ones(size(ca1g))];
for i = 1:size(ca1,1)
    group_time(i,:) =i;
end

[p table stats] = anovan(reshape([ca1 ca1g],size(ca1,1)*2*size(ca1,2),1),{reshape(group_time,size(group_time,1)*size(group_time,2),1) reshape(group_gamma,size(group_gamma,1)*size(group_gamma,2),1)},'model','full');
c = multcompare(stats);
%Main effect of time and gamma, interaction between two with slow gamma
%significantly larger than fast gamma, p<1e-5

[p table stats] = kruskalwallis(reshape(ca3g-ca3,size(ca1,1)*size(ca1,2),1),reshape(group_time,size(group_time,1)*size(group_time,2),1));
c = multcompare(stats);
%Main effect of time and gamma, interaction between two with slow gamma
%significantly larger than fast gamma, p<1e-5

%% PLOT LOW RIPPLE POWER

time = f(1).output{1}(1).time;
ca1 =[]; ca3 = [];
for an = 1:length(f)
    for d = 1%:length(f(an).output)
        for e = 1:length(f(an).output{d})
            ca1 = [ca1 f(an).output{d}(e).ca1_lowrip];
            ca3 = [ca3 f(an).output{d}(e).ca3_lowrip];
        end
    end
end


mean_1 = mean(ca1,2); se_1 = std(ca1,[],2)./sqrt(size(ca1,2)-1);
mean_3 = mean(ca3,2); se_3 = std(ca3,[],2)./sqrt(size(ca3,2)-1);

x = -0.4:0.1:0.4;

figure
hold on
bar(time,mean_1,'c')
bar(time,mean_3,'m')
errorbar2(time,mean_1,se_1,'k')
errorbar2(time,mean_3,se_3,'k')
legend([{'CA1 lowrip'},{'CA3 low rip'}],'location','NorthWest')
set(gca,'xtick',time,'xticklabel',x,'ylim',[-0.5 1.5],'xlim',[-0.5 0.5])
ylabel('Normalized gamma power')
xlabel('Time since ripple detection (s)')
box off

%Test for significance.
group = ones(size(ca1));
for i = 1:size(ca1,1)
    group(i,:) =i;
end

[p table stats] = kruskalwallis(reshape(ca1,size(ca1,1)*size(ca1,2),1),reshape(group,size(group,1)*size(group,2),1));
c = multcompare(stats);

[p table stats] = kruskalwallis(reshape(ca3,size(ca3,1)*size(ca3,2),1),reshape(group,size(group,1)*size(group,2),1));
c = multcompare(stats);

%Both CA1 and CA3 show a significant increase in low ripple power

%% RIPPLE PAPER FIGURES: COHERENCE

%Animals
animals = {'Bond','Frank','Ten'};

%Define epochs
epochfilter = [];
epochfilter{1} = 'isequal($type, ''run'')';

%Define iterator
iterator = 'epochbehaveanal';

%Define time filter
timefilter = {{'get2dstate', '($velocity < 4)'}};

f = createfilter('animal',animals,'epochs',epochfilter,'excludetime',timefilter,'iterator', iterator);
f = setfilterfunction(f, 'calchighgammacoherence', {'ripc'});
f = runfilter(f);

%% PLOT COHERENCE
time = f(1).output{1}(1).time;
gam = []; rip =[]; lowrip = []; baseline_gam = []; baseline_rip = []; baseline_lowrip = [];
for an = 1:length(f)
    for d = 1:length(f(an).output)
        for e = 1:length(f(an).output{d})
            gam = [gam f(an).output{d}(e).ca1_ca3];
            baseline_gam = [baseline_gam f(an).output{d}(e).ca1_ca3_baseline];
            rip = [rip f(an).output{d}(e).ca1_ca3_ripple];
            baseline_rip = [baseline_rip f(an).output{d}(e).ca1_ca3_ripple_baseline];
            lowrip = [lowrip f(an).output{d}(e).ca1_ca3_lowripple];
            baseline_lowrip = [baseline_lowrip f(an).output{d}(e).ca1_ca3_lowripple_baseline];
        end
    end
end

mean_gam = mean(gam-repmat(baseline_gam,length(time),1),2);
se_gam = std((gam-repmat(baseline_gam,length(time),1))')./sqrt(length(gam)-1);
mean_rip = mean(rip-repmat(baseline_rip,length(time),1),2);
se_rip = std((rip-repmat(baseline_rip,length(time),1))')./sqrt(length(gam)-1);
mean_lowrip = mean(lowrip-repmat(baseline_lowrip,length(time),1),2);
se_lowrip = std((lowrip-repmat(baseline_lowrip,length(time),1))')./sqrt(length(gam)-1);
x = -0.4:0.1:0.4;

figure
hold on
bar(time,mean_gam,'b')
errorbar2(time,mean_gam,se_gam,'k')
ylabel('Delta CA1-CA3 high gamma coherence')
set(gca,'xtick',time,'xticklabel',x,'ylim',[-0.01 0.01],'ytick',-0.01:0.005:0.01)
xlabel('Time since SWR detection (s)')
box off

%Test significance using kruskalwallis anova
group = ones(size(gam)); 
for i = 1:size(gam,1)
    group(i,:) =i;
end
[p table stats] = kruskalwallis(reshape(gam,size(gam,1)*size(gam,2),1)-repmat(baseline_gam',size(gam,1),1),reshape(group,size(group,1)*size(group,2),1));
c = multcompare(stats);
%RUN:
% No significant differences from baseline
% baseline = 0.54


% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_ca1ca3_highgamma_coherence.pdf', m, d, y);
print('-dpdf', savestring)

figure
hold on
bar(time,mean_rip,'c')
errorbar2(time,mean_rip,se_rip,'k')
ylabel('Delta CA1-CA3 ripple coherence')
set(gca,'xtick',time,'xticklabel',x,'ylim',[-0.025 0.025],'ytick',-0.05:0.005:0.05)
xlabel('Time since SWR detection (s)')
box off

%Test significance using kruskalwallis anova
group = ones(size(rip)); 
for i = 1:size(rip,1)
    group(i,:) =i;
end
[p table stats] = kruskalwallis(reshape(rip,size(rip,1)*size(rip,2),1)-repmat(baseline_rip',size(rip,1),1),reshape(group,size(group,1)*size(group,2),1));
c = multcompare(stats);
%RUN:
% Time -200,-100ms < baseline at p<0.05
% Time 0,100ms < baseline at p<1e-5
% baseline = 0.51

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_ca1ca3_ripple_coherence.pdf', m, d, y);
print('-dpdf', savestring)


figure
hold on
bar(time,mean_lowrip,'b')
errorbar2(time,mean_lowrip,se_lowrip,'k')
ylabel('Delta CA1-CA3 low ripple coherence')
set(gca,'xtick',time,'xticklabel',x,'ylim',[-0.025 0.025],'ytick',-0.05:0.005:0.05)
xlabel('Time since SWR detection (s)')
box off

%Test significance using kruskalwallis anova
group = ones(size(lowrip)); 
for i = 1:size(rip,1)
    group(i,:) =i;
end
[p table stats] = kruskalwallis(reshape(lowrip,size(lowrip,1)*size(lowrip,2),1)-repmat(baseline_rip',size(rip,1),1),reshape(group,size(group,1)*size(group,2),1));
c = multcompare(stats);
%RUN:
% Time -200, -100ms < baseline at p<0.05
% 0, 100ms < baseline at p<0.001
% baseline = 0.5

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_ca1ca3_lowripple_coherence.pdf', m, d, y);
print('-dpdf', savestring)

