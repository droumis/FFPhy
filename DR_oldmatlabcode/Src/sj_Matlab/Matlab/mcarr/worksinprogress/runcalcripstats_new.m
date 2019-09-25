%% RIPPLE PAPER FIGURES

%Animals
animals = {'Bond','Frank','Ten'};

%Define epochs
epochfilter = [];
epochfilter{1} = 'isequal($type, ''run'')';
%epochfilter{1} = 'isequal($type,''sleep'')';

%Define iterator
iterator = 'epochbehaveanal';

%Define time filter
timefilter = {{'get2dstate', '($velocity < 4)'}};
%timefilter = {{'get2dstate', '($velocity < 4) & $immobilitytime > 60'}};

f = createfilter('animal',animals,'epochs',epochfilter,'excludetime',timefilter,'iterator', iterator);
f = setfilterfunction(f, 'calcripstats_new', {'rip'});
f = runfilter(f);
% save('/data21/mcarr/RipplePaper/ripstats.mat','f')
% save('/data21/mcarr/RipplePaper/ripstats_sleep.mat','f')
time = f(1).output{1}(1).time;

%% Correlation between gamma power and peak ripple power
ca1 = []; ca3 = [];
for an = 1:length(f)
    for d =1:length(f(an).output)
        for e = 1:length(f(an).output{d})
            if ~isempty(f(an).output{d}(e).ca1_power_corr)
                ca1 = [ca1 f(an).output{d}(e).ca1_power_corr(:,1)];
                ca3 = [ca3 f(an).output{d}(e).ca3_power_corr(:,1)];
            end
        end
    end
end

mean_1 = nanmean(ca1,2); se_1 = nanstd(ca1,[],2)./sqrt(sum(~isnan(ca1),2)-1);
mean_3 = nanmean(ca3,2); se_3 = nanstd(ca3,[],2)./sqrt(sum(~isnan(ca3),2)-1);

x = [-0.4:0.1:0.4 0.6];

figure
hold on
bar(x,mean_1,'b')
legend('CA1','location','NorthWest')
errorbar2(x,mean_1,se_1,'k')
set(gca,'xtick',x,'xticklabel',[-0.4:0.1:0.4 100],'ylim',[-0.1 0.5],'xlim',[-0.5 0.65],'ytick',-0.1:0.1:0.5)
ylabel('Correlation between gamma power and SWR amplitude')
xlabel('Time since ripple detection (s)')
box off

%Test significance using kruskalwallis anova
group = ones(size(ca1)); 
for i = 1:size(ca1,1)
    group(i,:) =i;
end
[p table stats] = kruskalwallis(reshape(ca1,size(ca1,1)*size(ca1,2),1),reshape(group,size(group,1)*size(group,2),1));
c = multcompare(stats);
%RUN:
% Time 300ms > baseline at p<0.01
% Time 100, 200ms, peak > baseline at p<1e-5

%SLEEP:
% Time 100-400ms, peak > baseline at p<0.001

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_2f.pdf', m, d, y);
print('-dpdf', savestring)

figure
hold on
bar(x,mean_3,'r')
legend('CA3','location','NorthWest')
errorbar2(x,mean_3,se_3,'k')
set(gca,'xtick',x,'xticklabel',[-0.4:0.1:0.4 100],'ylim',[-0.1 0.5],'xlim',[-0.5 0.65])
ylabel('Correlation between gamma power and SWR amplitude')
xlabel('Time since ripple detection (s)')
box off

%Test significance using kruskalwallis anova
group = ones(size(ca3)); 
for i = 1:size(ca3,1)
    group(i,:) =i;
end
[p table stats] = kruskalwallis(reshape(ca3,size(ca3,1)*size(ca3,2),1),reshape(group,size(group,1)*size(group,2),1));
c = multcompare(stats);
%RUN:
% Time 100, 200ms,peak > baseline at p<0.001

%SLEEP:
% Time 100-400ms, peak > baseline at p<0.05

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_2g.pdf', m, d, y);
print('-dpdf', savestring)

%% Plot gamma power relative to ripple onset
ca1_power = []; ca3_power = [];
for an = 1:length(f)
    for d =1:length(f(an).output)
        for e = 1:length(f(an).output{d})
            for i = 1:length(f(an).output{d}(e).peak)
                peak_ind = find(f(an).output{d}(e).peak(i)==f(an).output{d}(e).ca1_ripple_power(:,i));
                ca1_power = [ca1_power f(an).output{d}(e).ca1_power([6:10:end peak_ind],i)];
                ca3_power = [ca3_power f(an).output{d}(e).ca3_power([6:10:end peak_ind],i)];
            end
        end
    end
end
mean_1 = mean(ca1_power,2); mean_3 = mean(ca3_power,2);
se_1 = std(ca1_power,[],2)./sqrt(size(ca1_power,2)-1);
se_3 = std(ca3_power,[],2)./sqrt(size(ca3_power,2)-1);

x = [-0.4:0.1:0.4 0.6];

figure
hold on
bar(x,mean_1,'b')
legend('CA1','location','NorthWest')
errorbar2(x,mean_1,se_1,'k')
set(gca,'xtick',x,'xticklabel',[-0.4:0.1:0.4 100],'ylim',[-0.1 1.5],...
    'xlim',[-0.5 0.65],'ytick',-0:0.5:1.5)
ylabel('Normalized gamma power')
xlabel('Time since ripple detection (s)')
box off

%Test significance using kruskalwallis anova
group = ones(size(ca1_power)); 
for i = 1:size(ca1_power,1)
    group(i,:) =i;
end
[p table stats] = kruskalwallis(reshape(ca1_power,size(ca1_power,1)*size(ca1_power,2),1),reshape(group,size(group,1)*size(group,2),1));
c = multcompare(stats);
%RUN:
% Time 0-400ms,peak > baseline at p<1e-5
% Time -100ms > baseline at p<0.01
%SLEEP:
% Time -100-400ms,peak > baseline at p<1e-5

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_2c.pdf', m, d, y);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_7a.pdf', m, d, y);
print('-dpdf', savestring)

figure
hold on
bar(x,mean_3,'r')
legend('CA3','location','NorthWest')
errorbar2(x,mean_3,se_3,'k')
set(gca,'xtick',x,'xticklabel',[-0.4:0.1:0.4 100],'ylim',[-0.1 1.5],...
    'xlim',[-0.5 0.65],'ytick',-0:0.5:1.5)
ylabel('Normalized gamma power')
xlabel('Time since ripple detection (s)')
box off

%Test significance using kruskalwallis anova
group = ones(size(ca3_power)); 
for i = 1:size(ca3_power,1)
    group(i,:) =i;
end
[p table stats] = kruskalwallis(reshape(ca3_power,size(ca3_power,1)*size(ca3_power,2),1),reshape(group,size(group,1)*size(group,2),1));
c = multcompare(stats);
%RUN:
% Time 0-400ms,peak > baseline at p<1e-5

%SLEEP:
% Time 0-400ms,peak > baseline at p<1e-5

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_2d.pdf', m, d, y);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_7b.pdf', m, d, y);
print('-dpdf', savestring)

%% Plot gamma coherence relative to ripple onset

ca11 = []; ca13 = []; ca33 = []; baseline_11 = []; baseline_13 = []; baseline_33 = [];
for an = 1:length(f)
    for d =2%:length(f(an).output)
        for e = 1:length(f(an).output{d})
            if ~isempty(f(an).output{d}(e).ca1_ca1_coherence)
            for i = 1:length(f(an).output{d}(e).peak)
                peak_ind = find(f(an).output{d}(e).peak(i)==f(an).output{d}(e).ca1_ripple_power(:,i));
                ca11 = [ca11 f(an).output{d}(e).ca1_ca1_coherence([6:10:end peak_ind],i)];
                ca13 = [ca13 f(an).output{d}(e).ca1_ca3_coherence([6:10:end peak_ind],i)];
                ca33 = [ca33 f(an).output{d}(e).ca3_ca3_coherence([6:10:end peak_ind],i)];
            end
            baseline_11 = [baseline_11 mean(f(an).output{d}(e).ca1_ca1_coherence(1:6,:))];
            baseline_13 = [baseline_13 mean(f(an).output{d}(e).ca1_ca3_coherence(1:6,:))];
            baseline_33 = [baseline_33 mean(f(an).output{d}(e).ca3_ca3_coherence(1:6,:))];

            end
        end
    end
end
invalid = isnan(ca11(1,:)); ca11(:,invalid) =[];  baseline_11(invalid) = [];
invalid = isnan(ca13(1,:)); ca13(:,invalid) =[];  baseline_13(invalid) = [];
invalid = isnan(ca33(1,:)); ca33(:,invalid) =[];  baseline_33(invalid) = [];

mean_11 = mean(ca11-repmat(baseline_11,size(ca11,1),1),2);
mean_13 = mean(ca13-repmat(baseline_13,size(ca13,1),1),2);
mean_33 = mean(ca33-repmat(baseline_33,size(ca33,1),1),2);
se_11 = std(ca11-repmat(baseline_11,size(ca11,1),1),[],2)./sqrt(size(ca11,2)-1);
se_13 = std(ca13-repmat(baseline_13,size(ca13,1),1),[],2)./sqrt(size(ca13,2)-1);
se_33 = std(ca33-repmat(baseline_33,size(ca33,1),1),[],2)./sqrt(size(ca33,2)-1);

x = [-0.4:0.1:0.4 0.6];

figure
hold on
bar(x,mean_11,'b')
legend('CA1-CA1','location','NorthWest')
errorbar2(x,mean_11,se_11,'k')
set(gca,'xtick',x,'xticklabel',[-0.4:0.1:0.4 100],'ylim',[-0.01 0.075],...
    'xlim',[-0.5 0.65],'ytick',0:0.025:0.075)
ylabel('delta CA1-CA1 coherence')
xlabel('Time since ripple detection (s)')
box off

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_3h.pdf', m, d, y);
print('-dpdf', savestring)

figure
hold on
bar(x,mean_13,'m')
legend('CA1-CA3','location','NorthWest')
errorbar2(x,mean_13,se_13,'k')
set(gca,'xtick',x,'xticklabel',[-0.4:0.1:0.4 100],'ylim',[-0.01 0.075],...
    'xlim',[-0.5 0.65],'ytick',0:0.025:0.075)
ylabel('delta CA1-CA3 coherence')
xlabel('Time since ripple detection (s)')
box off

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_3e.pdf', m, d, y);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_7e.pdf', m, d, y);
print('-dpdf', savestring)

figure
hold on
bar(x,mean_33,'r')
legend('CA3-CA3','location','NorthWest')
errorbar2(x,mean_33,se_33,'k')
set(gca,'xtick',x,'xticklabel',[-0.4:0.1:0.4 100],'ylim',[-0.01 0.075],...
    'xlim',[-0.5 0.65],'ytick',0:0.025:0.075)
ylabel('delta CA3-CA3 coherence')
xlabel('Time since ripple detection (s)')
box off

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_3h_bottom.pdf', m, d, y);
print('-dpdf', savestring)

%Test significance using kruskalwallis anova
group = ones(size(ca11)); 
for i = 1:size(ca11,1)
    group(i,:) =i;
end
[p table stats] = kruskalwallis(reshape(ca11,size(ca11,1)*size(ca11,2),1)-repmat(baseline_11',size(ca11,1),1),reshape(group,size(group,1)*size(group,2),1));
c = multcompare(stats);
%RUN:
% Time 0-100ms, peak > baseline at p<1e-5
% Time 200ms > baseline at p<0.05
% baseline = 0.69

%SLEEP:
% only peak > baseline at p>0.05
%basline =  0.7094
group = ones(size(ca13)); 
for i = 1:size(ca13,1)
    group(i,:) =i;
end
[p table stats] = kruskalwallis(reshape(ca13,size(ca13,1)*size(ca13,2),1),reshape(group,size(group,1)*size(group,2),1));
c = multcompare(stats);
%RUN:
% Time 0-400ms,peak > baseline at p<1e-5
% baseline = 0.58

%SLEEP:
% Time 200-400ms > baseline at p<0.05
% Time 100ms, peak > baseline at p<1e-5
% baseline = 0.62

group = ones(size(ca33)); 
for i = 1:size(ca33,1)
    group(i,:) =i;
end
[p table stats] = kruskalwallis(reshape(ca33,size(ca33,1)*size(ca33,2),1),reshape(group,size(group,1)*size(group,2),1));
c = multcompare(stats);
%RUN:
% Time 0-400ms,peak > baseline at p<1e-5
% baseline = 0.78

%SLEEP:
% Time 0-300ms,peak > baseline at p<1e-5
% 400ms > baseline at p<0.001
% baseline = 0.7977


%% How does coherence relate to gamma power?
ca111 = []; ca131 = []; ca331 = []; ca113 = []; ca133 = []; ca333 = [];
for an = 1:length(f)
    for d =1%:length(f(an).output)
        for e = 1:length(f(an).output{d})
            tca11 = []; tca13 = []; tca33 = []; tca1 = []; tca3 = [];
            for i = 1:length(f(an).output{d}(e).peak)
                peak_ind = find(f(an).output{d}(e).peak(i)==f(an).output{d}(e).ca1_ripple_power(:,i));
                tca11 = [tca11 f(an).output{d}(e).ca1_ca1_coherence([6:10:end peak_ind],i)];
                tca13 = [tca13 f(an).output{d}(e).ca1_ca3_coherence([6:10:end peak_ind],i)];
                tca33 = [tca33 f(an).output{d}(e).ca3_ca3_coherence([6:10:end peak_ind],i)];
                tca1 = [tca1 f(an).output{d}(e).ca1_power([6:10:end peak_ind],i)];
                tca3 = [tca3 f(an).output{d}(e).ca3_power([6:10:end peak_ind],i)];
            end
            
            r = corrcoef(tca11,tca1);
            ca111 = [ca111; r(1,2)];
        
            r = corrcoef(tca13,tca1);
            ca131 = [ca131; r(1,2)];
                
            r = corrcoef(tca33,tca1);
            ca331 = [ca331; r(1,2)];

            r = corrcoef(tca11,tca3);
            ca113 = [ca113; r(1,2)];
      
            r = corrcoef(tca13,tca3);
            ca133 = [ca133; r(1,2)];
                
            r = corrcoef(tca33,tca3);
            ca333 = [ca333; r(1,2)];
        end
    end
end

ca111(isnan(ca111)) = []; ca131(isnan(ca131)) = []; ca331(isnan(ca331)) = [];
ca113(isnan(ca113)) = []; ca133(isnan(ca133)) = []; ca333(isnan(ca333)) = [];

mean_111 = mean(ca111); se_111 = std(ca111)./sqrt(length(ca111)-1);
mean_131 = mean(ca131); se_131 = std(ca131)./sqrt(length(ca131)-1);
mean_331 = mean(ca331); se_331 = std(ca331)./sqrt(length(ca331)-1);
mean_113 = mean(ca113); se_113 = std(ca113)./sqrt(length(ca113)-1);
mean_133 = mean(ca133); se_133 = std(ca133)./sqrt(length(ca133)-1);
mean_333 = mean(ca333); se_333 = std(ca333)./sqrt(length(ca33)-1);

figure
hold on
bar(1:3, [mean_111 mean_131 mean_331],'b')
bar(4:6, [mean_113 mean_133 mean_333],'r')
legend([{'CA1 gamma power'},{'CA3 gamma power'}])
errorbar2(1:3, [mean_111 mean_131 mean_331],[se_111 se_131 se_331],'k')
errorbar2(4:6, [mean_113 mean_133 mean_333],[se_113 se_133 se_333],'k')
set(gca,'xtick',1:1:6,'xticklabel',[{'CA1-CA1'},{'CA1-CA3'},{'CA3-CA3'},{'CA1-CA1'},{'CA1-CA3'},{'CA3-CA3'}],...
    'ylim',[0 0.6],'ytick',0:0.1:0.6)
ylabel('Correlation between gamma coherence and gamma power')


% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_coherence_power_correlation.pdf', m, d, y);
print('-dpdf', savestring)

