%% RIPPLE PAPER FIGURES

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
f = setfilterfunction(f, 'calcripstats_new', {'rip'});
f = runfilter(f);

time = f(1).output{1}(1).time;

%% Correlation between gamma power and peak ripple power
ca1 = []; ca3 = [];
for an = 1:length(f)
    for d =1:length(f(an).output)
        for e = 1:length(f(an).output{d})
            ca1 = [ca1 f(an).output{d}(e).ca1_power_corr(:,1)];
            ca3 = [ca3 f(an).output{d}(e).ca3_power_corr(:,1)];
        end
    end
end

mean_1 = mean(ca1,2); se_1 = std(ca1,[],2)./sqrt(size(ca1,2)-1);
mean_3 = mean(ca3,2); se_3 = std(ca3,[],2)./sqrt(size(ca3,2)-1);

x = [-0.4:0.1:0.4 0.6];
[h p1] = ttest(ca1'); [h p3] = ttest(ca3');
p1 = p1*size(ca1,1); p3 = p3*size(ca3,1);

figure
hold on
bar(x,mean_1,'b')
legend('CA1','location','NorthWest')
errorbar2(x,mean_1,se_1,'k')
plot(x(p1<0.05),0.5,'k*')
set(gca,'xtick',x,'xticklabel',[-0.4:0.1:0.4 100],'ylim',[-0.1 0.6],'xlim',[-0.5 0.65])
ylabel('Correlation between gamma power and SWR amplitude')
xlabel('Time since ripple detection (s)')
box off
%Note: p1 is significant p<1e-5

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_ca1_gamma_power_correlation.pdf', m, d, y);
print('-dpdf', savestring)

figure
hold on
bar(x,mean_3,'r')
legend('CA3','location','NorthWest')
errorbar2(x,mean_3,se_3,'k')
plot(x(p3<0.05),0.5,'k*')
set(gca,'xtick',x,'xticklabel',[-0.4:0.1:0.4 100],'ylim',[-0.1 0.6],'xlim',[-0.5 0.65])
ylabel('Correlation between gamma power and SWR amplitude')
xlabel('Time since ripple detection (s)')
box off
%Note, p3 is significant p<0.01

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_ca3_gamma_power_correlation.pdf', m, d, y);
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
[h p1] = ttest(ca1_power'); [h p3] = ttest(ca3_power');
p1 = p1*size(ca1_power,1); p3 = p3*size(ca3_power,1);

figure
hold on
bar(x,mean_1,'b')
legend('CA1','location','NorthWest')
errorbar2(x,mean_1,se_1,'k')
plot(x(p1<0.05),1,'k*')
set(gca,'xtick',x,'xticklabel',[-0.4:0.1:0.4 100],'ylim',[-0.1 1.2],...
    'xlim',[-0.5 0.65],'ytick',-0.1:0.1:1)
ylabel('Normalized gamma power')
xlabel('Time since ripple detection (s)')
box off
%Note: p1 is significant p<1e-5 except at -0.2, p<0.05

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_ca1_gamma_power_summary.pdf', m, d, y);
print('-dpdf', savestring)


figure
hold on
bar(x,mean_3,'r')
legend('CA3','location','NorthWest')
errorbar2(x,mean_3,se_3,'k')
plot(x(p3<0.05),1,'k*')
set(gca,'xtick',x,'xticklabel',[-0.4:0.1:0.4 100],'ylim',[-0.1 1.2],...
    'xlim',[-0.5 0.65],'ytick',-0.1:0.1:1)
ylabel('Normalized gamma power')
xlabel('Time since ripple detection (s)')
box off
%Note: p3 is significant p<1e-5

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_ca3_gamma_power_summary.pdf', m, d, y);
print('-dpdf', savestring)

%% Plot gamma coherence relative to ripple onset

ca11 = []; ca13 = []; ca33 = [];
for an = 1:length(f)
    for d =1:length(f(an).output)
        for e = 1:length(f(an).output{d})
            for i = 1:length(f(an).output{d}(e).peak)
                peak_ind = find(f(an).output{d}(e).peak(i)==f(an).output{d}(e).ca1_ripple_power(:,i));
                ca11 = [ca11 f(an).output{d}(e).ca1_ca1_coherence([6:10:end peak_ind],i)];
                ca13 = [ca13 f(an).output{d}(e).ca1_ca3_coherence([6:10:end peak_ind],i)];
                ca33 = [ca33 f(an).output{d}(e).ca3_ca3_coherence([6:10:end peak_ind],i)];
            end
        end
    end
end
invalid = isnan(ca11(1,:)); ca11(:,invalid) =[];
invalid = isnan(ca13(1,:)); ca13(:,invalid) =[];
invalid = isnan(ca33(1,:)); ca33(:,invalid) =[];

mean_11 = mean(ca11-repmat(ca11(1,:),size(ca11,1),1),2);
mean_13 = mean(ca13-repmat(ca13(1,:),size(ca13,1),1),2);
mean_33 = mean(ca33-repmat(ca33(1,:),size(ca33,1),1),2);
se_11 = std(ca11-repmat(ca11(1,:),size(ca11,1),1),[],2)./sqrt(size(ca11,2)-1);
se_11(1) = std(ca11(1,:))./sqrt(length(ca11(1,:))-1);
se_13 = std(ca13-repmat(ca13(1,:),size(ca13,1),1),[],2)./sqrt(size(ca13,2)-1);
se_13(1) = std(ca13(1,:))./sqrt(length(ca13(1,:))-1);
se_33 = std(ca33-repmat(ca33(1,:),size(ca33,1),1),[],2)./sqrt(size(ca33,2)-1);
se_33(1) = std(ca33(1,:))./sqrt(length(ca33(1,:))-1);

x = [-0.4:0.1:0.4 0.6];
[h p11] = ttest((ca11-repmat(ca11(1,:),size(ca11,1),1))');
[h p13] = ttest((ca13-repmat(ca13(1,:),size(ca13,1),1))');
[h p33] = ttest((ca33-repmat(ca33(1,:),size(ca33,1),1))'); clear h
p11 = p11*size(ca11,1); p13 = p13*size(ca13,1); p33 = p33*size(ca33,1);

figure
hold on
bar(x,mean_11,'b')
legend('CA1-CA1','location','NorthWest')
errorbar2(x,mean_11,se_11,'k')
plot(x(p11<0.05),0.05,'k*')
set(gca,'xtick',x,'xticklabel',[-0.4:0.1:0.4 100],'ylim',[-0.01 0.06],...
    'xlim',[-0.5 0.65],'ytick',-0.01:0.01:0.06)
ylabel('delta CA1-CA1 coherence')
xlabel('Time since ripple detection (s)')
box off
%Note: p11 is significant p<1e-5 except at -0.1, p<0.01

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_ca1ca1_coherence_summary.pdf', m, d, y);
print('-dpdf', savestring)

figure
hold on
bar(x,mean_13,'m')
legend('CA1-CA3','location','NorthWest')
errorbar2(x,mean_13,se_13,'k')
plot(x(p13<0.05),0.05,'k*')
set(gca,'xtick',x,'xticklabel',[-0.4:0.1:0.4 100],'ylim',[-0.01 0.06],...
    'xlim',[-0.5 0.65],'ytick',-0.01:0.01:0.06)
ylabel('delta CA1-CA3 coherence')
xlabel('Time since ripple detection (s)')
box off
%Note: p13 is significant p<1e-5 except at -0.1, p<0.01

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_ca1ca3_coherence_summary.pdf', m, d, y);
print('-dpdf', savestring)

figure
hold on
bar(x,mean_33,'r')
legend('CA3-CA3','location','NorthWest')
errorbar2(x,mean_33,se_33,'k')
plot(x(p33<0.05),0.05,'k*')
set(gca,'xtick',x,'xticklabel',[-0.4:0.1:0.4 100],'ylim',[-0.01 0.06],...
    'xlim',[-0.5 0.65],'ytick',-0.01:0.01:0.06)
ylabel('delta CA3-CA3 coherence')
xlabel('Time since ripple detection (s)')
box off
%Note: p33 is significant p<1e-5

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_ca3ca3_coherence_summary.pdf', m, d, y);
print('-dpdf', savestring)


%% Correlation between gamma coherence and peak ripple power

ca11 = []; ca13 = []; ca33 = [];

for an = 1:length(f)
    for d =1:length(f(an).output)
        for e = 1:length(f(an).output{d})
            ca11 = [ca11 f(an).output{d}(e).ca1_ca1_coherence_corr(:,1)];
            ca13 = [ca13 f(an).output{d}(e).ca1_ca3_coherence_corr(:,1)];
            ca33 = [ca33 f(an).output{d}(e).ca3_ca3_coherence_corr(:,1)];
        end
    end
end
invalid = isnan(ca11(1,:)); ca11(:,invalid) =[];
invalid = isnan(ca13(1,:)); ca13(:,invalid) =[];
invalid = isnan(ca33(1,:)); ca33(:,invalid) =[];

mean_11 = mean(ca11,2); se_11 = std(ca11,[],2)./sqrt(size(ca11,2)-1);
mean_13 = mean(ca13,2); se_13 = std(ca13,[],2)./sqrt(size(ca13,2)-1);
mean_33 = mean(ca33,2); se_33 = std(ca33,[],2)./sqrt(size(ca33,2)-1);

x = [-0.4:0.1:0.4 0.6];
[h p11] = ttest(ca11'); [h p13] = ttest(ca13'); [h p33] = ttest(ca33');
p11 = p11*size(ca11,1); p13 = p13*size(ca13,1); p33 = p33*size(ca33,1);

figure
hold on
bar(x,mean_11,'b')
legend('CA1-CA1','location','NorthWest')
errorbar2(x,mean_11,se_11,'k')
plot(x(p11<0.05),0.2,'k*')
set(gca,'xtick',x,'xticklabel',[-0.4:0.1:0.4 100],'ylim',[-0.1 0.3],'xlim',[-0.5 0.65])
ylabel('Correlation between CA1-CA1 gamma coherence and SWR amplitude')
xlabel('Time since ripple detection (s)')
box off
%Note: p11 is significant p<0.01 at 0.1 0.2 0.3 0.4 peak and p<0.05 at 0

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_ca1ca1_gamma_coherence_correlation.pdf', m, d, y);
print('-dpdf', savestring)

figure
hold on
bar(x,mean_13,'m')
legend('CA1-CA3','location','NorthWest')
errorbar2(x,mean_13,se_13,'k')
plot(x(p13<0.05),0.2,'k*')
set(gca,'xtick',x,'xticklabel',[-0.4:0.1:0.4 100],'ylim',[-0.1 0.3],'xlim',[-0.5 0.65])
ylabel('Correlation between CA1-CA3 gamma coherence and SWR amplitude')
xlabel('Time since ripple detection (s)')
box off
%Note: p13 is significant p<0.01

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_ca1ca3_gamma_coherence_correlation.pdf', m, d, y);
print('-dpdf', savestring)

figure
hold on
bar(x,mean_33,'r')
legend('CA3-CA3','location','NorthWest')
errorbar2(x,mean_33,se_33,'k')
plot(x(p33<0.05),0.2,'k*')
set(gca,'xtick',x,'xticklabel',[-0.4:0.1:0.4 100],'ylim',[-0.1 0.3],'xlim',[-0.5 0.65])
ylabel('Correlation between CA3-CA3 gamma coherence and SWR amplitude')
xlabel('Time since ripple detection (s)')
box off
%Note: p33 is significant p<0.01


% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_ca3ca3_gamma_coherence_correlation.pdf', m, d, y);
print('-dpdf', savestring)


%Test significance using Anova instead
group = ones(size(ca11)); 
for i = 1:size(ca11,1)
    group(i,:) =i;
end
[p table stats] = kruskalwallis(reshape(ca11,size(ca11,1)*size(ca11,2),1),reshape(group,size(group,1)*size(group,2),1));
c = multcompare(stats);

group = ones(size(ca13)); 
for i = 1:size(ca13,1)
    group(i,:) =i;
end
[p table stats] = kruskalwallis(reshape(ca13,size(ca13,1)*size(ca13,2),1),reshape(group,size(group,1)*size(group,2),1));
c = multcompare(stats);

group = ones(size(ca33)); 
for i = 1:size(ca33,1)
    group(i,:) =i;
end
[p table stats] = kruskalwallis(reshape(ca33,size(ca33,1)*size(ca33,2),1),reshape(group,size(group,1)*size(group,2),1));
c = multcompare(stats);

%% How does coherence relate to gamma power?
ca111 = []; ca131 = []; ca331 = []; ca113 = []; ca133 = []; ca333 = [];
for an = 1:length(f)
    for d =1:length(f(an).output)
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

