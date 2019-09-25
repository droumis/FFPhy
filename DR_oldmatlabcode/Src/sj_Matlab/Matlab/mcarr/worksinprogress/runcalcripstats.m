%% RIPPLE PAPER FIGURES: POWER

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
f = setfilterfunction(f, 'calcripsstats', {'rips'});
f = runfilter(f);

time = f(1).output{1}(1).time;
frequency = f(1).output{1}(1).frequency;
%% CA1 and CA3 spectrum
ca1 = zeros(length(time),length(frequency)); ca3 = zeros(length(time),length(frequency)); count1 = 0; count3 = 0;
for an = 1:length(f)
    for d = 1:length(f(an).output)
        for e = 1:length(f(an).output{d})
            if ~isempty(f(an).output{d}(e).ca1_spectrum)
                ca1 = ca1 + f(an).output{d}(e).ca1_spectrum;
                count1 = count1+f(an).output{d}(e).nrips;
            end
            if ~isempty(f(an).output{d}(e).ca3_spectrum)
                ca3 = ca3 + f(an).output{d}(e).ca3_spectrum;
                count3 = count3+f(an).output{d}(e).nrips;
            end
        end
    end
end
ca1 = ca1./count1; ca3 = ca3./count3;

%Plot CA1 rip triggered spectrum
figure
surf(time,frequency,ca1')
view(0,90); shading interp; colormap hot
set(gca,'clim',[0 2.5],'xlim',[time(26) time(86)],'ylim',f(1).output{1}(1).frequency([1 end]),'xtick',-0.2:0.1:0.4)
colorbar('ytick',0:0.25:1.5)

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_ca1_riptriggered_spectrum_allruns.png', m, d, y);
print('-dpng', savestring)

%Plot CA3 rip triggered spectrum
figure
surf(time,frequency,ca3')
view(0,90); shading interp;  colormap hot
set(gca,'clim',[0 1.5],'xlim',[time(26) time(86)],'ylim',f(1).output{1}(1).frequency([1 end]),'xtick',-0.2:0.1:0.4)
colorbar('ytick',0:0.25:1.5)

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_ca3_riptriggered_spectrum_allruns.png', m, d, y);
print('-dpng', savestring)

%% Plot examples
load '/data13/mcarr/Ten/tenrips04.mat'
rips = rips{4}{6};

time = rips.time;
frequency = rips.frequency;

load '/data13/mcarr/Ten/tenpos04.mat'
speed = pos{4}{6}.data(lookup(rips.starttime,pos{4}{6}.data(:,1)),8);
minVelocity = 4;

%CA1 spectrum
ca1 = mean(rips.ca1_spectrum(:,:,speed<minVelocity),3);
figure
surf(time,frequency,ca1')
view(0,90); shading interp; colormap hot
set(gca,'clim',[-0.5 3],'xlim',[time(6) time(86)],'xtick',time(6:10:86),'xticklabel',-0.4:0.1:0.4,'ylim',frequency([1 end]))
colorbar('ytick',-0.5:0.5:3)

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_2a.png', m, d, y);
print('-dpng', savestring)


%CA3 spectrum
ca3 = mean(rips.ca3_spectrum(:,:,speed<minVelocity),3);
figure
surf(time,frequency,ca3')
view(0,90); shading interp; colormap hot
set(gca,'clim',[-0.5 2],'xlim',[time(6) time(86)],'xtick',time(6:10:86),'xticklabel',-0.4:0.1:0.4,'ylim',frequency([1 end]))
colorbar('ytick',-0.5:0.5:3)

%Save Figure
[y,m,d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_2b.png', m, d, y);
print('-dpng', savestring)


%% Plot Single SWR examples
load '/data13/mcarr/Ten/tenrips04.mat'
rips = rips{4}{6};

time = rips.time;
frequency = rips.frequency;

load '/data13/mcarr/Ten/tenpos04.mat'
speed = pos{4}{6}.data(lookup(rips.starttime,pos{4}{6}.data(:,1)),8);
minVelocity = 4;

%CA1 1
figure
event = [6 11 16];
for i = 1:length(event)
    subplot(1,3,i)
    surf(time,frequency,rips.ca1_spectrum(:,:,event(i))');
    view(0,90); shading interp; colormap hot
    set(gca,'clim',[-1 8],'xlim',[time(6) time(86)],'xtick',time(6:10:86),'xticklabel',-0.4:0.1:0.4,'ylim',frequency([1 end]))
end


% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_ca1_spectrum_sampleSWRs.png', m, d, y);
print('-dpng', savestring)



%% RIPPLE PAPER FIGURES: COHERENCE

%Animals
animals = {'Bond','Frank','Ten'};

%Define epochs
epochfilter = [];
%epochfilter{1} = 'isequal($type, ''run'')';
epochfilter{1} = 'isequal($type,''sleep'')';

%Define iterator
iterator = 'epochbehaveanal';

%Define time filter
%timefilter = {{'get2dstate', '($velocity < 4)'}};
timefilter = {{'get2dstate', '($velocity < 4) & $immobilitytime > 60'}};

f = createfilter('animal',animals,'epochs',epochfilter,'excludetime',timefilter,'iterator', iterator);
f = setfilterfunction(f, 'calcripcstats', {'ripc'});
f = runfilter(f);

time = f(1).output{1}(1).time;
frequency = f(1).output{1}(1).frequency;
%save('/data13/mcarr/RipplePaper/ripple_phaselocking.mat','f')

%% CA1 and CA3 phase coherence
p13 = [];
for an =1:length(f)
    for d = 1:length(f(an).output)
        for e = 1:length(f(an).output{d})
            if ~isempty(f(an).output{d}(e).ca1_ca1_phase_std)
                p13 = [p13 f(an).output{d}(e).ca1_ca3_phase_std(:,2)];
            end
        end
    end
end

baseline13 = repmat(nanmean(p13(1:5,:)),size(p13,1),1); 

x = -0.4:0.1:0.4;
baseline13 = baseline13(6:10:end,:); p13 = p13(6:10:end,:);
mean_13 = nanmean(p13-baseline13,2); se_13 = nanstd(p13-baseline13,[],2)./sqrt(sum(~isnan(p13),2)-1);


figure
hold on
bar(1:length(mean_13),mean_13,'b')
errorbar2(1:length(mean_13),mean_13,se_13,'k')
set(gca,'xtick',1:9,'xticklabel',x,'ylim',[-0.01 0.075],'ytick',0:0.025:0.075)
ylabel('delta phase locking')
xlabel('Time since ripple detection (s)')
box off

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_3b.pdf', m, d, y);
print('-dpdf', savestring)

group = ones(size(p13));
for i = 1:size(p13,1)
    group(i,:) = i;
end
[p table stats]=kruskalwallis(reshape(p13-baseline13,size(p13,1)*size(p13,2),1),reshape(group,size(group,1)*size(group,2),1));
multcompare(stats);

%RUN:
% 100-200ms >baseline at p<1e-5
% 0,300-400ms > baseline at p<0.001
% -100ms > baseline at p < 0.05
%mean baseline for phase locking: 0.82

figure
hold on
bar(1:length(mean_13),mean_13,'b')
errorbar2(1:length(mean_13),mean_13,se_13,'k')
set(gca,'xtick',1:9,'xticklabel',x,'ylim',[-0.02 0.03],'ytick',-0.02:0.01:0.04)
ylabel('delta phase locking')
xlabel('Time since ripple detection (s)')
box off

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_7d.pdf', m, d, y);
print('-dpdf', savestring)

%SLEEP:
% None significantly different
%mean baseline: 0.87

%% CA1 -CA1 and CA3-CA3 phase coherence

p11 = []; p33 = [];
for an = 1:length(f)
    for d = 1:length(f(an).output)
        for e = 1:length(f(an).output{d})
            if ~isempty(f(an).output{d}(e).ca1_ca1_phase_std)
                p11 = [p11 f(an).output{d}(e).ca1_ca1_phase_std(:,2)];
            end
            if ~isempty(f(an).output{d}(e).ca3_ca3_phase_std)
                p33 = [p33 f(an).output{d}(e).ca3_ca3_phase_std(:,2)];
            end
        end
    end
end

baseline11 = repmat(nanmean(p11(1:6,:)),size(p11,1),1); baseline11 = baseline11(6:10:end,:);
baseline33 = repmat(nanmean(p33(1:6,:)),size(p33,1),1); baseline33 = baseline33(6:10:end,:);
p11 = p11(6:10:end,:);
p33 = p33(6:10:end,:);

mean_11= nanmean(p11-baseline11,2); se_11 = nanstd(p11-baseline11,[],2)./sqrt(sum(~isnan(p11),2)-1);
mean_33= nanmean(p33-baseline33,2); se_33 = nanstd(p33-baseline33,[],2)./sqrt(sum(~isnan(p33),2)-1);

x = -0.4:0.1:0.4;

figure
hold on
bar(1:length(mean_11),mean_11,'b')
legend({'CA1-CA1'})
errorbar2(1:length(mean_11),mean_11,se_11,'k')
set(gca,'xtick',1:9,'xticklabel',x,'ylim',[-0.01 0.04],'ytick',0:0.02:0.04)
ylabel('delta phase locking')
xlabel('Time since ripple detection (s)')
box off

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_3g.pdf', m, d, y);
print('-dpdf', savestring)

group = ones(size(p11));
for i = 1:size(p11,1)
    group(i,:) = i;
end
[p table stats]=kruskalwallis(reshape(p11-baseline11,size(p11,1)*size(p11,2),1),reshape(group,size(group,1)*size(group,2),1));
multcompare(stats,0.01);
%RUN:
% 100ms >baseline at p<1e-5
% 0ms > baseline at p<0.001
% 200-300ms > baseline at p < 0.05
%mean baseline for phase locking: 0.91

figure
hold on
bar(1:length(mean_33),mean_33,'b')
legend({'CA3-CA3'})
errorbar2(1:length(mean_33),mean_33,se_33,'k')
set(gca,'xtick',1:9,'xticklabel',x,'ylim',[-0.01 0.04],'ytick',0:0.02:0.04)
ylabel('delta phase locking')
xlabel('Time since ripple detection (s)')
box off

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_3g_bottom.pdf', m, d, y);
print('-dpdf', savestring)

group = ones(size(p33));
for i = 1:size(p33,1)
    group(i,:) = i;
end
[p table stats]=kruskalwallis(reshape(p33-baseline33,size(p33,1)*size(p33,2),1),reshape(group,size(group,1)*size(group,2),1));
multcompare(stats,0.01);
%RUN:
% 100ms >baseline at p<1e-5
% 0ms,200-400ms > baseline at p<0.05
%mean baseline for phase locking: 0.98

%% Histogram of frequency with peak coherence

p11 = []; p13 = []; p33 = [];
for an =1:length(f)
    for d = 1:length(f(an).output)
        for e = 1:length(f(an).output{d})
            if ~isempty(f(an).output{d}(e).peak11)
                p11 = [p11; f(an).output{d}(e).peak11];
                p13 = [p13; f(an).output{d}(e).peak13];
                p33 = [p33; f(an).output{d}(e).peak33];
            end
            
        end
    end
end

plot(frequency,hist(p33,frequency)./length(p33),'r')
hold on
plot(frequency,hist(p13,frequency)./length(p13),'c')
plot(frequency,hist(p11,frequency)./length(p11),'b')


%% Plot examples
load '/data13/mcarr/Ten/tenripc04.mat'
ripc = ripc{4}{6};

time = ripc.time;
frequency = ripc.frequency;

load '/data13/mcarr/Ten/tenpos04.mat'
speed = pos{4}{6}.data(lookup(ripc.starttime,pos{4}{6}.data(:,1)),8);
minVelocity = 4;

%CA1-CA3 coherence magnitude
ca13 = mean(ripc.ca1_ca3_coherence(:,:,speed<minVelocity),3);
figure
surf(time,frequency,ca13')
view(0,90); shading interp; colormap hot
set(gca,'clim',[0.4 0.7],'xlim',[time(6) time(86)],'xtick',time(6:10:86),'xticklabel',-0.4:0.1:0.4,'ylim',frequency([1 end]))
colorbar('ytick',0.4:0.05:0.7)

%Save Figure
[y,m,d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_3d.png', m, d, y);
print('-dpng', savestring)

%Phase coherence distributions:
bin = -pi+pi/40:pi/20:pi-pi/40;
phase_distribution = zeros(length(time),length(bin));
phase_locking = zeros(length(time),1);
for i = 1:length(time)
    phase_distribution(i,:) = hist(ripc.ca1_ca3_phase(i,speed<minVelocity),bin);
    [m r] = anglemean(ripc.ca1_ca3_phase(i,speed<minVelocity));
    phase_locking(i) = r;
end
baseline = repmat(mean(phase_locking(1:6)),size(time,2),1); baseline = baseline(6:10:end);
phase_locking = phase_locking(6:10:end);

figure
imagesc(time,bin,phase_distribution'./sum(speed<minVelocity))
set(gca,'clim',[0 0.2],'xlim',[-0.4 0.4],'ytick',-pi:pi/2:pi,'yticklabel',[{'-pi'},{'-pi/2'},{'0'},{'pi/2'},{'pi'}])
axis xy
colorbar; colormap hot
ylabel('Gamma phase')
xlabel('Time since ripple detection (s)')
box off

%Save figure
[y,m,d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_3a.png', m, d, y);
print('-dpng', savestring)

figure
bar(time(6:10:end),phase_locking-baseline,'k')
set(gca,'xlim',[-0.5 0.5],'xtick',-.4:0.1:0.4,'ytick',0:0.025:0.075,'ylim',[-0.01 0.075])
ylabel('delta phase locking')
xlabel('Time since ripple detection (s)')
box off

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_3a_top.pdf', m, d, y);
print('-dpdf', savestring)

