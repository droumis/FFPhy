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
f = setfilterfunction(f, 'calcrippleduration', {'ripples','cellinfo'});
f = runfilter(f);

%% HISTOGRAM OF SWR DURATION

g = numericgroupcombine(f);
g = g{1}; g(isnan(g)) = [];

figure
hold on
subs = 0:0.01:0.5;
val = lookup(g(1,:),subs,1);
boxplot(g(1,:),'colors','k','orientation','horizontal','notch','on','whisker',0,'symbol','k')
plot(subs,hist(val,1:length(subs))./length(val))
set(gca,'xlim',[0 0.5],'ylim',[0 1.1],'ytick',0:0.05:0.2,'ytickLabel',0:0.05:0.2)
xlabel('SWR duration (s)')
ylabel('Proportion')
box off

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_distribution_SWR_duration.pdf', m, d, y);
print('-dpdf', savestring)

