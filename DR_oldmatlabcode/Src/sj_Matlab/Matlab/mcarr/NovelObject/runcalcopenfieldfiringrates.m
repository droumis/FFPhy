%% RUN FILTER TO CALCULATE OCCUPANCY NORMALIZED PLACE FIELDS

% Animal Selection
animals = {'CorrianderNO','Cyclops','Dunphy','Fafnir','Godzilla','Grendel','Cml','Nico'};

% Epoch selection
epochfilter{1} = 'isequal($session,''familiar'')';
epochfilter{2} = 'isequal($session,''novel'')';
epochfilter{3} = 'isequal($session,''supernovel'')';

%Define time filter
timefilter = {{'get2dstate', '($velocity > 1)'},{'getriptimes', '($nripples == 0)', [], 'cellfilter', 'isequal($area, ''CA1'')'}};

%RUN FOR CA1
% Cell Filter
cellfilter = 'isequal($area, ''CA1'') & $meanrate < 7 & $numspikes > 100';
iterator = 'multicellanal';

f = createfilter('animal',animals,'epochs',epochfilter,'cells',cellfilter,'excludetimefilter', timefilter, 'iterator', iterator);
f = setfilterfunction(f, 'calcopenfieldfiringrates', {'spikes','pos','task'});
f = runfilter(f);

%RUN FOR CA3
% Cell Filter
cellfilter = 'isequal($area, ''CA3'') & $meanrate < 7 & $numspikes > 100';
iterator = 'multicellanal';

g = createfilter('animal',animals,'epochs',epochfilter,'cells',cellfilter,'excludetimefilter', timefilter, 'iterator', iterator);
g = setfilterfunction(g, 'calcopenfieldfiringrates', {'spikes','pos','task'});
g = runfilter(g);

%% INITIALIZE VARIABLES
ca1 = calcfiringratestats(f); 
ca3 = calcfiringratestats(g);

save('/data6/monster/analysis/ca1_openfieldfiringrates.mat','ca1')
save('/data6/monster/analysis/ca3_openfieldfiringrates.mat','ca3')

%% COMPARE REMAPPING BETWEEN CA1 AND CA3

figure
hold on
[y1 x1] = ecdf([ca1.corrn; ca1.corrs]);
[y3 x3] = ecdf([ca3.corrn; ca3.corrs]);
stairs(x1,y1,'r')
stairs(x3,y3,'k')
xlabel('Similarity between sessions')
ylabel('Cummulative proportion')
legend('CA1','CA3','location','NorthWest')

figure
subplot(1,2,1)
hold on
[y1 x1] = ecdf([ca1.corrn]);
[y3 x3] = ecdf([ca3.corrn]);
stairs(x1,y1,'r')
stairs(x3,y3,'k')
xlabel('Similarity between sessions')
ylabel('Cummulative proportion')
legend('CA1','CA3','location','NorthWest')
set(gca,'xlim',[-0.1 1])

subplot(1,2,2)
hold on
[y1 x1] = ecdf([ca1.corrs]);
[y3 x3] = ecdf([ca3.corrs]);
stairs(x1,y1,'r')
stairs(x3,y3,'k')
xlabel('Similarity between sessions')
ylabel('Cummulative proportion')
legend('CA1','CA3','location','NorthWest')
set(gca,'xlim',[-0.1 1])

%The median similarity is higher for CA3 place cells than for CA1
p = ranksum([ca1.corrn ca1.corrs],[ca3.corrn ca3.corrs]); %pvalue < 0.01
p = ranksum(ca1.corrn,ca3.corrn); %pvalue < 0.05
p = ranksum(ca1.corrs,ca3.corrs); %pvalue < 0.05

figure
hold on
bar([1:2 4:5],[median(ca1.corrn) median(ca3.corrn) median(ca1.corrs) median(ca3.corrs)],'b')
plot([1:2 4:5],[prctile(ca1.corrn,[25 50 75]); prctile(ca3.corrn,[25 50 75]); prctile(ca1.corrs,[25 50 75]); prctile(ca3.corrs,[25 50 75])],'ko')
set(gca,'xtick',[1 2 4 5],'xticklabel',[{'CA1'},{'CA3'}])
ylabel('Similarity between sessions')

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/Users/maggiecarr/Desktop/%d_%d_%d_placefieldcorrelation_bargraph.pdf', m, d, y);
print('-dpdf', savestring)

%% SUBSAMPLE CA1 NEURONS TO #OF CA3 NEURONS RECORDED
nboot = 5000;
Zn = size(ca3.corrn,1); Zs = size(ca3.corrs,1);
Xn = zeros(nboot,1); Xs = zeros(nboot,1);
for s = 1:nboot
    boot = ceil(size(ca1.corrn,1)*rand(Zn,1));
    Xn(s) = median(ca1.corrn(boot));
    boot = ceil(size(ca1.corrs,1)*rand(Zs,1));
    Xs(s) = median(ca1.corrs(boot));
end
%1-sum(median(ca3.corrn)>Xn)./nboot; p < 0.05
%1-sum(median(ca3.corrs)>Xs)./nboot; p < 0.01

%% COMPARE MEAN RATES IN NOVEL VS. FAMILIAR SESSIONS FOR CA3 AND CA1

n1 = ca1.nmean; n3 = ca3.nmean; s1 = ca1.smean; s3 = ca3.smean;
invalid = isnan(n1); n1(invalid) = []; invalid = isnan(n3); n3(invalid) = [];
invalid = isnan(s1); s1(invalid) = []; invalid = isnan(s3); s3(invalid) = []; clear invalid
meandiff = [n1; s1]; meandiff3 = [n3; s3];

figure
hold on
bar(1:2:5,[mean(n1) mean(s1) mean(meandiff)],'b')
bar(2:2:6,[mean(n3) mean(s3) mean(meandiff3)],'r')
legend('CA1','CA3')
errorbar(1:6,[mean(n1) mean(n3) mean(s1) mean(s3) mean(meandiff) mean(meandiff3)],[stderr(n1) stderr(n3) stderr(s1) stderr(s3) stderr(meandiff) stderr(meandiff3)],'k')
set(gca,'xtick',1.5:2:5.5,'xticklabel',[{'Novel vs. Familiar'},{'Super vs. Familiar'},{'All Novelty'}],'ylim',[0 1.75])
ylabel('Novel/Familiar mean rate ratio')
box off

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/Users/maggiecarr/Desktop/%d_%d_%d_meanrate.pdf', m, d, y);
print('-dpdf', savestring)


%The mean firing rate is significantly higher in novel sessions for CA1, no difference for CA3
p = signrank(n1,1); %pvalue < 0.001
p = signrank(n3,1); %pvalue > 0.3
p = signrank(s1,1); %pvalue < 0.05
p = signrank(s3,1); %pvalue > 0.3
p = signrank(meandiff,1); %pvalue<1e-5
p = signrank(meandiff3,1); %pvalue>0.15

%% COMPARE MEAN RATE BY QUADRANT

n1 = [ca1.nmean ca1.typen]; n3 = [ca3.nmean ca3.typen];
s1 = [ca1.smean ca1.types]; s3 = [ca3.smean ca3.types];
invalid = isnan(n1(:,1)); n1(invalid,:) = []; invalid = isnan(n3(:,1)); n3(invalid,:) = [];
invalid = isnan(s1(:,1)); s1(invalid,:) = []; invalid = isnan(s3(:,1)); s3(invalid,:) = []; clear invalid
meandiff = [n1; s1]; meandiff3 = [n3; s3];

figure
hold on
bar([1 2],[median(meandiff(meandiff(:,2)==1|meandiff(:,2)==3,1)) median(meandiff(meandiff(:,2)==2|meandiff(:,2)==4,1))],'b')
bar([3 4],[median(meandiff3(meandiff3(:,2)==1|meandiff3(:,2)==3,1)) median(meandiff3(meandiff3(:,2)==2|meandiff3(:,2)==4,1))],'r')
legend('CA1','CA3')
plot([1 3 2 4],[prctile(meandiff(meandiff(:,2)==1|meandiff(:,2)==3,1),[25 50 75]); prctile(meandiff3(meandiff3(:,2)==1|meandiff3(:,2)==3,1),[25 50 75]); ...
    prctile(meandiff(meandiff(:,2)==2|meandiff(:,2)==4,1),[25 50 75]); prctile(meandiff3(meandiff3(:,2)==2|meandiff3(:,2)==4,1),[25 50 75])],'ko')
set(gca,'xtick',1:4,'xticklabel',[{'NovelQ'},{'FamiliarQ'}],'ylim',[0 2])
ylabel('Novel/Familiar mean rate ratio')
box off

%No significant difference in the rate remapping between novel object and
%familiar object quadrant for eather CA1 or CA3
ranksum(meandiff(meandiff(:,2)==1|meandiff(:,2)==3),meandiff(meandiff(:,2)==2|meandiff(:,2)==4)); %p>0.4
ranksum(meandiff3(meandiff3(:,2)==1|meandiff3(:,2)==3),meandiff3(meandiff3(:,2)==2|meandiff3(:,2)==4)); %p>0.15


% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/Users/maggiecarr/Desktop/%d_%d_%d_meanrate_byquad.pdf', m, d, y);
print('-dpdf', savestring)

%% ARE UNSTABLE CELLS MORE LIKELY TO FIRE MORE?

c1 = [ca1.corrn; ca1.corrs]; c3 = [ca3.corrn; ca3.corrs];
m1 = [ca1.nmean; ca1.smean]; m3 = [ca3.nmean; ca3.smean];
invalid = isnan(c1) | isnan(m1); c1(invalid) = []; m1(invalid) = [];
invalid = isnan(c3) | isnan(m3); c3(invalid) = []; m3(invalid) = [];

%statistics
[r p] = corr(c1,m1); %r = -0.25, p<0.01
[r p] = corr(c3,m3); %r = 0.08, p>0.5
[r p] = corr(c1,m1,'type','Spearman'); %r = -0.28, p<0.001
[r p] = corr(c3,m3,,'type','Spearman'); %r = 0.09, p>0.5

%Plot
bin = [0.2 0.4 0.6 1];
subs = lookup(c1,bin);
a1 = accumarray(subs,m1,[length(bin) 1],@(x) mean(x));
a1_e = accumarray(subs,m1,[length(bin) 1],@(x) stderr(x));
subs = lookup(c3,bin);
a3 = accumarray(subs,m3,[length(bin) 1],@(x) mean(x));
a3_e = accumarray(subs,m3,[length(bin) 1],@(x) stderr(x));

figure
hold on
bar(bin,a1,'b')
bar(bin+0.1, a3,'r')
errorbar(bin,a1,a1_e,'k')
errorbar(bin+0.1,a3,a3_e,'k')


%% How is stability and firing rate related to distance from novel object

c1 = [ca1.corrn; ca1.corrs]; c3 = [ca3.corrn; ca3.corrs];
m1 = [ca1.nmean; ca1.smean]; m3 = [ca3.nmean; ca3.smean];
d1 = [ca1.nlocation; ca1.slocation]; d3 = [ca3.nlocation; ca3.slocation];

invalid = isnan(c1) | isnan(m1) | isnan(d1); c1(invalid) = []; m1(invalid) = []; d1(invalid) = [];
invalid = isnan(c3) | isnan(m3) | isnan(d3); c3(invalid) = []; m3(invalid) = []; d3(invalid) = [];

%statistics
[r p] = corr(c1,m1); %r = -0.25, p<0.01
[r p] = corr(c3,m3); %r = 0.08, p>0.5
[r p] = corr(c1,m1,'type','Spearman'); %r = -0.28, p<0.001
[r p] = corr(c3,m3,,'type','Spearman'); %r = 0.09, p>0.5

%Plot
bin = [min([d1;d3]) max([d1;d3])];
subs = lookup(d1,bin);
a1 = accumarray(subs,c1,[length(bin) 1],@(x) mean(x));
a1_e = accumarray(subs,c1,[length(bin) 1],@(x) stderr(x));
subs = lookup(d3,bin);
a3 = accumarray(subs,c3,[length(bin) 1],@(x) mean(x));
a3_e = accumarray(subs,c3,[length(bin) 1],@(x) stderr(x));

figure
hold on
bar(bin,a1,'b')
bar(bin+0.1, a3,'r')
errorbar(bin,a1,a1_e,'k')
errorbar(bin+0.1,a3,a3_e,'k')

%% PLOT EXAMPLES

load('/data6/monster/Cyc/Cycpos03.mat');
load('/data6/monster/Cyc/Cycspikes03.mat');
load('/data6/monster/Cyc/Cycripples03.mat');
load('/data6/monster/Cyc/Cyccellinfo.mat');

day =3; epochs = [2 4 6]; index = epochs'; index(:,2) = day;
p = pos{day}; s = spikes{day};
riptimes = getriptimes('/data6/monster/Cyc/','Cyc',[index(:,2) index(:,1)],[],'cellfilter','isequal($area,''CA1'')');
    
cells = [5 1; 11 2];
close all
for i = 1:length(epochs)
    figure(i)
    hold on
    
    %Plot valid positions
    valid = p{epochs(i)}.data(:,5) >= 1;
    plot(p{epochs(i)}.data(valid,2),p{epochs(i)}.data(valid,3),'k.')
    
    %Plot valid spikes for CA3 cell
    spikerip = lookup(s{epochs(i)}{cells(1,1)}{cells(1,2)}.data(:,1),riptimes{day}{epochs(i)}.time);
    spikerip = riptimes{day}{epochs(i)}.nripples(spikerip);
    
    spikepos = lookup(s{epochs(i)}{cells(1,1)}{cells(1,2)}.data(:,1),p{epochs(i)}.data(:,1));
    spikepos = p{epochs(i)}.data(spikepos,5);

    valid = logical(spikerip==0) & logical(spikepos > 1)';
    plot(s{epochs(i)}{cells(1,1)}{cells(1,2)}.data(valid,2),s{epochs(i)}{cells(1,1)}{cells(1,2)}.data(valid,3),'r.')
    
    %Plot valid spikes for CA1 cell
    spikerip = lookup(s{epochs(i)}{cells(2,1)}{cells(2,2)}.data(:,1),riptimes{day}{epochs(i)}.time);
    spikerip = riptimes{day}{epochs(i)}.nripples(spikerip);
    
    spikepos = lookup(s{epochs(i)}{cells(2,1)}{cells(2,2)}.data(:,1),p{epochs(i)}.data(:,1));
    spikepos = p{epochs(i)}.data(spikepos,5);

    valid = logical(spikerip==0) & logical(spikepos > 1)';
    plot(s{epochs(i)}{cells(2,1)}{cells(2,2)}.data(valid,2),s{epochs(i)}{cells(2,1)}{cells(2,2)}.data(valid,3),'b.')

    legend('Position','CA3','CA1')
    
    % Save figure
    %[y, m, d] = datevec(date);
    %savestring = sprintf('/home/mcarr/Figures/NovelObjectPaper/%d_%d_%d_example_ca1ca3_placefields_%d.pdf', m, d, y,epochs(i));
    %print('-dpdf', savestring)
end
