%% RUN FILTER TO CALCULATE OCCUPANCY NORMALIZED PLACE FIELDS

% Animal Selection
animals = {'CorrianderNO','Cyclops','Dunphy','Fafnir','Godzilla','Grendel'};

% Epoch selection
epochfilter{1} = 'isequal($type,''run'') & isequal($session,''familiar'')';
epochfilter{2} = 'isequal($type,''run'') & isequal($session,''novel'')';
epochfilter{3} = 'isequal($type,''run'') & isequal($session,''supernovel'')';

%Define time filter
timefilter = {{'get2dstate', '($velocity > 1)'},{'getriptimes', '($nripples == 0)', [], 'cellfilter', 'isequal($area, ''CA1'')'}};

%RUN FOR CA1
% Cell Filter
cellfilter = 'isequal($area, ''CA1'') & $meanrate < 7 & $numspikes > 1';
iterator = 'multicellanal';

f = createfilter('animal',animals,'epochs',epochfilter,'cells',cellfilter,'excludetimefilter', timefilter, 'iterator', iterator);
f = setfilterfunction(f, 'calcopenfieldfiringrates', {'spikes','pos','task'});
f = runfilter(f);

%RUN FOR CA3
% Cell Filter
cellfilter = 'isequal($area, ''CA3'') & $meanrate < 7 & $numspikes > 1';
iterator = 'multicellanal';
% 
g = createfilter('animal',animals,'epochs',epochfilter,'cells',cellfilter,'excludetimefilter', timefilter, 'iterator', iterator);
g = setfilterfunction(g, 'calcopenfieldfiringrates', {'spikes','pos','task'});
g = runfilter(g);

%% GET ALL RIPPLE EVENTS FOR RUNS

timefilter = {{'get2dstate', '$velocity < 4'}};
cellfilter = 'isequal($area, ''CA1'') & $meanrate < 7 & $numspikes > 1';

decodefilter = createfilter('animal',animals,'epochs',epochfilter,'cells',cellfilter,'excludetime', timefilter,'iterator',iterator);
decodefilter = setfilterfunction(decodefilter, 'getSWRactivity', {'spikes','task','ripples','cellinfo'},'cellfilter','isequal($area,''CA1'')');
decodefilter = runfilter(decodefilter);

cellfilter = 'isequal($area, ''CA3'') & $meanrate < 7 & $numspikes > 1';

decodefilter3 = createfilter('animal',animals,'epochs',epochfilter,'cells',cellfilter,'excludetime', timefilter,'iterator',iterator);
decodefilter3 = setfilterfunction(decodefilter3, 'getSWRactivity', {'spikes','task','ripples','cellinfo'},'cellfilter','isequal($area,''CA1'')');
decodefilter3 = runfilter(decodefilter3);

%% INITIALIZE VARIABLES

ca1 = calcsinglecellreactivationstats(f,decodefilter); 
ca3 = calcsinglecellreactivationstats(g,decodefilter3);

%% ACTIVATION PROBABILITY CA1: NOVEL vs. FAMILIAR, PLACE CELL vs. NON PLACE FIELD

%Initialize variables
fam = ca1.activationf(~isnan(ca1.activationf) & ca1.peakf > 3 & ca1.fmean > 0.1);
nov = ca1.activationn(~isnan(ca1.activationn) & ca1.peakn > 3 & ca1.nmean > 0.1);
sup = ca1.activations(~isnan(ca1.activations) & ca1.peaks > 3 & ca1.smean > 0.1);
famn = ca1.activationf(~isnan(ca1.activationf) & (ca1.peakf < 3 | ca1.fmean < 0.1));
novn = ca1.activationn(~isnan(ca1.activationn) & (ca1.peakn < 3 | ca1.nmean < 0.1));
supn = ca1.activations(~isnan(ca1.activations) & (ca1.peaks < 3 | ca1.smean < 0.1));

figure
hold on
bar([1 5],[mean(fam) mean(famn)],'b')
bar([2 6],[mean(nov) mean(novn)],'c')
bar([3 7],[mean(sup) mean(supn)],'m')
legend('Familiarization','Novel Location','Novel Object & Location')
errorbar2([1 2 3 5 6 7],[mean(fam) mean(nov) mean(sup) mean(famn) mean(novn) mean(supn)],...
    [stderr(fam) stderr(nov) stderr(sup) stderr(famn) stderr(novn) stderr(supn)],'k')
set(gca,'xtick',[2 6],'xticklabel',[{'cells with place fields'},{'cells without place fields'}])
ylabel('Activation probability per SWR')

% 2 way Anova
X = [fam; nov; sup; famn; novn; supn];
group1 = [ones(size(fam)); 2*ones(size(nov)); 3*ones(size(sup)); ones(size(famn)); 2*ones(size(novn)); 3*ones(size(supn))];
group2 = [ones(size(fam)); ones(size(nov)); ones(size(sup)); 2*ones(size(famn)); 2*ones(size(novn)); 2*ones(size(supn))];
anovan(X,{group1 group2},'model','full')
%Main effect of novelty (p<0.05) and main effect of place cell (p<1e-5)

%Kruskalwallis anova
X = [fam; nov; sup];
group = [ones(size(fam)); 2*ones(size(nov)); 3*ones(size(sup))];
[p table stats] = kruskalwallis(X,group);
a = multcompare(stats);
%Main effect of novelty p<0.001, fam vs. nov p<0.05, fam vs. sup p<0.01

%Rank sum
p = ranksum(fam,nov); %p<0.01
p = ranksum(fam,sup); %p<0.001
p = ranksum(famn,fam); %p<1e-5
p = ranksum(novn,nov); %p<1e-5
p = ranksum(supn,sup); %p<1e-5

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/NovelObjectPaper/%d_%d_%d_singlecell_reactivation.pdf', m, d, y);
print('-dpdf', savestring)


% Compare average change within cell
nov = ca1.activationn(ca1.peakn > 3 & ca1.nmean > 0.1) - ca1.activationf(ca1.peakn > 3 & ca1.nmean > 0.1);
sup = ca1.activations(ca1.peaks > 3 & ca1.smean > 0.1) - ca1.activationf(ca1.peaks > 3 & ca1.smean > 0.1);
invalid = isnan(nov)|isinf(nov); nov(invalid) = []; invalid = isnan(sup)|isinf(sup); sup(invalid) = [];

figure
hold on
bar(1,mean(nov),'c')
bar(2,mean(sup),'m')
errorbar2([1 2],[mean(nov) mean(sup)],[stderr(nov) stderr(sup)],'k')
set(gca,'xtick',[1 2],'xticklabel',[{'Novel Location'},{'Novel Object & Location'}])
ylabel('Change in activation probabiltiy per SWR')

%Sign rank
p = signrank(nov,0); %p<0.001
p = signrank(sup,0); %p<1e-5

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/NovelObjectPaper/%d_%d_%d_singlecell_ratio_reactivation.pdf', m, d, y);
print('-dpdf', savestring)

%% ACTIVATION PROBABILITY CA3: NOVEL vs. FAMILIAR, PLACE CELL vs. NON PLACE FIELD

%Initialize variables
fam = ca3.activationf(~isnan(ca3.activationf) & ca3.peakf > 3 & ca3.fmean > 0.1);
nov = ca3.activationn(~isnan(ca3.activationn) & ca3.peakn > 3 & ca3.nmean > 0.1);
sup = ca3.activations(~isnan(ca3.activations) & ca3.peaks > 3 & ca3.smean > 0.1);
famn = ca3.activationf(~isnan(ca3.activationf) & (ca3.peakf < 3 | ca3.fmean < 0.1));
novn = ca3.activationn(~isnan(ca3.activationn) & (ca3.peakn < 3 | ca3.nmean < 0.1));
supn = ca3.activations(~isnan(ca3.activations) & (ca3.peaks < 3 | ca3.smean < 0.1));

figure
hold on
bar([1 5],[mean(fam) mean(famn)],'b')
bar([2 6],[mean(nov) mean(novn)],'c')
bar([3 7],[mean(sup) mean(supn)],'m')
legend('Familiarization','Novel Location','Novel Object & Location')
errorbar2([1 2 3 5 6 7],[mean(fam) mean(nov) mean(sup) mean(famn) mean(novn) mean(supn)],...
    [stderr(fam) stderr(nov) stderr(sup) stderr(famn) stderr(novn) stderr(supn)],'k')
set(gca,'xtick',[2 6],'xticklabel',[{'cells with place fields'},{'cells without place fields'}],'ylim',[0 0.2])
ylabel('Activation probability per SWR')

% 2 way Anova
X = [fam; nov; sup; famn; novn; supn];
group1 = [ones(size(fam)); 2*ones(size(nov)); 3*ones(size(sup)); ones(size(famn)); 2*ones(size(novn)); 3*ones(size(supn))];
group2 = [ones(size(fam)); ones(size(nov)); ones(size(sup)); 2*ones(size(famn)); 2*ones(size(novn)); 2*ones(size(supn))];
anovan(X,{group1 group2},'model','full')
%main effect of place cell (p<1e-5), no main effect of novelty

%Kruskalwallis anova
X = [fam; nov; sup];
group = [ones(size(fam)); 2*ones(size(nov)); 3*ones(size(sup))];
[p table stats] = kruskalwallis(X,group);
a = multcompare(stats);
%No main effect of novelty p>0.06

%Rank sum
p = ranksum(fam,nov); %p<0.03
p = ranksum(fam,sup); %p>0.05
p = ranksum(famn,fam); %p<0.001
p = ranksum(novn,nov); %p<1e-5
p = ranksum(supn,sup); %p<0.001

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/NovelObjectPaper/%d_%d_%d_CA3_singlecell_reactivation.pdf', m, d, y);
print('-dpdf', savestring)

% Compare average change within cell
nov = ca3.activationn(ca3.peakn > 3 & ca3.nmean > 0.1) - ca3.activationf(ca3.peakn > 3 & ca3.nmean > 0.1);
sup = ca3.activations(ca3.peaks > 3 & ca3.smean > 0.1) - ca3.activationf(ca3.peaks > 3 & ca3.smean > 0.1);
invalid = isnan(nov)|isinf(nov); nov(invalid) = []; invalid = isnan(sup)|isinf(sup); sup(invalid) = [];

figure
hold on
bar(1,mean(nov),'c')
bar(2,mean(sup),'m')
errorbar2([1 2],[mean(nov) mean(sup)],[stderr(nov) stderr(sup)],'k')
set(gca,'xtick',[1 2],'xticklabel',[{'Novel Location'},{'Novel Object & Location'}])
ylabel('Change in activation probabiltiy per SWR')
set(gca,'ylim',[-0.01 0.1])
%Sign rank
p = signrank(nov,0); %p>0.2
p = signrank(sup,0); %p>0.1

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/NovelObjectPaper/%d_%d_%d_CA3_singlecell_ratio_reactivation.pdf', m, d, y);
print('-dpdf', savestring)

%% ACTIVATION PROBABILITY BY QUADRANT: CA1

%Initialize variables
nov = ca1.activationn(~isnan(ca1.activationn) & ca1.peakn > 3 & ca1.nmean > 0.1 & ca1.typen == 1)- ca1.activationf(~isnan(ca1.activationn) & ca1.peakn > 3 & ca1.nmean > 0.1 & ca1.typen == 1);
famn = ca1.activationn(~isnan(ca1.activationn) & ca1.peakn > 3 & ca1.nmean > 0.1 & ca1.typen == 2) -  ca1.activationf(~isnan(ca1.activationn) & ca1.peakn > 3 & ca1.nmean > 0.1 & ca1.typen == 2);
sup = ca1.activations(~isnan(ca1.activations) & ca1.peaks > 3 & ca1.smean > 0.1 & ca1.types == 1) - ca1.activationf(~isnan(ca1.activations) & ca1.peaks > 3 & ca1.smean > 0.1 & ca1.types == 1);
fams = ca1.activations(~isnan(ca1.activations) & ca1.peaks > 3 & ca1.smean > 0.1 & ca1.types == 2) - ca1.activationf(~isnan(ca1.activations) & ca1.peaks > 3 & ca1.smean > 0.1 & ca1.types == 2);
invalid = isnan(nov); nov(invalid) = []; invalid = isnan(famn); famn(invalid) = [];
invalid = isnan(sup); sup(invalid) = []; invalid = isnan(fams); fams(invalid) = [];

figure
hold on
bar([1 2],[mean(nov) mean(famn)],'c')
bar([4 5],[mean(sup) mean(fams)],'m')
errorbar2([1 2 4 5],[mean(nov) mean(famn) mean(sup) mean(fams)],...
    [stderr(nov) stderr(famn) stderr(sup) stderr(fams)],'k')
set(gca,'xtick',[1 2 4 5],'xticklabel',[{'Novel Object'},{'Familiar Object'}],'ylim',[-0.02 0.08])
ylabel('Change in activation probability per SWR')

%Rank sum
p = ranksum(nov,famn); %p>0.5
p = ranksum(sup,fams); %p>0.6

%No difference in single cell reactivation for the Novel Quadrant vs. Familiar quadrant
% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/NovelObjectPaper/%d_%d_%d_singlecell_reactivation_quadrant.pdf', m, d, y);
print('-dpdf', savestring)

%% ACTIVATION PROBABILITY BY QUADRANT: CA3

%Initialize variables
nov = ca3.activationn(~isnan(ca3.activationn) & ca3.peakn > 3 & ca3.nmean > 0.1 & ca3.typen == 1)- ca3.activationf(~isnan(ca3.activationn) & ca3.peakn > 3 & ca3.nmean > 0.1 & ca3.typen == 1);
famn = ca3.activationn(~isnan(ca3.activationn) & ca3.peakn > 3 & ca3.nmean > 0.1 & ca3.typen == 2) -  ca3.activationf(~isnan(ca3.activationn) & ca3.peakn > 3 & ca3.nmean > 0.1 & ca3.typen == 2);
sup = ca3.activations(~isnan(ca3.activations) & ca3.peaks > 3 & ca3.smean > 0.1 & ca3.types == 1) - ca3.activationf(~isnan(ca3.activations) & ca3.peaks > 3 & ca3.smean > 0.1 & ca3.types == 1);
fams = ca3.activations(~isnan(ca3.activations) & ca3.peaks > 3 & ca3.smean > 0.1 & ca3.types == 2) - ca3.activationf(~isnan(ca3.activations) & ca3.peaks > 3 & ca3.smean > 0.1 & ca3.types == 2);
invalid = isnan(nov); nov(invalid) = []; invalid = isnan(famn); famn(invalid) = [];
invalid = isnan(sup); sup(invalid) = []; invalid = isnan(fams); fams(invalid) = [];

figure
hold on
bar([1 2],[mean(famn) mean(nov)],'c')
bar([4 5],[mean(fams) mean(sup)],'m')
errorbar2([1 2 4 5],[mean(famn) mean(nov) mean(fams) mean(sup)],...
    [stderr(famn) stderr(nov) stderr(fams) stderr(sup)],'k')
set(gca,'xtick',[1 2 4 5],'xticklabel',[{'Familiar Quadrant'},{'Novel Quadrant'}])
ylabel('Change in activation probability per SWR')

%Rank sum
p = ranksum(nov,famn); %p>0.7
p = ranksum(sup,fams); %p>0.4

%No difference in single cell reactivation for the Novel Quadrant vs. Familiar quadrant
% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/NovelObjectPaper/%d_%d_%d_CA3_singlecell_reactivation_quadrant.pdf', m, d, y);
print('-dpdf', savestring)