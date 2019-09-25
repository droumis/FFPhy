%% RUN FILTER TO CALCULATE OCCUPANCY NORMALIZED PLACE FIELDS

% Animal Selection% Animal Selection
animals = {'CorrianderNO','Cyclops','Dunphy','Fafnir','Godzilla','Grendel'};

% Epoch selection
epochfilter{1} = 'isequal($type,''run'') & isequal($session,''familiar'')';
epochfilter{2} = 'isequal($type,''run'') & isequal($session,''novel'')';
epochfilter{3} = 'isequal($type,''run'') & isequal($session,''supernovel'')';

%Define time filter
timefilter = {{'get2dstate', '($velocity < 4)'}};

%RUN FOR CA1
% Cell Filter
cellfilter = 'isequal($area, ''CA1'')|isequal($area,''CA3'') & $numspikes > 100';

iterator = 'multicellanal';

f = createfilter('animal',animals,'epochs',epochfilter,'cells',cellfilter,'excludetimefilter', timefilter, 'iterator', iterator);
f = setfilterfunction(f, 'getcoactiveripples_and_overlap', {'ripples','pos','task','spikes','cellinfo'},'binsize',2);
f = runfilter(f);

%% Initialize variables

g= numericgroupcombine(f, 1); 

% First column: animal
% Second column: overlap
% Third colunn: coactive ripples
% Fourth column: coactive z-score
% Fifth column: number of ripples
% Sixth column: joint surprise
% Seventh column: 1st cell's peak quadrant
% Eighth column: 2nd cell's peak quadrant

%% Coactivation probability per SWR
fam = g{1}(:,3); nov = g{2}(:,3); sup = g{3}(:,3);
invalid = isnan(fam); fam(invalid) = [];
invalid = isnan(nov); nov(invalid) = [];
invalid = isnan(sup); sup(invalid) = [];

figure
hold on
bar(1:3,[mean(fam) mean(nov) mean(sup)],'b');
errorbar2(1:3,[mean(fam) mean(nov) mean(sup)],[std(fam)./sqrt(length(fam)-1) std(nov)./sqrt(length(nov)-1) std(sup)./sqrt(length(sup)-1)],'k')
set(gca,'xtick',1:3,'xticklabel',[{'Familiar Session'},{'Novel Session'},{'Super Novel Session'}])
ylabel('Coactivation probability per SWR')

%Use kruskalwalis anova to test differences
X = [fam; nov; sup];
group = [ones(size(fam)); 2*ones(size(nov)); 3*ones(size(sup))];
[p table stats] = kruskalwallis(X,group);
a = multcompare(stats);

%Fam is different than nov and sup at p<1e-5
%Nov and sup are different at p<1e-5

% % Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/NovelObjectPaper/%d_%d_%d_coactivation_probability.pdf', m, d, y);
print('-dpdf', savestring)

%% Coactivity z-score

fam = g{1}(:,4); nov = g{2}(:,4); sup = g{3}(:,4);
invalid = isnan(fam); fam(invalid) = [];
invalid = isnan(nov); nov(invalid) = [];
invalid = isnan(sup); sup(invalid) = [];

figure
hold on
bar(1:3,[mean(fam) mean(nov) mean(sup)],'b');
errorbar2(1:3,[mean(fam) mean(nov) mean(sup)],[std(fam)./sqrt(length(fam)-1) std(nov)./sqrt(length(nov)-1) std(sup)./sqrt(length(sup)-1)],'k')
set(gca,'xtick',1:3,'xticklabel',[{'Familiar Session'},{'Novel Session'},{'Super Novel Session'}])
ylabel('Coactivity z-score')

%Use kruskalwalis anova to test differences
X = [fam; nov; sup];
group = [ones(size(fam)); 2*ones(size(nov)); 3*ones(size(sup))];
[p table stats] = kruskalwallis(X,group);
a = multcompare(stats);

%Fam is different than nov and sup at p<0.01
%Nov and sup are different at p<0.05

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/NovelObjectPaper/%d_%d_%d_coactivity_zscore.pdf', m, d, y);
print('-dpdf', savestring)

%% Is enhancement during novelty due to coordinated reactivation of novel object cells?

%No difference in coactivation z-score for cell pairs with containing a
%novel obj. vs. a familiar object
nov = [g{2}(:,[3 7 8])];%; g{3}(:,[4 7 8])];
invalid = isnan(nov(:,1)); nov(invalid,:) = [];

FAM = nov(nov(:,2) == 2 | nov(:,3)==2,1);
NOV = nov(nov(:,2) == 1 | nov(:,3)==1,1);
p = ranksum(FAM,NOV);
%pvalue > 0.26

%There is more reactivation of cells with place fields in the novel object
%quadrant than cells with place fields in the familiar object quadrant
FAM = nov(nov(:,2) == 2 | nov(:,3)==2,1);
%length(FAM)./size(nov,1) = 0.29

FAM = nov(nov(:,2) == 2 & nov(:,3)==2,1);
%length(FAM)./size(nov,1) = 0.06

NOV = nov(nov(:,2) == 1 & nov(:,3)==1,1);
%length(NOV)./size(nov,1) = 0.11
ztestprop([length(NOV) size(nov,1)-length(NOV)],1/16)


%% IS COACTIVATION CONSISTANT WITH REPLAY?

fam = g{1}(:,[2 4]); nov = g{2}(:,[2 4]); sup = g{3}(:,[2 4]);
invalid = isnan(fam(:,1))|isnan(fam(:,2)); fam(invalid,:) = [];
invalid = isnan(nov(:,1))|isnan(nov(:,2)); nov(invalid,:) = [];
invalid = isnan(sup(:,1))|isnan(sup(:,2)); sup(invalid,:) = [];

figure
hold on
bin = [0.001 0.2 0.4 0.6];
subs = lookup(fam(:,1),bin,1);
af = accumarray(subs,fam(:,2),[length(bin) 1],@(x) mean(x));
afe= accumarray(subs,fam(:,2),[length(bin) 1],@(x) std(x)./sqrt(length(x)-1));
subs = lookup(nov(:,1),bin,1);
an = accumarray(subs,nov(:,2),[length(bin) 1],@(x) mean(x));
ane= accumarray(subs,nov(:,2),[length(bin) 1],@(x) std(x)./sqrt(length(x)-1));
subs = lookup(sup(:,1),bin,1);
as = accumarray(subs,sup(:,2),[length(bin) 1],@(x) mean(x));
ase= accumarray(subs,sup(:,2),[length(bin) 1],@(x) std(x)./sqrt(length(x)-1));

bar(bin,af,'b')
bar(bin+1,an,'r')
bar(bin+2,as,'c')
errorbar2(bin,af,afe,'k')
errorbar2(bin+1,an,ane,'k')
errorbar2(bin+2,as,ase,'k')

bin(1) = 0;
set(gca,'xtick',[bin bin+1 bin+2],'xtickLabel',[bin bin bin])
ylabel('Coactivity z-score')
xlabel('Place field overlap')
legend('Familiar','Novel','SuperNovel')
box off
% No difference between sessions in the relationship between coactivation
% and overlab as a function of novelty
[r p] = corrcoef(fam(:,1),fam(:,2)); %p <1e-5
[r p] = corrcoef(nov(:,1),nov(:,2)); %p <1e-5
[r p] = corrcoef(sup(:,1),sup(:,2)); %p <1e-5

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/NovelObjectPaper/%d_%d_%d_coactivity_zscore_vs_placefieldoverlap.pdf', m, d, y);
print('-dpdf', savestring)
