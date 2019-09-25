%% RUN FILTER TO CALCULATE OCCUPANCY NORMALIZED PLACE FIELDS

% Animal Selection
animals = {'CorrianderNO','Cyclops','Dunphy','Fafnir','Godzilla','Grendel'};

% Epoch selection
epochfilter = [];
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
f = setfilterfunction(f, 'calcopenfieldfiringrates_new', {'spikes','pos','task'});
f = runfilter(f);

%RUN FOR CA3
% Cell Filter
cellfilter = 'isequal($area, ''CA3'') & $meanrate < 7 & $numspikes > 1';
iterator = 'multicellanal';
% 
g = createfilter('animal',animals,'epochs',epochfilter,'cells',cellfilter,'excludetimefilter', timefilter, 'iterator', iterator);
g = setfilterfunction(g, 'calcopenfieldfiringrates_new', {'spikes','pos','task'});
g = runfilter(g);

%% GET ALL RIPPLE EVENTS FOR SLEEPS

% Animal Selection
animals = {'CorrianderNO','Cyclops','Dunphy','Fafnir','Godzilla','Grendel'};

% Epoch selection
epochfilter{1} = 'isequal($type,''sleep'') & isequal($runbefore,''baseline'')';
epochfilter{2} = 'isequal($type,''sleep'') & isequal($runbefore,''familiar'')';
epochfilter{3} = 'isequal($type,''sleep'') & isequal($runbefore,''novel'')';
epochfilter{4} = 'isequal($type,''sleep'') & isequal($runbefore,''supernovel'')';

%Define time filter
timefilter = {{'get2dstate', '$velocity < 4'}};
cellfilter = 'isequal($area, ''CA1'') & $meanrate < 7 & $numspikes > 1';

decodefilter = createfilter('animal',animals,'epochs',epochfilter,'cells',cellfilter,'excludetime', timefilter,'iterator',iterator);
decodefilter = setfilterfunction(decodefilter, 'getSWRactivity_sleep', {'spikes','pos','ripples','cellinfo'},'cellfilter','isequal($area,''CA1'')');
decodefilter = runfilter(decodefilter);

cellfilter = 'isequal($area, ''CA3'') & $meanrate < 7 & $numspikes > 1';

decodefilter3 = createfilter('animal',animals,'epochs',epochfilter,'cells',cellfilter,'excludetime', timefilter,'iterator',iterator);
decodefilter3 = setfilterfunction(decodefilter3, 'getSWRactivity_sleep', {'spikes','pos','ripples','cellinfo'},'cellfilter','isequal($area,''CA1'')');
decodefilter3 = runfilter(decodefilter3);

%% INITIALIZE VARIABLES

ca1 = calcsinglecellreactivationstats_sleep(f,decodefilter); 
ca3 = calcsinglecellreactivationstats_sleep(g,decodefilter3);

%% PLOT ACTIVATION PROBABILITY CA1

%Initialize variables
bas = ca1.activationb(~isnan(ca1.activationb));
fam = ca1.activationf(~isnan(ca1.activationf));
nov = ca1.activationn(~isnan(ca1.activationn));
sup = ca1.activations(~isnan(ca1.activations));

figure
hold on
bar([1 2 3 4],[mean(bas) mean(fam) mean(nov) mean(sup)],'b')
errorbar2([1 2 3 4],[mean(bas) mean(fam) mean(nov) mean(sup)],...
    [stderr(bas) stderr(fam) stderr(nov) stderr(sup)],'k')
set(gca,'xtick',1:4,'xticklabel',[{'Sleep Before'},{'Familiar'},{'Novel'},{'SuperNovel'}])
ylabel('Activation probability per SWR')

%Use kruskalwalis anova to test differences
X = [bas; fam; nov; sup];
group = [ones(size(bas)); 2*ones(size(fam)); 3*ones(size(nov)); 4*ones(size(sup))];
[p table stats] = anova1(X,group);
a = multcompare(stats);


%% PLOT ACTIVATION PROBABILITY CA3

%Initialize variables
bas = ca3.activationb(~isnan(ca3.activationb));
fam = ca3.activationf(~isnan(ca3.activationf));
nov = ca3.activationn(~isnan(ca3.activationn));
sup = ca3.activations(~isnan(ca3.activations));

figure
hold on
bar([1 2 3 4],[mean(bas) mean(fam) mean(nov) mean(sup)],'b')
errorbar2([1 2 3 4],[mean(bas) mean(fam) mean(nov) mean(sup)],...
    [stderr(bas) stderr(fam) stderr(nov) stderr(sup)],'k')
set(gca,'xtick',1:4,'xticklabel',[{'Sleep Before'},{'Familiar'},{'Novel'},{'SuperNovel'}])
ylabel('Activation probability per SWR')

%Use kruskalwalis anova to test differences
X = [bas; fam; nov; sup];
group = [ones(size(bas)); 2*ones(size(fam)); 3*ones(size(nov)); 4*ones(size(sup))];
[p table stats] = anova1(X,group);
a = multcompare(stats);

%% ACTIVATION PROBABILITY BY QUADRANT

%Initialize variables
nov1 = [ca1.activationn(ca1.typen==1)./ca1.activationb(ca1.typen==1); ca1.activations(ca1.types==1)./ca1.activationb(ca1.types==1)];
fam1 = [ca1.activationn(ca1.typen==2)./ca1.activationb(ca1.typen==2); ca1.activations(ca1.types==2)./ca1.activationb(ca1.types==2)];
nov1(isnan(nov1) | isinf(nov1)) = []; fam1(isnan(fam1) | isinf(fam1)) = [];
nov3 = [ca3.activationn(ca3.typen==1)./ca3.activationb(ca3.typen==1); ca3.activations(ca3.types==1)./ca3.activationb(ca3.types==1)];
fam3 = [ca3.activationn(ca3.typen==2)./ca3.activationb(ca3.typen==2); ca3.activations(ca3.types==2)./ca3.activationb(ca3.types==2)];
nov3(isnan(nov3) | isinf(nov3)) = []; fam3(isnan(fam3) | isinf(fam3)) = [];

figure
hold on
bar([1 2 5 6],[mean(fam1) mean(nov1) mean(fam3) mean(nov3)],'b')
errorbar2([1 2 5 6],[mean(fam1) mean(nov1) mean(fam3) mean(nov3)],...
    [stderr(fam1) stderr(nov1) stderr(fam3) stderr(nov3)],'k')
set(gca,'xtick',1:2,'xticklabel',[{'FamiliarObject'},{'NovelObject'}])
ylabel('Activation probability per SWR')

%Use kruskalwalis anova to test differences
X = [fam1; nov1; fam3; nov3];
group = [ones(size(fam1)); 2*ones(size(nov1)); 3*ones(size(fam1)); 4*ones(size(nov3))];
[p table stats] = kruskalwallis(X,group);
a = multcompare(stats);