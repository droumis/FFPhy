
%Animal selection
%-----------------------------------------------------
%animals = {'Dudley','Miles','Conley','Bond','Frank','Nine','Ten'};
%animals = {'Dudley','Miles','Conley','Alex','Bond','Frank','Nine'};
%animals = {'Ten'};
%animals = {'Dudley','Miles','Conley', 'Ten'};
animals = {'Bond','Frank','Nine'};
%animals = {'Bond'};
%-----------------------------------------------------


%Filter creation
%--------------------------------------------------------
epochfilter = [];
for i = 1:14
    %epochfilter{i} = ['($exposureday == ',num2str(i),') & ($dailyexposure == 1)'];
    epochfilter{i} = ['($experimentday == ',num2str(i),') & ($dailyexposure == 1) & isequal($environment, ''TrackA'')'];
end

ca1cellfilter = '(isequal($area, ''CA1'') && ($meanrate < 7))';
ca3cellfilter = '(isequal($area, ''CA3'') && ($meanrate < 7))';

timefilter = { {'getlinstate', '(($traj ~= -1) & (abs($velocity) >= 3))', 6} };
%timefilter = {};

iterator = 'singlecellanal';
%iterator = 'multicellanal';

ca1f = createfilter('animal',animals,'epochs',epochfilter,'cells',ca1cellfilter,'excludetime', timefilter, 'iterator', iterator);
ca3f = createfilter('animal',animals,'epochs',epochfilter,'cells',ca3cellfilter,'excludetime', timefilter, 'iterator', iterator);
%-----------------------------------------------------------




%non occ-normalized mean rate
%--------------------------------------------
figure
ca1f = setfilterfunction(ca1f, 'calctotalmeanrate', {'spikes'});
ca3f = setfilterfunction(ca3f, 'calctotalmeanrate', {'spikes'});

ca1f = runfilter(ca1f);
ca3f = runfilter(ca3f);

ca1groups = numericgroupcombine(ca1f);
ca3groups = numericgroupcombine(ca3f);

plotgroups(ca1groups);
hold on
plotgroups(ca3groups,'r');
%-------------------------------------------------

%occ normalized mean rate
%---------------------------------------------------
figure
ca1f = setfilterfunction(ca1f, 'calcoccnormmeanrate', {'spikes', 'linpos'});
ca3f = setfilterfunction(ca3f, 'calcoccnormmeanrate', {'spikes', 'linpos'});
disp('ca1');
ca1f = runfilter(ca1f);
disp('ca3');
ca3f = runfilter(ca3f);

ca1groups = numericgroupcombine(ca1f);
ca3groups = numericgroupcombine(ca3f);

plotgroups(ca1groups);
hold on
plotgroups(ca3groups,'r');
%----------------------------------------------------

%population mean rate
%---------------------------------------------------
figure
ca1f = setfilterfunction(ca1f, 'calcpopulationrate', {'spikes'},60);
ca3f = setfilterfunction(ca3f, 'calcpopulationrate', {'spikes'},60);
disp('ca1');
ca1f = runfilter(ca1f);
disp('ca3');
ca3f = runfilter(ca3f);

ca1groups = numericgroupcombine(ca1f);
ca3groups = numericgroupcombine(ca3f);

plotgroups(ca1groups);
hold on
plotgroups(ca3groups,'r');
%----------------------------------------------------

