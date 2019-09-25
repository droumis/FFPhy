
ca1f = [];
ca3f = [];
ca1groups = [];
ca3groups = [];


%common filter attributes
%---------------------------------------------------------
ca1cellfilter = '(isequal($area, ''CA1'') && ($meanrate < 7))';
ca3cellfilter = '(isequal($area, ''CA3'') && ($meanrate < 7))';
timefilter = { {'getlinstate', '(($traj ~= -1) & (abs($velocity) >= 3))', 6}, {'getriptimes','($nripples == 0)', [], 'cellfilter', '((isequal($area, ''CA3''))|(isequal($area, ''CA1'')))'}};
iterator = 'singlecellanal';
filterfunction = {'calctotalmeanrate', {'spikes'}};
%--------------------------------------------------------------


j = 1;
animals = {'Bond','Frank','Nine'};
epochfilter = [];
for i = 1:14
    epochfilter{i} = ['($experimentday == ',num2str(i),') & ($dailyexposure == 1) & isequal($environment, ''TrackA'')'];
end
ca1f{j} = createfilter('animal',animals,'epochs',epochfilter,'cells',ca1cellfilter,'excludetime', timefilter,'iterator',iterator,'filterfunction',filterfunction);
ca3f{j} = createfilter('animal',animals,'epochs',epochfilter,'cells',ca3cellfilter,'excludetime', timefilter,'iterator',iterator,'filterfunction',filterfunction);
ca1f{j} = runfilter(ca1f{j});
ca3f{j} = runfilter(ca3f{j});
ca1groups{j} = numericgroupcombine(ca1f{j},0);
ca3groups{j} = numericgroupcombine(ca3f{j},0);
plotgroups(ca1groups{j},'c');
hold on
plotgroups(ca3groups{j},'m');

j = 2;
animals = {'Bond','Frank','Nine'};
epochfilter = [];
for i = 1:14
    epochfilter{i} = ['($experimentday == ',num2str(i),') & ($dailyexposure == 1) & isequal($environment, ''TrackB'')'];
end
ca1f{j} = createfilter('animal',animals,'epochs',epochfilter,'cells',ca1cellfilter,'excludetime', timefilter,'iterator',iterator,'filterfunction',filterfunction);
ca3f{j} = createfilter('animal',animals,'epochs',epochfilter,'cells',ca3cellfilter,'excludetime', timefilter,'iterator',iterator,'filterfunction',filterfunction);
ca1f{j} = runfilter(ca1f{j});
ca3f{j} = runfilter(ca3f{j});
ca1groups{j} = numericgroupcombine(ca1f{j},0);
ca3groups{j} = numericgroupcombine(ca3f{j},0);
plotgroups(ca1groups{j},'b');
plotgroups(ca3groups{j},'r');

j = 3;
animals = {'Dudley','Miles','Conley', 'Ten'};
epochfilter = [];
for i = 1:14
    epochfilter{i} = ['($experimentday == ',num2str(i),') & ($dailyexposure == 1) & isequal($environment, ''TrackA'')'];
end
ca1f{j} = createfilter('animal',animals,'epochs',epochfilter,'cells',ca1cellfilter,'excludetime', timefilter,'iterator',iterator,'filterfunction',filterfunction);
ca3f{j} = createfilter('animal',animals,'epochs',epochfilter,'cells',ca3cellfilter,'excludetime', timefilter,'iterator',iterator,'filterfunction',filterfunction);
ca1f{j} = runfilter(ca1f{j});
ca3f{j} = runfilter(ca3f{j});
ca1groups{j} = numericgroupcombine(ca1f{j},0);
ca3groups{j} = numericgroupcombine(ca3f{j},0);
plotgroups(ca1groups{j},'c');
plotgroups(ca3groups{j},'m');

j = 4;
animals = {'Dudley','Miles','Conley', 'Ten'};
epochfilter = [];
for i = 1:14
    epochfilter{i} = ['($experimentday == ',num2str(i),') & ($dailyexposure == 1) & isequal($environment, ''TrackB'')'];
end
ca1f{j} = createfilter('animal',animals,'epochs',epochfilter,'cells',ca1cellfilter,'excludetime', timefilter,'iterator',iterator,'filterfunction',filterfunction);
ca3f{j} = createfilter('animal',animals,'epochs',epochfilter,'cells',ca3cellfilter,'excludetime', timefilter,'iterator',iterator,'filterfunction',filterfunction);
ca1f{j} = runfilter(ca1f{j});
ca3f{j} = runfilter(ca3f{j});
ca1groups{j} = numericgroupcombine(ca1f{j},0);
ca3groups{j} = numericgroupcombine(ca3f{j},0);
plotgroups(ca1groups{j},'b');
plotgroups(ca3groups{j},'r');











