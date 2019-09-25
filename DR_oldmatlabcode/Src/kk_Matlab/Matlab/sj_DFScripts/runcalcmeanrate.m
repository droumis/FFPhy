%Animal selection
animals = {'Bond','Conley','Corriander','Dudley','Eight','Five','Frank','Miles','Ten'};

%Filter creation
epochfilter = [];
for i = 1:14
    %epochfilter{i} = ['(isequal($type,''sleep'')) & ($sleepnum == 2) & ($runbefore.exposureday == ',num2str(i),') & ($runbefore.dailyexposure == 1)'];  
    epochfilter{i} = ['($exposureday == ',num2str(i),') & ($dailyexposure == 1)'];
end

ca1cellfilter = '(isequal($area, ''CA1'') && ($meanrate < 7))';
ca3cellfilter = '(isequal($area, ''CA3'') && ($meanrate < 7))';

%timefilter = { {'getlinstate', 'abs($velocity) >= 3)', 6} };
%timefilter = { {'getlinstate', 'abs($velocity) >= 3)', 6},{'getriptimes','($nripples == 0)', [], 'cellfilter', '((isequal($area, ''CA3''))|(isequal($area, ''CA1'')))'}};
%timefilter = { {'getlinstate', 'abs($velocity) < 3)', 6}, {'getriptimes','($nripples > 2)', [], 'cellfilter', '((isequal($area, ''CA3''))|(isequal($area, ''CA1'')))'}};
%timefilter = {{'getriptimes', '($nripples >= 1)', [], 'cellfilter', '(isequal($area, ''CA1''))|(isequal($area, ''CA3''))'}};
timefilter = {{'getgammatimes', '($ngamma == 0)', [], 'high','tetfilter', 'isequal($area,''CA1'')','minthresh',2}, ...
     {'getgammatimes', '($ngamma == 0)', [], 'low','tetfilter', 'isequal($area,''CA1'')','minthresh',2}};

ca1f = createfilter('animal',animals,'epochs',epochfilter,'cells',ca1cellfilter,'excludetime', timefilter);
ca3f = createfilter('animal',animals,'epochs',epochfilter,'cells',ca3cellfilter,'excludetime', timefilter);

%run function- single cells
iterator = 'singlecellanal';

ca1f = setfilteriterator(ca1f,iterator);
ca3f = setfilteriterator(ca3f,iterator);

ca1f = setfilterfunction(ca1f, 'calctotalmeanrate', {'spikes'});
ca3f = setfilterfunction(ca3f, 'calctotalmeanrate', {'spikes'});

ca1f = runfilter(ca1f);
ca3f = runfilter(ca3f);

%plot
ca1groups = numericgroupcombine(ca1f,0);
ca3groups = numericgroupcombine(ca3f,0);

figure
plotgroups(ca1groups);
hold on
plotgroups(ca3groups,'r'); 