%Animal selection
%-----------------------------------------------------
%animals = {'Dudley','Miles','Conley','Bond','Frank','Nine','Ten'};
%animals = {'Miles','Conley','Ten','Bond','Frank','Nine'};
%animals = {'Dudley','Miles','Conley','Alex','Bond','Frank','Nine'};
%animals = {'Ten'};
%animals = {'Dudley','Miles','Conley', 'Ten'};
animals = {'Bond','Frank','Nine'};
%animals = {'Bond'};
%-----------------------------------------------------


%Filter creation
%--------------------------------------------------------
epochfilter = [];
epochfilter{1} = ['($dailyexposure == 1) & isequal($environment, ''TrackA'')'];
epochfilter{2} = ['($dailyexposure == 1) & isequal($environment, ''TrackB'')'];

%epochfilter{1} = ['($dailyexposure == 1) & isequal($environment, ''TrackB'') & ($experimentday > 3)'];
%epochfilter{2} = ['($dailyexposure == 2) & isequal($environment, ''TrackB'') & ($experimentday > 3)']; 

ca1cellfilter = '(isequal($area, ''CA1'') && ($meanrate < 7))';
ca3cellfilter = '(isequal($area, ''CA3'') && ($meanrate < 7))';

timefilter = { {'getlinstate', '(($traj ~= -1) & (abs($velocity) >= 3))', 6} };
%timefilter = {};

ca1f = createfilter('animal',animals,'epochs',epochfilter,'cells',ca1cellfilter,'excludetime', timefilter);
ca3f = createfilter('animal',animals,'epochs',epochfilter,'cells',ca3cellfilter,'excludetime', timefilter);
%-----------------------------------------------------------

%Pairwise comparison
%-----------------------------------------------------------
iterator = 'singlecellanal';

ca1f = setfilteriterator(ca1f,iterator);
ca3f = setfilteriterator(ca3f,iterator);

ca1f = setfilterfunction(ca1f, 'calctotalmeanrate', {'spikes'},'appendindex',1);
ca3f = setfilterfunction(ca3f, 'calctotalmeanrate', {'spikes'},'appendindex',1);

%ca1f = setfilterfunction(ca1f, 'funcSwitchBox', {'spikes','linpos'},'calctotalmeanrate','calcoutfieldfiring');
%ca3f = setfilterfunction(ca3f, 'funcSwitchBox', {'spikes','linpos'},'calctotalmeanrate','calcoutfieldfiring');

ca1f = runfilter(ca1f);
ca3f = runfilter(ca3f);

ca1groups = numericgroupcombine(ca1f,1);
ca3groups = numericgroupcombine(ca3f,1);

indexcolumns = [1 2 4 5];
ca1cmp = indexmatch(ca1groups{1},ca1groups{2},indexcolumns);
ca3cmp = indexmatch(ca3groups{1},ca3groups{2},indexcolumns);

%--------------------------------------------------------------


