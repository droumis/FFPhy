%Animal selection
%-----------------------------------------------------
%animals = {'Miles','Conley','Bond','Frank','Nine','Ten'};

animals = {'Miles','Nine','Ten'};
%animals = {'Conley','Bond','Frank'};


%animals = {'Bond'};
%animals = {'Miles','Conley','Ten','Bond','Frank','Nine'};
%animals = {'Dudley','Miles','Conley','Alex','Bond','Frank','Nine'};
%animals = {'Ten'};
%animals = {'Dudley','Miles','Conley', 'Ten'};
%animals = {'Bond','Frank','Nine'};
%animals = {'Bond'};
%-----------------------------------------------------


%Filter creation
%--------------------------------------------------------
epochfilter1 = [];
epochfilter2 = [];
epochfilter1{1} = ['($dailyexposure == 1) & isequal($environment, ''TrackB'') & ($exposureday <= 4)'];
epochfilter2{1} = ['($dailyexposure == 1) & isequal($environment, ''TrackA'')'];

cellfilter1 = '((isequal($area, ''CA3'')) && ($meanrate < 7) && ($meanrate < .1))';
cellfilter2 = '((isequal($area, ''CA3'')) && ($meanrate < 7) && ($meanrate > .1))';



% epochfilter1{1} = ['($dailyexposure == 1) & isequal($environment, ''TrackA'') & ($exposureday < 4)'];
% epochfilter2{1} = ['($dailyexposure == 1) & isequal($environment, ''TrackA'') & ($exposureday < 4)'];
% 
% cellfilter1 = '((isequal($area, ''CA1'')|isequal($area, ''CA3'')) && ($meanrate < 7) && ($meanrate > .1))';
% cellfilter2 = '((isequal($area, ''CA1'')|isequal($area, ''CA3'')) && ($meanrate < 7) && ($meanrate > .1))';


%timefilter1 = {{'getriptimes', '($nripples > 0)', [], 'cellfilter', '(isequal($area, ''CA1''))','maxcell',1,'minstd',3}};
%timefilter1 = { {'getlinstate', '(($traj ~= -1) & (abs($velocity) >= 3))',6} };
%timefilter2 = {{'getriptimes', '($nripples > 0)', [], 'cellfilter', '(isequal($area, ''CA1''))','maxcell',1,'minstd',3}};
%timefilter2 = { {'getlinstate', '(($traj ~= -1) & (abs($velocity) >= 3))',6} };

timefilter1 = {};
timefilter2 = {};

f1 = createfilter('animal',animals,'epochs',epochfilter1,'cells',cellfilter1,'excludetime', timefilter1);
f2 = createfilter('animal',animals,'epochs',epochfilter2,'cells',cellfilter2,'excludetime', timefilter2);

%-----------------------------------------------------------

%Pairwise comparison
%-----------------------------------------------------------
iterator = 'singlecellanal';

f1 = setfilteriterator(f1,iterator);
f2 = setfilteriterator(f2,iterator);

%f1 = setfilterfunction(f1, 'calctotalmeanrate', {'spikes'},'appendindex',1);
%f1 = setfilterfunction(f1, 'calcpeakrate', {'spikes','linpos'},'appendindex',1);
%f2 = setfilterfunction(f2, 'calctotalmeanrate', {'spikes'},'appendindex',1);
%f2 = setfilterfunction(f2, 'calcpeakrate', {'spikes','linpos'},'appendindex',1);

f1 = setfilterfunction(f1, 'calcrippleactivationprob', {'spikes','ripples','cellinfo'},'appendindex',1,'ratio',1);
f2 = setfilterfunction(f2, 'calcrippleactivationprob', {'spikes','ripples','cellinfo'},'appendindex',1,'ratio',1);

f1 = runfilter(f1);
f2 = runfilter(f2);

f1groups = numericgroupcombine(f1,1);
f2groups = numericgroupcombine(f2,1);

indexcolumns = [1 2 4 5];
cmp = indexmatch(f1groups{1},f2groups{1},indexcolumns);


%--------------------------------------------------------------


