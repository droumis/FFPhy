%Animal selection
%-----------------------------------------------------
animals = {'JYB'};


%Filter creation
%--------------------------------------------------------

% epoch filter
epochfilter{1} = ['isequal($environment, ''run1'')'];


% cell filter
MScellfilter = '(isequal($area, ''MS'') && ($meanrate < 100))';

% time filter
timefilter = { {'getlinstate', '(($traj ~= -1) & (abs($velocity) >= 0))', 6} };

% iterator
iterator = 'singlecellanal';

% filter creation
MSf = createfilter('animal',animals,'epochs',epochfilter,'cells',MScellfilter,'excludetime', timefilter, 'iterator', iterator);

% set analysis function
MSf = setfilterfunction(MSf, 'calctotalmeanrate', {'spikes'});

% run analysis
MSf = runfilter(MSf);
