%Animal selection
%-----------------------------------------------------
animals = {'Miles','Frank','Nine'};


%Filter creation
%--------------------------------------------------------

% epoch filter
epochfilter{1} = ['isequal($environment, ''TrackA'')'];
epochfilter{2} = ['isequal($environment, ''TrackB'')'];


% cell filter
ca1cellfilter = '(isequal($area, ''CA1'') && ($meanrate < 7))';

% time filter
timefilter = { {'getlinstate', '(($traj ~= -1) & (abs($velocity) >= 0))', 6} };

% iterator
iterator = 'singlecellanal';

% filter creation
ca1f = createfilter('animal',animals,'epochs',epochfilter,'cells',ca1cellfilter,'excludetime', timefilter, 'iterator', iterator);

% set analysis function
ca1f = setfilterfunction(ca1f, 'calctotalmeanrate', {'spikes'});

% run analysis
ca1f = runfilter(ca1f);
