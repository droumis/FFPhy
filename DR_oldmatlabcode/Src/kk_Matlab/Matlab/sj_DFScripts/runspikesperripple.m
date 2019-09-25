% Animal Selection
animals = {'Fear'};

% Epoch Filter
epochfilter{1} = ['($exposure == 1)'];

% Time Filter
timefilter = { {'getlinstate', '(($traj ~= -1) & abs($velocity) <= 2)', 6} };

% Cell Filter
cellfilter = '(isequal($area, ''CA1'') & $numspikes >= 1)';

% Iterator
iterator = 'multicellanal';

% Filter Creation
f = createfilter('animal', animals, 'epochs', epochfilter, 'cells', cellfilter, 'excludetime', timefilter, 'iterator', iterator);

% Set Analysis Function
f = setfilterfunction(f, 'getspikesperripple', {'ripples','spikes'},'minthresh',5,'binsize',60);

% Run Analysis
f = runfilter(f);