% Animal Selection
animals = {'Fear'};

% Epoch Filter
epochfilter{1} = ['($exposure == 1)'];

% Time Filter
timefilter = { {'getlinstate', '(($traj ~= -1) & abs($velocity) <= 2)', 6} };

% Tetrode Filter
tetrodefilter = '(isequal($area, ''CA1''))';

% Iterator
iterator = 'multitetrodeanal';

% Filter Creation
f = createfilter('animal', animals, 'epochs', epochfilter, 'eegtetrodes', tetrodefilter, 'excludetime', timefilter, 'iterator', iterator);

% Set Analysis Function
f = setfilterfunction(f, 'getripratehist', {'ripples'},'binsize',30, 'numtetrodes', 1,'minthresh',7);

% Run Analysis
f = runfilter(f);