% Animal Selection
animals = {'Five'};

% Epoch Filter
epochfilter = '(isequal($description, ''TrackA''))';

% Time Filter
timefilter = { {'getlinstate', '(($traj ~= -1) & abs($velocity) >=0)', 2} };

% Tetrode Filter
tetrodefilter = '(isequal($area, ''CA1''))';

% Iterator
iterator = 'multitetrodeanal';

% Filter Creation
f = createfilter('animal', animals, 'epochs', epochfilter, 'eegtetrodes', tetrodefilter, 'excludetime', timefilter, 'iterator', iterator);

% Set Analysis Function
f = setfilterfunction(f, 'getriprate', {'ripples'});

% Run Analysis
f = runfilter(f);