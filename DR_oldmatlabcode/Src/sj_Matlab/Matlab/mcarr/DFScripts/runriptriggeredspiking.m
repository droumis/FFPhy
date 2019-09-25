% Animal Selection
animals = {'Eight'};

% Epoch Filter
epochfilter = '(isequal($description, ''TrackA''))';
%epochfilter = ['($exposureday >= 4)'];
%epochfilter = ['($exposureday == 1)'];

% Time Filter
timefilter = { {'getlinstate', '(($traj ~= -1) & abs($velocity) >=0)', 6} };

% Cell Pair Filter
cellpairfilter = {'difftet','(isequal($area, ''CA1''))', '(isequal($area, ''MEC''))'};

% Iterator
iterator = 'multicellanal';

% Filter Creation
f = createfilter('animal', animals, 'epochs', epochfilter, 'cellpairs', cellpairfilter, 'excludetime', timefilter, 'iterator', iterator);

% Set Analysis Function
f = setfilterfunction(f, 'getriptriggeredspiking', {'ripples','spikes'},'minthresh',2,'window',[0.5 0.5]);

% Run Analysis
f = runfilter(f);