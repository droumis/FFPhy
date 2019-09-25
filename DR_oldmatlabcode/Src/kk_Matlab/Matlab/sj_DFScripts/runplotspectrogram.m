% Animal selection
animals = {'Eight'};

% Epoch selection
epochfilter{1} = ['isequal($description, ''TrackA'') & ($exposure == 1)'];
epochfilter{2} = ['isequal($description, ''TrackA'') & ($exposure == 7)'];

% Tetrode selection
tetrodefilter = '(isequal($area, ''MEC''))';

% Time selection
timefilter = {{'getlinstate', '(($traj ~= -1) & abs($velocity) >= 0)', 6}};

% Iterator selection
iterator = 'eeganal';


% Create and Run Filter
f = createfilter('animal',animals,'epochs',epochfilter,'excludetime',timefilter,'eegtetrodes', tetrodefilter, 'iterator', iterator);

f = setfilterfunction(f, 'calcspectrogram', {'eeg'},'appendindex',1,'fpass',[0 150]);

f = runfilter(f);