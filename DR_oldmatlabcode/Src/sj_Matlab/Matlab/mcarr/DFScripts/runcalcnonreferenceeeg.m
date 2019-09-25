% Animal selection
animals = {'Corriander'};

% Epoch selection
epochfilter = [];
epochfilter{1} = ['isequal($type, ''run'')'];

% Tetrode selection
tetrodepairfilter = {'(~isequal($area, ''Reference''))', '(isequal($area, ''Reference''))'};

% Iterator selection
iterator = 'eeganal';

% Create and Run Filter
f = createfilter('animal',animals,'epochs',epochfilter,'eegtetrodepairs', tetrodepairfilter, 'iterator', iterator);

f = setfilterfunction(f, 'calcnonreferenceeeg', {'eeg'},'/data13/mcarr/Cor/','Cor');

f = runfilter(f);