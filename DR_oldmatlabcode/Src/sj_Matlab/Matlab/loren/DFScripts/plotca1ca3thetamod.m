
%Animal selection
%-----------------------------------------------------
animals = {'Frank'};
%-----------------------------------------------------


%Filter creation
%--------------------------------------------------------
% look at the first exposure
epochfilter{1} = ['($exposure == 1)'];

% create two separate tetrode filters for CA1 and CA3
ca1tetfilter = '(isequal($area, ''CA1''))';
ca3tetfilter = '(isequal($area, ''CA3''))';

% create two separate cell filters for CA1 and CA3
ca1cellfilter = '(isequal($area, ''CA1'') && ($meanrate < 7) && ($numspikes > 200))';
ca3cellfilter = '(isequal($area, ''CA3'') && ($meanrate < 7) && ($numspikes > 200))';


% exclude ripples and low velocity periods
timefilter = {{'getriptimes', '($nripples == 0)', [], 'cellfilter', '(isequal($area, ''CA1''))'},{'getlinstate', '(($traj ~= -1) & (abs($velocity) >= 5))', 6}};

% create two eeg filters to get the maximum variance CA1 or CA3 theta tetrode
ca1eegfilter = {'geteegtet', 'theta', 'maxvar', 1, 'tetfilter', ca1tetfilter};
ca3eegfilter = {'geteegtet', 'theta', 'maxvar', 1, 'tetfilter', ca3tetfilter};

iterator = 'singlecelleeganal';

% create the ca1 and ca3 filtersj
ca1f = createfilter('animal',animals,'epochs',epochfilter, 'excludetimefilter', timefilter, 'cells', ca1cellfilter,'eegtetrodes', ca1eegfilter, 'iterator', iterator);

ca3f = createfilter('animal',animals,'epochs',epochfilter, 'excludetimefilter', timefilter, 'cells', ca3cellfilter,'eegtetrodes', ca3eegfilter, 'iterator', iterator);


% plot theta modulation for each cell from CA1 or CA3
ca1f = setfilterfunction(ca1f, 'plotthetamod', {'spikes', 'theta'});
ca1f = runfilter(ca1f);

ca3f = setfilterfunction(ca3f, 'plotthetamod', {'spikes', 'theta'});
ca3f = runfilter(ca3f);
