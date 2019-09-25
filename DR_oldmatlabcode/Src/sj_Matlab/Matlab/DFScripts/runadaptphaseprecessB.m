
%Animal selection
%-----------------------------------------------------
animals = {'Dudley','Miles','Conley','Bond','Frank','Nine','Ten'};
%animals = {'Dudley','Miles','Conley','Alex','Bond','Frank','Nine'};
%animals = {'Ten'};
%animals = {'Dudley','Miles','Conley', 'Ten'};
%animals = {'Bond','Frank','Nine'};
%animals = {'Frank'};
%animals = {'Bond'};
%-----------------------------------------------------


%Filter creation
%--------------------------------------------------------
%epochfilter{1} = ['isequal($environment, ''TrackA'') | isequal($environment, ''TrackB'')'];
%epochfilter{1} = ['isequal($environment, ''TrackA'')'];
%epochfilter{1} = ['($experimentday == 7)'];
for i = 1:14
    epochfilter{i} = ['isequal($environment, ''TrackB'') & ($exposure == ', num2str(i), ')'];
end

ca1cellfilter = '(isequal($area, ''CA1'') && ($meanrate < 7) && ($numspikes > 100))';
ca3cellfilter = '(isequal($area, ''CA3'') && ($meanrate < 7) && ($numspikes > 100))';

%timefilter = { {'getlinstate', '(($traj ~= -1) & (abs($velocity) >= 0))', 3} };
timefilter = { {'getlinstate', '(($traj ~= -1) & (abs($velocity) >= 3))', 6},{'getriptimes','($nripples == 0)', [], 'minthresh', 3, 'cellfilter', '((isequal($area, ''CA3''))|(isequal($area, ''CA1'')))'} };

ca3tetfilter = '(isequal($area, ''CA3''))';
ca1tetfilter = '(isequal($area, ''CA1''))';

ca3eegfilter = {'geteegtet', 'theta', 'maxvar', 1, 'tetfilter', ca3tetfilter};
ca1eegfilter = {'geteegtet', 'theta', 'maxvar', 1, 'tetfilter', ca1tetfilter};

iterator = 'singlecelleeganal';
%iterator = 'multicellanal';

ca3f = createfilter('animal',animals,'epochs',epochfilter,'cells',ca3cellfilter,'excludetime', timefilter, 'eegtetrodes', ca3eegfilter, 'iterator', iterator);
%-----------------------------------------------------------
ca1f = createfilter('animal',animals,'epochs',epochfilter,'cells',ca1cellfilter,'excludetime', timefilter, 'eegtetrodes', ca1eegfilter, 'iterator', iterator);
%-----------------------------------------------------------




ca3f = setfilterfunction(ca3f, 'calcphaseprecessadapt', {'spikes', 'linpos', 'theta'});
ca1f = setfilterfunction(ca3f, 'calcphaseprecessadapt', {'spikes', 'linpos', 'theta'});

ca3f = runfilter(ca3f);

save /data/mkarlsso/Adaptive/ca3trackB.mat  ca3f;
clear ca3f;

ca1f = runfilter(ca1f);

save /data/mkarlsso/Adaptive/ca1trackB.mat  ca1f;
clear ca1f;
