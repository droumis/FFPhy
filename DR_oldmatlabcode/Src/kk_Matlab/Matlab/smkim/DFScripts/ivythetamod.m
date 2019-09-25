
%Animal selection
%-----------------------------------------------------
%animals = {'Bond','Frank','Ten'};
animals = {'Alex', 'Dudley','Miles','Conley','Bond','Frank','Nine'};
%animals = {'Ten'};
%animals = {'Dudley','Miles','Conley', 'Ten'};
%animals = {'Bond','Frank','Nine'};
%animals = {'Frank'};
%animals = {'Bond'};
%-----------------------------------------------------


%Filter creation
%--------------------------------------------------------
epochfilter{1} = ['isequal($environment, ''TrackA'') | isequal($environment, ''TrackB'')'];
%for i = 1:14
    %epochfilter{i} = ['isequal($environment, ''TrackA'') | isequal($environment, ''TrackB'')'];
%   epochfilter{i} = ['isequal($environment, ''TrackA'')'];
%nd

ca1cellfilter = '(isequal($area, ''CA1'') && ($meanrate < 7) && ($propbursts < 0.04) && ($csi < .05) && ($numspikes > 100))';
ca3cellfilter = '(isequal($area, ''CA3'') && ($meanrate < 7) && ($propbursts < 0.04) && ($csi < .05) && ($numspikes > 100))';

timefilter = { {'getlinstate', '(($traj ~= -1) & (abs($velocity) >= 0))', 6} };
%timefilter = {{'getriptimes', '($nripples > 1)', [], 'cellfilter', '(isequal($area, ''CA1''))'}};

%timefilter = {};

iterator = 'singlecellanal';
%iterator = 'multicellanal';

ca1f = createfilter('animal',animals,'epochs',epochfilter,'cells',ca1cellfilter,'excludetime', timefilter, 'iterator', iterator);
ca3f = createfilter('animal',animals,'epochs',epochfilter,'cells',ca3cellfilter,'excludetime', timefilter, 'iterator', iterator);
%-----------------------------------------------------------




%ca1f = setfilterfunction(ca1f, 'filtercalclinfields', {'spikes', 'linpos'});
ca1f = setfilterfunction(ca1f, 'plotthetamod', {'spikes'});
ca3f = setfilterfunction(ca3f, 'plotthetamod', {'spikes'});
%ca3f = setfilterfunction(ca1f, 'filtercalclinfields', {'spikes', 'linpos'});

'ca1'
ca1f = runfilter(ca1f);
'ca3'
ca3f = runfilter(ca3f);
 
