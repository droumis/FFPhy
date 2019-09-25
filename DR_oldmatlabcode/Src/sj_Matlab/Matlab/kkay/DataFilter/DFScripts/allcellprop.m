
%Animal selection
%-----------------------------------------------------
%animals = {'Bond','Frank','Ten'};
%animals = {'Dudley','Miles','Conley','Bond','Frank','Nine','Ten'};
%animals = {'Miles','Conley','Bond','Frank','Nine','Ten'};
%animals = {'Ten'};
%animals = {'Dudley','Miles','Conley', 'Ten'};
%animals = {'Bond','Frank','Nine'};
%animals = {'Frank'};
animals = {'Bond'};
%-----------------------------------------------------


%Filter creation
%--------------------------------------------------------
epochfilter{1} = ['isequal($environment, ''TrackB'') && ($dailyexposure == 1)'];
%for i = 1:14
    %epochfilter{i} = ['isequal($environment, ''TrackA'') | isequal($environment, ''TrackB'')'];
%   epochfilter{i} = ['isequal($environment, ''TrackA'')'];
%nd

ca1cellfilter = '(isequal($area, ''CA1'') && ($numspikes > 100))';
ca3cellfilter = '(isequal($area, ''CA3'') && ($numspikes > 100))';

% theta periods
timefilter = { {'getlinstate', '(($traj ~= -1) & (abs($velocity) >= 4))', 6} };
% ripple periods --- FIX this w/ tetlist
%timefilter = {{'getriptimes', '($nripples > 1)', [], 'cellfilter', '(isequal($area, ''CA1''))'}};

%timefilter = {};

iterator = 'singlecellanal';
%iterator = 'multicellanal';

ca1f = createfilter('animal',animals,'epochs',epochfilter,'cells',ca1cellfilter,'iterator', iterator, 'excludetime', timefilter);
%ca3f = createfilter('animal',animals,'epochs',epochfilter,'cells',ca3cellfilter,'iterator', iterator);
%-----------------------------------------------------------

keyboard



%ca1f = setfilterfunction(ca1f, 'cellprop', {'spikes'});
%ca3f = setfilterfunction(ca3f, 'cellprop', {'spikes'});
ca3f = setfilterfunction(ca1f, 'filtercalclinfields', {'spikes', 'linpos'});
ca1f = setfilterfunction(ca1f, 'filtercalclinfields', {'spikes', 'linpos'});

'ca1'
ca1f = runfilter(ca1f);

ca1 = numericgroupcombine(ca1f);
'ca3'
ca3f = runfilter(ca3f);
ca3 = numericgroupcombine(ca3f);
c3 = ca3{1};
 
