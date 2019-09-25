
%Animal selection
%-----------------------------------------------------
%animals = {'Dudley','Miles','Conley','Bond','Frank','Nine','Ten'};
%animals = {'Dudley','Miles','Conley','Alex','Bond','Frank','Nine'};
%animals = {'Ten'};
%animals = {'Dudley','Miles','Conley', 'Ten'};
animals = {'Frank'};
%animals = {'Bond'};
%-----------------------------------------------------


%Filter creation
%--------------------------------------------------------
epochfilter = [];
epochfilter{1} = ['(isequal($type, ''run'')) & ($exposure == 19) & (isequal($environment, ''TrackA''))'];

cellfilter = '((isequal($area, ''CA3'') || isequal($area, ''CA1'')) && ($meanrate < 7))';

timefilter = {{'getriptimes', '($nripples == 0)', [], 'cellfilter', '(isequal($area, ''CA1''))'}, {'getlinstate', '(($traj ~= -1) & (abs($velocity) >= 5))', 6}};

iterator = 'singlecellsave';

f = createfilter('animal',animals,'epochs',epochfilter,'cells',cellfilter,'excludetimefilter', timefilter, 'iterator', iterator);
%-----------------------------------------------------------

% keep only the spikes outside of ripples
%--------------------------------------------
f = setfilterfunction(f, 'filterspikes', {'spikes'});

f = runfilter(f);
