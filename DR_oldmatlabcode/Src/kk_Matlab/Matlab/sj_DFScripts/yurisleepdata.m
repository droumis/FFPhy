
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
epochfilter{1} = ['(isequal($type, ''sleep'')) & ($experimentday == 13)'];

cellfilter = '((isequal($area, ''CA3'') || isequal($area, ''CA1'')) && ($meanrate < 7))';

timefilter = {{'getriptimes', '($nripples > 1)', [], 'cellfilter', '(isequal($area, ''CA1''))'}};

iterator = 'singlecellsave';

f = createfilter('animal',animals,'epochs',epochfilter,'cells',cellfilter,'excludetimefilter', timefilter, 'iterator', iterator);
%-----------------------------------------------------------

% keep only the spikes within ripples
%--------------------------------------------
f = setfilterfunction(f, 'filterspikes', {'spikes'});

f = runfilter(f);
