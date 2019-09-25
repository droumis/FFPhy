
%Animal selection
%-----------------------------------------------------
%animals = {'Dudley','Miles','Conley','Bond','Frank','Nine','Ten'};
%animals = {'Dudley','Miles','Conley','Bond','Frank','Nine'};
animals = {'Miles'};
%animals = {'Dudley','Miles','Conley', 'Ten'};
%animals = {'Frank'};
%animals = {'Bond'};
%-----------------------------------------------------


%Filter creation
%--------------------------------------------------------
epochfilter = [];
epochfilter{1} = ['($exposure == 1)'];
%epochfilter{2} = ['($exposure == 6)'];

cellpairfilter = {'allcomb', '(isequal($area, ''CA1'') && ($meanrate < 7))', '(isequal($area, ''CA1'') && ($meanrate < 7))'};

timefilter = {{'getriptimes', '($nripples >= 0)', [], 'cellfilter', '(isequal($area, ''CA1''))'}};


iterator = 'singlecellanal';

f = createfilter('animal',animals,'epochs',epochfilter,'cellpairs',cellpairfilter,'excludetimefilter', timefilter, 'iterator', iterator);

f = setfilterfunction(f, 'calcxcorrmeasures', {'spikes'}, 'edgespikes', 1);
f = runfilter(f);

