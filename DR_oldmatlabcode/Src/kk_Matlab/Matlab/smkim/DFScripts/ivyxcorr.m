
%Animal selection
%-----------------------------------------------------
%animals = {'Dudley','Miles','Conley','Bond','Frank','Nine','Ten'};
%animals = {'Bond','Frank'};
%animals = {'Dudley','Miles','Conley','Bond','Frank','Nine'};
%animals = {'Ten'};
%animals = {'Dudley','Miles','Conley', 'Ten'};
animals = {'Frank'};
%animals = {'Bond'};
%-----------------------------------------------------


%Filter creation
%--------------------------------------------------------
epochfilter = [];
epochfilter{1} = ['($exposure > 0)'];

%epochfilter{2} = ['($exposure == 6)'];

%ca1cellfilter = '(isequal($area, ''CA1'') && ($meanrate < 7) && ($propbursts < 0.) && ($csi < .05) && ($numspikes > 100))';
%ca3cellfilter = '(isequal($area, ''CA3'') && ($meanrate < 7) && ($propbursts < 0.1) && ($csi < .05) && ($numspikes > 100))';
cellpairfilter = {'allcomb', '(isequal($area, ''CA3'') && ($meanrate < 7) && ($propbursts > .2) && ($csi > 0.1) && ($numspikes > 200))', '(isequal($area, ''CA3'') && ($meanrate < 7) && ($propbursts < 0.05) && ($csi < 0.05) && ($numspikes > 200))'};

%timefilter = {{'getriptimes', '($nripples > 1)', [], 'cellfilter', '(isequal($area, ''CA1''))'}};
timefilter = { {'getlinstate', '(($traj ~= -1) & (abs($velocity) >= 0))', 6} };


iterator = 'singlecellanal';

f = createfilter('animal',animals,'epochs',epochfilter,'cellpairs',cellpairfilter,'excludetimefilter', timefilter, 'iterator', iterator);

f = setfilterfunction(f, 'calcxcorrmeasures', {'spikes', 'linfields'}, 'edgespikes', 1);
f = runfilter(f);

