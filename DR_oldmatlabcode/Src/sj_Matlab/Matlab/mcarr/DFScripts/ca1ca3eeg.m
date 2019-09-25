
%Animal selection
%-----------------------------------------------------
%animals = {'Dudley','Miles','Conley','Bond','Frank','Nine','Ten'};
%animals = {'Dudley','Miles','Conley','Bond','Frank','Ten'};
%animals = {'Dudley'};
%animals = {'Dudley','Miles','Conley','Bond','Frank','Nine'};
%animals = {'Ten'};
%animals = {'Dudley','Miles','Conley', 'Ten'};
animals = {'Frank'};
%animals = {'Bond'};
%-----------------------------------------------------


%Filter creation
%--------------------------------------------------------
epochfilter = [];
epochfilter{1} = ['($exposure == 1)'];
%epochfilter{2} = ['($exposure == 6)'];

tetrodepairfilter = {'(isequal($area, ''CA3''))', '(isequal($area, ''CA3''))'};

timefilter = {{'getriptimes', '($nripples == 0)', [], 'cellfilter', '(isequal($area, ''CA1''))'},{'getlinstate', '(($traj ~= -1) & (abs($velocity) >= 5))', 6}};


iterator = 'eeganal';

f = createfilter('animal',animals,'epochs',epochfilter,'tetrodepairs',tetrodepairfilter,'excludetimefilter', timefilter, 'iterator', iterator);

f = setfilterfunction(f, 'calcphasediff', {'theta'});
f = runfilter(f);
