
%Animal selection
%-----------------------------------------------------
animals = {'Dudley','Miles','Conley','Bond','Frank','Nine','Ten'};
%animals = {'Dudley','Miles','Conley','Bond','Frank','Ten'};
%animals = {'Dudley'};
%animals = {'Dudley','Miles','Conley','Bond','Frank','Nine'};
%animals = {'Ten'};
%animals = {'Dudley','Miles','Conley', 'Ten'};
%animals = {'Frank'};
%animals = {'Bond'};
%-----------------------------------------------------


%Filter creation
%--------------------------------------------------------
epochfilter = [];
epochfilter{1} = ['($exposure == 1)'];
%epochfilter{2} = ['($exposure == 6)'];

ca1tetfilter = '(isequal($area, ''CA1''))';
ca3tetfilter = '(isequal($area, ''CA3''))';

ca1cellfilter = '(isequal($area, ''CA1'') && ($meanrate < 7) && ($numspikes > 200))';
ca3cellfilter = '(isequal($area, ''CA3'') && ($meanrate < 7) && ($numspikes > 200))';


timefilter = {{'getriptimes', '($nripples == 0)', [], 'cellfilter', '(isequal($area, ''CA1''))'},{'getlinstate', '(($traj ~= -1) & (abs($velocity) >= 5))', 6}};

eegfilter = {'geteegtet', 'theta', 'maxvar', 1, 'tetfilter', ca1tetfilter};
%eegfilter = {'geteegtet', 'theta', 'sametet', 1};

iterator = 'singlecelleeganal';

f = createfilter('animal',animals,'epochs',epochfilter, 'excludetimefilter', timefilter, 'cells', ca3cellfilter,'eegtetrodes', eegfilter, 'iterator', iterator);


f = setfilterfunction(f, 'plotthetamod', {'spikes', 'theta'});
f = runfilter(f);

d = numericgroupcombine(f);
