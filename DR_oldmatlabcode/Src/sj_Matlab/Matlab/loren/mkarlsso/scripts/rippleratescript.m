%Animal selection
%-----------------------------------------------------
%animals = {'Dudley','Miles','Conley','Ten','Bond','Frank','Nine'};
%animals = {'Dudley','Miles','Conley','Alex','Bond','Frank','Nine','Ten'};


animals = {'Conley','Bond','Frank'};
%animals = {'Miles','Nine', 'Ten'};

%animals = {'Ten'};
%animals = {'Dudley','Miles','Conley', 'Ten'};
%animals = {'Bond','Frank','Nine'};
%animals = {'Bond'};
%-----------------------------------------------------


%Filter creation for training data
%--------------------------------------------------------
epochfilter = [];
%epochfilter{1} = ['($dailyexposure == 1) & isequal($environment, ''TrackB'')'];  
epochfilter{1} = ['($dailyexposure == 1) & isequal($environment, ''TrackA'') & ($exposureday > 3)'];  
%epochfilter{1} = ['(isequal($type,''sleep'')) & isequal($runbefore.environment, ''TrackA'') & ($runbefore.exposureday > 3)'];  
%epochfilter{1} = ['(isequal($type,''sleep'')) & isequal($runbefore.environment, ''TrackB'')'];  



%epochfilter{1} = ['($dailyexposure == 1) & isequal($environment, ''TrackB'') & ($exposureday <= 4)'];  
%epochfilter{2} = ['($dailyexposure == 1) & isequal($environment, ''TrackA'') & ($exposureday > 3)'];  
%cellfilter = '(isequal($area, ''CA1'') && ($meanrate < 7))';
%cellfilter = '(($meanrate < 7) && ($meanrate > .1))';
%timefilter = { {'getlinstate', '(($traj ~= -1) & (abs($velocity) >= 3))', 6},{'getriptimes','($nripples == 0)', [], 'cellfilter', '((isequal($area, ''CA3''))|(isequal($area, ''CA1'')))'} };
%timefilter = { {'getriptimes','($nripples == 0)', [], 'cellfilter', '(isequal($area, ''CA1''))'} };
timefilter = {{'get2dstate', '($immobilitytime > 2.4)'}};

%timefilter = {};
filter = createfilter('animal',animals,'epochs',epochfilter,'excludetime', timefilter);
%-----------------------------------------------------------

%create training data by calulating the linearized rates of all cells 
%--------------------------------------------
iterator = 'singleepochanal';
filter = setfilteriterator(filter,iterator);
filter = setfilterfunction(filter, 'calcripplerate', {'pos','ripples','cellinfo'});
filter = runfilter(filter);
%-------------------------------------------------

groups = numericgroupcombine(filter,0);