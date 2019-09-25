%Animal selection
%-----------------------------------------------------
%animals = {'Dudley','Miles','Conley','Ten','Bond','Frank','Nine'};
%animals = {'Dudley','Miles','Conley','Alex','Bond','Frank','Nine','Ten'};
%animals = {'Conley','Bond','Frank','Miles','Nine', 'Ten'};
%animals = {'Miles','Nine', 'Ten'};
animals = {'Conley','Bond','Frank'};
%animals = {'Ten'};
%animals = {'Dudley','Miles','Conley', 'Ten'};
%animals = {'Bond','Frank','Nine'};
%animals = {'Bond'};
%-----------------------------------------------------

%Filter creation
%--------------------------------------------------------

%epochfilter{1} = ['($dailyexposure == 1) & isequal($environment, ''TrackB'')'];  
%epochfilter{1} = ['($dailyexposure == 1) & isequal($environment, ''TrackA'') & ($exposureday > 3)'];  
%epochfilter{1} = ['(isequal($type,''sleep'')) & isequal($runbefore.environment, ''TrackA'') & ($runbefore.exposureday > 3)'];  
%epochfilter{1} = ['(isequal($type,''sleep'')) & isequal($runbefore.environment, ''TrackB'')'];  
%epochfilter{1} = ['(isequal($type,''sleep'')) &
epochfilter{1} = ['(isequal($type,''sleep'')) & isequal($runbefore.environment, ''TrackB'') & isequal($runafter.environment, ''TrackB'')'];  

timefilter = {};
filter = createfilter('animal',animals,'epochs',epochfilter,'excludetime', timefilter);

%----------------------------------------------------------

%create training data by calulating the linearized rates of all cells 
%--------------------------------------------
iterator = 'singleepochanal';
filter = setfilteriterator(filter,iterator);
%filter = setfilterfunction(filter, 'calcmobility', {'pos'}, 'maxtime', 2.4);
%filter = setfilterfunction(filter, 'calcmobility', {'pos'}, 'maxtime', [1:20 inf],'mintime',0:20);
filter = setfilterfunction(filter, 'calcmobility', {'pos'}, 'maxtime', [1 inf],'mintime',[-1 1]);
filter = runfilter(filter);
%-------------------------------------------------

groups = numericgroupcombine(filter,0);


