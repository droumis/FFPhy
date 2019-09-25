
animals = {'Putin'};

epochfilter = [];
epochfilter{1} = {'task','isequal($task,''markov'')'};
epochfilter{2} = {'task','isequal($task,''markov'')'};

datafilter = [];
datafilter{1} = {'spikes','isequal($area, ''cg1'')','trials','$rewarded == 1'};
datafilter{2} = {'spikes','isequal($area, ''cg1'')','trials','$rewarded == 0'};

timefilter = [];
timefilter{1} =  {'pos', '($vel > 1)','trials', '$rewarded == 1'};

%timefilter{1} = {'<function> get2dstate <argname> immobilecutoff <argval> 1','($immobilitytime > 5)'};

filterfunction = {'calctotalmeanrate',{'spikes'},'appendindex',1};

f = createfilter('animal', animals, 'epochs', epochfilter, 'data', datafilter, 'excludetime', timefilter,'function',filterfunction);
f = runfilter(f);


%Change the function and run again 
f = setfilterfunction(f, 'calcspatialrate', {'spikes','pos'});
f = runfilter(f);




