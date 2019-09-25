%Animal selection
%-----------------------------------------------------
%animals = {'Dudley','Miles','Conley','Ten','Bond','Frank','Nine'};
%animals = {'Dudley','Miles','Conley','Alex','Bond','Frank','Nine','Ten'};


animals = {'Conley','Bond','Frank','Miles','Nine', 'Ten'};
%animals = {'Miles','Nine', 'Ten'};

%animals = {'Ten'};
%animals = {'Dudley','Miles','Conley', 'Ten'};
%animals = {'Bond','Frank','Nine'};
%animals = {'Bond'};
%-----------------------------------------------------


%Filter creation for training data
%--------------------------------------------------------
epochfilter = [];
epochfilter{1} = ['($dailyexposure == 1) & isequal($environment, ''TrackB'')'];  
epochfilter{2} = ['($dailyexposure == 1) & isequal($environment, ''TrackA'') & ($exposureday > 3)'];  
%cellfilter = '(isequal($area, ''CA1'') && ($meanrate < 7))';
cellfilter = '(($meanrate < 7) && ($meanrate > .1))';
%timefilter = { {'getlinstate', '(($traj ~= -1) & (abs($velocity) >= 3))', 6},{'getriptimes','($nripples == 0)', [], 'cellfilter', '((isequal($area, ''CA3''))|(isequal($area, ''CA1'')))'} };
%timefilter = { {'getriptimes','($nripples == 0)', [], 'cellfilter', '((isequal($area, ''CA3''))|(isequal($area, ''CA1'')))'} };
timefilter = { {'getlinstate', '(abs($velocity) < 3)', 6},{'getriptimes','($nripples >= 2)', [], 'cellfilter', '((isequal($area, ''CA3''))|(isequal($area, ''CA1'')))'} };
filter = createfilter('animal',animals,'epochs',epochfilter,'cells',cellfilter,'excludetime', timefilter);
%-----------------------------------------------------------

%Get the binned spike counts for all cells
%-----------------------------------------------------------
window = .5;
cellcountthresh = 1;
iterator = 'multicellanal';
filter = setfilteriterator(filter,iterator);
filter = setfilterfunction(filter, 'getpopulationevents', {'spikes','linpos'},window,cellcountthresh);
filter = runfilter(filter);
%----------------------------------------------------------

%list = [1 1;1 2;2 1;2 2;2 3;3 2;3 3];
list = [1 1;2 1;3 1;3 2;3 3];
propactiveNov = [];
propactiveFam = [];
for i = 1:size(list,1)
    %out = [out; calcepochreplaystats(list(i,:), 1, trainingfilter, decodefilter)];
    out = calcepochactivation(list(i,:), filter);
    propactiveNov = [propactiveNov; out{1}];
    propactiveFam = [propactiveFam; out{2}];
    
end
%---------------------------------------------------------