%Animal selection
%-----------------------------------------------------
%animals = {'Miles','Conley','Bond','Frank','Nine','Ten'};

%animals = {'Miles','Nine','Ten'};
animals = {'Conley','Bond','Frank'};


%animals = {'Bond'};
%animals = {'Miles','Conley','Ten','Bond','Frank','Nine'};
%animals = {'Dudley','Miles','Conley','Alex','Bond','Frank','Nine'};
%animals = {'Ten'};
%animals = {'Dudley','Miles','Conley', 'Ten'};
%animals = {'Bond','Frank','Nine'};
%animals = {'Bond'};
%-----------------------------------------------------
epochfilters1 = [];
epochfilters2 = [];

epochfilters1{1} = ['($dailyexposure == 1) & isequal($environment, ''TrackB'') & ($exposureday <= 12)'];
% epochfilters1{2} = ['($dailyexposure == 1) & isequal($environment, ''TrackA'')'];
epochfilters2{1} = ['($dailyexposure == 1) & isequal($environment, ''TrackA'')'];
% epochfilters2{2} = ['($dailyexposure == 1) & isequal($environment, ''TrackB'') & ($exposureday <= 4)'];

%epochfilters1{1} = ['($dailyexposure == 1) & isequal($environment, ''TrackB'') & ($exposureday <= 12)'];
%epochfilters2{1} = ['isequal($type, ''sleep'') & ($runbefore.dailyexposure == 2) & isequal($runbefore.environment, ''TrackB'')'];
%epochfilters2{1} = ['isequal($type, ''sleep'') & ($runafter.dailyexposure == 1) & isequal($runafter.environment, ''TrackB'')'];


ca1cmp = [];
ca3cmp = [];

for i = 1:length(epochfilters1)
    %Filter creation
    %--------------------------------------------------------
    epochfilter1 = [];
    epochfilter2 = [];
    epochfilter1{1} = epochfilters1{i};
    epochfilter2{1} = epochfilters2{i};

    CA1cellfilter1 = '((isequal($area, ''CA1'')) && ($meanrate < 7) && ($meanrate >= 0))';
    CA1cellfilter2 = '((isequal($area, ''CA1'')) && ($meanrate < 7) && ($meanrate < .1))';
    %CA1cellfilter2 = '((isequal($area, ''CA1'')) && ($meanrate < 7))';
    

    CA3cellfilter1 = '((isequal($area, ''CA3'')) && ($meanrate < 7) && ($meanrate > .1))';
    %CA3cellfilter2 = '((isequal($area, ''CA3'')) && ($meanrate < 7) && ($meanrate < .1))';
    CA3cellfilter2 = '((isequal($area, ''CA3'')) && ($meanrate < 7))';

    %timefilter1 = { {'getlinstate', '(abs($velocity) >= 3)',6},{'getriptimes', '($nripples == 0)', [], 'cellfilter', '(isequal($area, ''CA1''))|(isequal($area, ''CA3''))','minstd',3} };
    timefilter1 = {{'getriptimes', '($nripples > 2)', [], 'cellfilter', '(isequal($area, ''CA1''))','minstd',3}};
    %timefilter2 = {{'getriptimes', '($nripples > 0)', [], 'cellfilter', '(isequal($area, ''CA1''))|(isequal($area, ''CA3''))','minstd',3}};
    timefilter2 = {{'getriptimes', '($nripples > 2)', [], 'cellfilter', '(isequal($area, ''CA1''))','minstd',3}};
    %timefilter1 = {};
    %timefilter2 = {};

    CA1f1 = createfilter('animal',animals,'epochs',epochfilter1,'cells',CA1cellfilter1,'excludetime', timefilter1);
    CA1f2 = createfilter('animal',animals,'epochs',epochfilter2,'cells',CA1cellfilter2,'excludetime', timefilter2);

    %CA3f1 = createfilter('animal',animals,'epochs',epochfilter1,'cells',CA3cellfilter1,'excludetime', timefilter1);
    %CA3f2 = createfilter('animal',animals,'epochs',epochfilter2,'cells',CA3cellfilter2,'excludetime', timefilter2);

    %-----------------------------------------------------------

    %Pairwise comparison
    %-----------------------------------------------------------
    iterator = 'singlecellanal';

    CA1f1 = setfilteriterator(CA1f1,iterator);
    CA1f2 = setfilteriterator(CA1f2,iterator);
    %CA3f1 = setfilteriterator(CA3f1,iterator);
    %CA3f2 = setfilteriterator(CA3f2,iterator);

    CA1f1 = setfilterfunction(CA1f1, 'calctotalmeanrate', {'spikes'},'appendindex',1);
    CA1f2 = setfilterfunction(CA1f2, 'calctotalmeanrate', {'spikes'},'appendindex',1);
    %CA3f1 = setfilterfunction(CA3f1, 'calctotalmeanrate', {'spikes'},'appendindex',1);
    %CA3f2 = setfilterfunction(CA3f2, 'calctotalmeanrate', {'spikes'},'appendindex',1);


    %CA1f1 = setfilterfunction(CA1f1, 'calcrippleactivationprob', {'spikes','ripples','cellinfo'},'appendindex',1,'ratio',1);
    %CA1f2 = setfilterfunction(CA1f2, 'calcrippleactivationprob', {'spikes','ripples','cellinfo'},'appendindex',1,'ratio',1);

    CA1f1 = runfilter(CA1f1);
    CA1f2 = runfilter(CA1f2);
    %CA3f1 = runfilter(CA3f1);
    %CA3f2 = runfilter(CA3f2);

    CA1f1groups = numericgroupcombine(CA1f1,1);
    CA1f2groups = numericgroupcombine(CA1f2,1);
    %CA3f1groups = numericgroupcombine(CA3f1,1);
    %CA3f2groups = numericgroupcombine(CA3f2,1);

    indexcolumns = [1 2 4 5];
    ca1cmp{i} = indexmatch(CA1f1groups{1},CA1f2groups{1},indexcolumns);
    %ca3cmp{i} = indexmatch(CA3f1groups{1},CA3f2groups{1},indexcolumns);

    %--------------------------------------------------------------
end