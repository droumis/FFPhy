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
%animals = {'Frank'};
%-----------------------------------------------------


%Filter creation for training data
%--------------------------------------------------------
epochfilter = [];
epochfilter{1} = ['($dailyexposure == 1) & isequal($environment, ''TrackB'')'];  
epochfilter{2} = ['($dailyexposure == 1) & isequal($environment, ''TrackA'') & ($exposureday > 3)']; 

%cellfilter = '(isequal($area, ''CA1'') && ($meanrate < 7))';
cellfilter = '(($meanrate < 7))';
timefilter = { {'getlinstate', '(($traj ~= -1) & (abs($velocity) >= 3))', 6},{'getriptimes','($nripples == 0)', [], 'cellfilter', '((isequal($area, ''CA3''))|(isequal($area, ''CA1'')))'} };
%timefilter = { {'getriptimes','($nripples == 0)', [], 'cellfilter', '(isequal($area, ''CA1''))'} };
trainingfilter = createfilter('animal',animals,'epochs',epochfilter,'cells',cellfilter,'excludetime', timefilter);
%-----------------------------------------------------------


%create training data by calulating the linearized rates of all cells 
%--------------------------------------------
iterator = 'multicellanal';
trainingfilter = setfilteriterator(trainingfilter,iterator);
trainingfilter = setfilterfunction(trainingfilter, 'calcpopulationlinfields', {'spikes','linpos'},2,3);
trainingfilter = runfilter(trainingfilter);
%-------------------------------------------------

track = 1;
numoverlap1 = [];
numoverlap2 = [];
for a = 1:length(trainingfilter)
    for e = 1:length(trainingfilter(a).output{track})
	t1 = trainingfilter(a).output{1}(e).index;   
	t2 = trainingfilter(a).output{2}(e).index;   
	corresp1 = rowfind(t1(:,[1 3 4]), t2(:,[1 3 4]));
	corresp2 = rowfind(t2(:,[1 3 4]), t1(:,[1 3 4]));
	numoverlap1 = [numoverlap1; [length(corresp1) length(find(corresp1))]];
	numoverlap2 = [numoverlap2; [length(corresp2) length(find(corresp2))]];
    end
end

[154 374]
[154 313]
