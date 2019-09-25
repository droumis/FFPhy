%Animal selection
%-----------------------------------------------------
animals = {'Conley','Bond','Frank'};
%-----------------------------------------------------


%Filter for getting cluster quality
%--------------------------------------------------------
epochfilter = [];
epochfilter{1} = ['($sleepnum == 1)'];  
epochfilter{2} = ['($dailyexposure == 1) & isequal($environment, ''TrackB'')']; 
epochfilter{3} = ['($sleepnum == 2)'];  
epochfilter{4} = ['($dailyexposure == 2) & isequal($environment, ''TrackB'')']; 
epochfilter{5} = ['($sleepnum == 3)'];  
epochfilter{6} = ['($dailyexposure == 1) & isequal($environment, ''TrackA'') & ($exposureday > 3)']; 
epochfilter{7} = ['($sleepnum == 4)'];  

cellfilter = '(($meanrate < 7))';
timefilter = {};

iterator = 'singlecellanal';

qf = createfilter('animal',animals,'epochs',epochfilter,'cells',...
		 cellfilter,'excludetime', timefilter, 'iterator', iterator);
qf = setfilterfunction(qf, 'getclusterqual', {'clustqual'}, 'appendindex', 1);
qf = runfilter(qf);



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

cq = [];
% Go through each set of cells from training filter and find each of the epochs
% where that cell has a set of cluster quality measures.
for a = 1:length(trainingfilter)
    for env = 1:length(trainingfilter(a).output)
	for epoch = 1:length(trainingfilter(a).output{env})
	    ind = trainingfilter(a).output{env}(epoch).index;
	    % find these cells in each epoch of qf
	    nqfe = length(qf(a).output)
	    q = ones(size(ind,1), nqfe) * -1;
	    for i = 1:length(qf(a).output)
		tmpi = rowfind(ind(:, [1 3 4]), qf(a).output{i}(:,[1 3 4]));
		q(i,:) = qf(a).output{i}(tmpi,6);
	    end
	end
	%cq = [cq q];
    end
end


	        
