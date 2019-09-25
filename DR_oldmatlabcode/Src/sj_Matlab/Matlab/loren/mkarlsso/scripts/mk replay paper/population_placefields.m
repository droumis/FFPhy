%Animal selection
%-----------------------------------------------------
animals = {'Conley','Bond','Frank'};
%-----------------------------------------------------


%Filter creation for training data
%--------------------------------------------------------
epochfilter = [];
epochfilter{1} = ['($dailyexposure == 1) & isequal($environment, ''TrackB'')'];  
epochfilter{2} = ['($dailyexposure == 1) & isequal($environment, ''TrackA'') & ($exposureday > 3)']; 

%cellfilter = '(isequal($area, ''CA1'') && ($meanrate < 7))';
cellfilter = '(($meanrate < 7))';
timefilter = {{'getlinstate', '(($traj ~= -1) & (abs($velocity) >= 3))',6}, ...
   {'getriptimes','($nripples == 0)', [], 'cellfilter', ...
   '((isequal($area, ''CA3''))|(isequal($area, ''CA1'')))'}};

timefilter1 = {{'getlinstate', '(($traj ~= -1) & (abs($velocity) >= 3))',6}, ...
   {'getriptimes','($nripples == 0)', [], 'cellfilter', ...
   '((isequal($area, ''CA3''))|(isequal($area, ''CA1'')))'}, ...
   {'getsplittimes', '($includeseg == 1)', 2, 1 }};

timefilter2 = {{'getlinstate', '(($traj ~= -1) & (abs($velocity) >= 3))',6}, ...
   {'getriptimes','($nripples == 0)', [], 'cellfilter', ...
   '((isequal($area, ''CA3''))|(isequal($area, ''CA1'')))'}, ...
   {'getsplittimes', '($includeseg == 1)', 2, 2 }};

pf = createfilter('animal',animals,'epochs',epochfilter,'cells',cellfilter,'excludetime', timefilter);
pf1 = createfilter('animal',animals,'epochs',epochfilter,'cells',cellfilter,'excludetime', timefilter1);
pf2 = createfilter('animal',animals,'epochs',epochfilter,'cells',cellfilter,'excludetime', timefilter2);
%-----------------------------------------------------------


%create training data by calulating the linearized rates of all cells 
%--------------------------------------------
iterator = 'multicellanal';
pf = setfilteriterator(pf,iterator);
pf = setfilterfunction(pf, 'calcpopulationlinfields', {'spikes','linpos'},2,3);
pf = runfilter(pf);
pf1 = setfilteriterator(pf1,iterator);
pf1 = setfilterfunction(pf1, 'calcpopulationlinfields', {'spikes','linpos'},2,0);
pf1 = runfilter(pf1);
pf2 = setfilteriterator(pf2,iterator);
pf2 = setfilterfunction(pf2, 'calcpopulationlinfields', {'spikes','linpos'},2,0);
pf2 = runfilter(pf2);
%-------------------------------------------------


% go through each field from the training filter (pf) and calculate the overlap
% for the pf1 and pf2 fields
for track = 1:2
    o{track} = [];
%    for a = 1:length(pf)
    for a = [1 3]
	for e = 1:length(pf(a).output{track})
	    ind = pf(a).output{track}(e).index;
	    p1 = pf1(a).output{track}(e);
	    p2 = pf2(a).output{track}(e);
	    % go through the list of cells
	    for i = 1:size(ind,1)
		t1ind = rowfind(ind(i,:), p1.index);
		t2ind = rowfind(ind(i,:), p2.index);
		o{track}(end+1,:) = [overlap(p1.rates(t1ind,:),  ...
				    p2.rates(t2ind,:)) ...
		       overlap(p1.rates(t1ind,:) ./ max(p1.rates(t1ind,:)),  ...
		               p2.rates(t2ind,:) ./ max(p2.rates(t2ind,:))),  ...
			    max(pf(a).output{track}(e).rates(i,:))];
	    end
	end
    end
end



