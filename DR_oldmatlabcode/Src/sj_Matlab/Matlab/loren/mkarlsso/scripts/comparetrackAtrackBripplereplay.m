%Animal selection
%-----------------------------------------------------
animals = {'Dudley''Miles','Conley','Bond','Frank','Nine','Ten'};

%animals = {'Dudley','Miles','Conley','Alex','Bond','Frank','Nine','Ten'};
%animals = {'Alex'};
%animals = {'Conley', 'Ten'};
%animals = {'Bond','Frank','Nine'};
%animals = {'Bond'};
%-----------------------------------------------------


%Filter creation - get across-day groups
%--------------------------------------------------------
epochfilter = [];
for i = 1:14
    epochfilter{i} = ['($exposureday == ',num2str(i),') & ($dailyexposure == 1)']; %we want the index to the first daily run on each track
   
end

ca1cellfilter = '(isequal($area, ''CA1'') && ($meanrate < 7))';
ca3cellfilter = '(isequal($area, ''CA3'') && ($meanrate < 7))';

%timefilter = { {'getlinstate', '(($traj ~= -1) & (abs($velocity) >= 3))', 6} };
timefilter = {};

ca1f = createfilter('animal',animals,'epochs',epochfilter,'cells',ca1cellfilter,'excludetime', timefilter);
ca3f = createfilter('animal',animals,'epochs',epochfilter,'cells',ca3cellfilter,'excludetime', timefilter);
%-----------------------------------------------------------


%run function to get day group cell indices
%--------------------------------------------
iterator = 'singlecellanal';

ca1f = setfilteriterator(ca1f,iterator);
ca3f = setfilteriterator(ca3f,iterator);

%Get the indices the the cells (first run of the day)
ca1f = setfilterfunction(ca1f, 'getindex', {});
ca3f = setfilterfunction(ca3f, 'getindex', {});
%ca1f = setfilterfunction(ca1f, 'funcSwitchBox', {'spikes','linpos'},'calcpeakrate','calctotalmeanrate');
%ca3f = setfilterfunction(ca3f, 'funcSwitchBox', {'spikes','linpos'},'calcpeakrate','calctotalmeanrate');

ca1f = runfilter(ca1f);
ca3f = runfilter(ca3f);

ca1groups = numericgroupcombine(ca1f,1);
ca3groups = numericgroupcombine(ca3f,1);

%-------------------------------------------------


%Filter creation (for pairs)
%--------------------------------------------------------
epochfilter = [];

epochfilter{1} = '( ($dailyexposure == 1) & isequal($environment, ''TrackB'') )';
epochfilter{2} = '( ($dailyexposure == 2) & isequal($environment, ''TrackB'') )';
epochfilter{3} = '( ($dailyexposure == 1) & isequal($environment, ''TrackA'') )'; 

ca1cellfilter = '(isequal($area, ''CA1'') && ($meanrate < 7))';
ca3cellfilter = '(isequal($area, ''CA3'') && ($meanrate < 7))';

timefilter = { {'getlinstate', '(($traj ~= -1) & (abs($velocity) >= 3))', 6} };
%timefilter = {};

ca1f = createfilter('animal',animals,'epochs',epochfilter,'cells',ca1cellfilter,'excludetime', timefilter);
ca3f = createfilter('animal',animals,'epochs',epochfilter,'cells',ca3cellfilter,'excludetime', timefilter);
%-----------------------------------------------------------

%Pairwise comparison
%-----------------------------------------------------------
iterator = 'singlecellanal';

ca1f = setfilteriterator(ca1f,iterator);
ca3f = setfilteriterator(ca3f,iterator);

ca1f = setfilterfunction(ca1f, 'filtercalclinfields', {'spikes','linpos'});
ca3f = setfilterfunction(ca3f, 'filtercalclinfields', {'spikes','linpos'});

ca1f = runfilter(ca1f);
ca3f = runfilter(ca3f);

%ca1pairs = getSameCellLinearFields(ca1f);
%ca3pairs = getSameCellLinearFields(ca3f);

ca1pairs = getSameCellLinearFields_threeepochs(ca1f);
ca3pairs = getSameCellLinearFields_threeepochs(ca3f);

%get a list of all the first run index for each pair
ca1pairindex = [];
ca3pairindex = [];
for i = 1:length(ca1pairs)
    ca1pairindex = [ca1pairindex; ca1pairs(i).runindex(1,:)];
end
for i = 1:length(ca3pairs)
    ca3pairindex = [ca3pairindex; ca3pairs(i).runindex(1,:)];
end
%------------------------------------------------------------------

%compare place field shapes
%-------------------------------------------------------------------
ca1daycmp = [];
ca3daycmp = [];
%comparefunction = 'ratevschange'; %calulate rate vs change across bins 
comparefunction = 'overlap_threeepochs'; %calulate rate vs change for the peak only
%comparefunction = 'meanratevschange'; %calulate rate vs change for mean rate

%For each data group, pick out the pairs that match the indices in the
%group, and run COMPAREFUNCTION on the pairs
% CA1
for i = 1:length(ca1groups)
    groupcellindex = ca1groups{i}(:,[1 2 3 4 5]);
    pairindex = rowfind(groupcellindex,ca1pairindex);
    pairindex = pairindex(find(pairindex > 0));
    tmppairdata = ca1pairs(pairindex);
    ca1daycmp{i} = feval(comparefunction,tmppairdata);
end

% CA3
for i = 1:length(ca3groups)
    groupcellindex = ca3groups{i}(:,[1 2 3 4 5]);
    pairindex = rowfind(groupcellindex,ca3pairindex);
    pairindex = pairindex(find(pairindex > 0));
    tmppairdata = ca3pairs(pairindex);
    ca3daycmp{i} = feval(comparefunction,tmppairdata);
end
%--------------------------------------------------


ca1combineddata = [];
ca3combineddata = [];
for n = 1:14 %which days to combine
    %tmpdata = [ca1daycmp{n}.overlapBB ca1daycmp{n}.overlapBA];
    tmpdata = [ca1daycmp{n}.ratesB1 ca1daycmp{n}.ratesB2 ca1daycmp{n}.ratesA];
    ca1combineddata = [ca1combineddata; tmpdata];
    %tmpdata = [ca3daycmp{n}.overlapBB ca3daycmp{n}.overlapBA];
    tmpdata = [ca3daycmp{n}.ratesB1 ca3daycmp{n}.ratesB2 ca3daycmp{n}.ratesA];
    ca3combineddata = [ca3combineddata; tmpdata];
end
ratediff = abs(ca3combineddata(:,2)-ca3combineddata(:,3)) - abs(ca3combineddata(:,1)-ca3combineddata(:,2));

denindex = find(ca3combineddata(:,1)>0);
normoverlap = ca3combineddata(denindex,2)./ca3combineddata(denindex,1);











