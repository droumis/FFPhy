%Animal selection
%-----------------------------------------------------
animals = {'Dudley','Miles','Conley','Bond','Frank','Nine','Ten'};
%animals = {'Dudley','Miles','Conley','Alex','Bond','Frank','Nine','Ten'};
%animals = {'Ten'};
%animals = {'Dudley','Miles','Conley', 'Ten'};
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

timefilter = { {'getlinstate', '(($traj ~= -1) & (abs($velocity) >= 3))', 6} };

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
ca1f = setfilterfunction(ca1f, 'funcSwitchBox', {'spikes','linpos'},'calcpeakrate','calctotalmeanrate');
ca3f = setfilterfunction(ca3f, 'funcSwitchBox', {'spikes','linpos'},'calcpeakrate','calctotalmeanrate');

ca1f = runfilter(ca1f);
ca3f = runfilter(ca3f);

ca1groups = numericgroupcombine(ca1f,1);
ca3groups = numericgroupcombine(ca3f,1);

%-------------------------------------------------


%Filter creation (for pairs)
%--------------------------------------------------------
epochfilter = [];

epochfilter{1} = '( (($dailyexposure == 1) & isequal($environment, ''TrackB'') & ($experimentday > 3)) | (($dailyexposure == 1) & isequal($environment, ''TrackA'') & ($experimentday <= 3)) )';
epochfilter{2} = '( (($dailyexposure == 2) & isequal($environment, ''TrackB'') & ($experimentday > 3)) | (($dailyexposure == 2) & isequal($environment, ''TrackA'') & ($experimentday <= 3)) )'; 

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

ca1pairs = getSameCellLinearFields(ca1f);
ca3pairs = getSameCellLinearFields(ca3f);

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
comparefunction = 'peakratevschange'; %calulate rate vs change for the peak only
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


%when COMPAREFUNCTION = 'ratevschange'
%------------------------------------------
ratebinsize = 5;
ca1combineddata = [];
ca3combineddata = [];
for n = 3:8 %which days to combine
    tmpdata = [ca1daycmp{n}.rates ca1daycmp{n}.ratechange];
    ca1combineddata = [ca1combineddata; tmpdata];
    tmpdata = [ca3daycmp{n}.rates ca3daycmp{n}.ratechange];
    ca3combineddata = [ca3combineddata; tmpdata];
end


%plot binned data
%ca1binneddata = calcbinnedtrend(ca1combineddata(:,1),ca1combineddata(:,2)./ca1combineddata(:,1),ratebinsize);
%ca3binneddata = calcbinnedtrend(ca3combineddata(:,1),ca3combineddata(:,2)./ca3combineddata(:,1),ratebinsize);
ca1binneddata = calcbinnedtrend(ca1combineddata(:,1),ca1combineddata(:,2),ratebinsize);
ca3binneddata = calcbinnedtrend(ca3combineddata(:,1),ca3combineddata(:,2),ratebinsize);


figure
plotgrouppoints(ca1binneddata.binneddata(1:10), ca1binneddata.xbins(1:10));
hold on
plot(ca1binneddata.xbins(1:10), ca1binneddata.median(1:10));
%errorbar(ca1binneddata.xbins(1:10), ca1binneddata.mean(1:10),ca1binneddata.error(1:10,1),ca1binneddata.error(1:10,2))


figure
errorbar(ca3binneddata.xbins, ca3binneddata.mean,ca3binneddata.error(:,1),ca3binneddata.error(:,2),'r.')

%Smoothed trend plotting
[ca1trenddata, ca1confbounds] = calctrend(ca1combineddata(:,1),ca1combineddata(:,2)./ca1combineddata(:,1),ratebinsize,10);
[ca3trenddata, ca3confbounds] = calctrend(ca3combineddata(:,1),ca3combineddata(:,2)./ca3combineddata(:,1),ratebinsize,10);

figure
plot(ca1trenddata(:,1),ca1trenddata(:,2))
hold on
plot(ca1trenddata(:,1),ca1confbounds(:,1),'--')
plot(ca1trenddata(:,1),ca1confbounds(:,2),'--')

figure
plot(ca3trenddata(:,1),ca3trenddata(:,2),'r')
hold on
plot(ca3trenddata(:,1),ca3confbounds(:,1),'r--')
plot(ca3trenddata(:,1),ca3confbounds(:,2),'r--')



%------------------------------------------------------
xdistvalues = [0:.01:70];
%xdistvalues = [0:.01:7];
groups = {[1 2], [3 4 5 6], [7 8 9 10 11 12 13 14]};

ca1perioddata = [];
ca3perioddata = [];
ca1fields = [];
for groupnum = 1:length(groups)
    ca1perioddata{groupnum} = [];
    ca3perioddata{groupnum} = [];
    ca1fields{groupnum} = [];
    for n = groups{groupnum} %which days to combine for end data
        %tmpdata = [ca1groups{n}(:,6:7)];
        tmpdata = [ca1daycmp{n}.rates ca1daycmp{n}.meanrate];
        ca1perioddata{groupnum} = [ca1perioddata{groupnum}; tmpdata];
        ca1fields{groupnum} = stack(ca1fields{groupnum},ca1daycmp{n}.fields);
        %tmpdata = [ca3groups{n}(:,6:7)];
        tmpdata = [ca3daycmp{n}.rates ca3daycmp{n}.meanrate];
        ca3perioddata{groupnum} = [ca3perioddata{groupnum}; tmpdata];
    end
end

ca1combineddata = [];
ca3combineddata = [];
for n = [3:8] %which days to combine for rule
    tmpdata = [ca1daycmp{n}.rates ca1daycmp{n}.ratechange];
    ca1combineddata = [ca1combineddata; tmpdata];
    tmpdata = [ca3daycmp{n}.rates ca3daycmp{n}.ratechange];
    ca3combineddata = [ca3combineddata; tmpdata];
end
ca1changedata = [ca1combineddata(:,1) ca1combineddata(:,2)];
%ca1changedata(:,2) = ca1changedata(randperm(size(ca1changedata,1)),2);
ca3changedata = [ca3combineddata(:,1) ca3combineddata(:,2)];

%[out, distances] = simulateRateModel(ca1perioddata{1}(:,1), ca1perioddata{3}(:,1), 'linear');

P0 = calcprobdist(ca1perioddata{1}(:,1),xdistvalues);
Pmid = calcprobdist(ca1perioddata{2}(:,1),xdistvalues);
Pend = calcprobdist(ca1perioddata{3}(:,1),xdistvalues);
%Pend2 = calcprobdist(ca1perioddata{4}(:,1),xdistvalues);
figure
plot(xdistvalues,P0);
hold on
plot(xdistvalues,Pmid,'g');
plot(xdistvalues,Pend,'r');
%plot(xdistvalues,Pend2,'c');

figure
plot(xdistvalues,P0);
hold on

%the first column of RATES is the mean rate, the 2nd is the initial peak
%rate, the thris is the 1st iteration, and so on...
rates = bootstrp(10000, @simulateRateDecline,ca1changedata,{ca1perioddata{1}},100);
P = [];
for i = 1:100
    P(i,1:length(xdistvalues)) = calcprobdist(rates(:,i+1),xdistvalues);
end
%plot(xdistvalues,P(2,:),'r');
plot(xdistvalues,P(20,:),'g');
plot(xdistvalues,P(100,:),'c');


%rates2 = bootstrp(5000, @simulateRateDecline2,ca1changedata,{ca1fields{1}},100);


origmeanrate = rates(:,1);
finalmeanrate = [];
for i = 1:size(rates,1)
    if (rates(i,1) > 0)
        finalmeanrate(i,1) = rates(i,1)+(rates(i,1)*((rates(i,end)-rates(i,2))/rates(i,2)));
    else
        finalmeanrate(i,1) = 0;
    end
end


%Plot CA3
CA3P0 = calcprobdist(ca3perioddata{1}(:,1),xdistvalues);
CA3Pmid = calcprobdist(ca3perioddata{2}(:,1),xdistvalues);
CA3Pend = calcprobdist(ca3perioddata{3}(:,1),xdistvalues);

figure
plot(xdistvalues,CA3P0);
hold on
plot(xdistvalues,CA3Pmid,'g');
plot(xdistvalues,CA3Pend,'r');

figure
plot(xdistvalues,CA3P0);
hold on
rates = bootstrp(10000, @simulateRateDecline,ca3changedata,{ca3perioddata{1}},100);
P = [];
for i = 1:100
    P(i,1:length(xdistvalues)) = calcprobdist(rates(:,i),xdistvalues);
end

%plot(xdistvalues,P{5},'r');
plot(xdistvalues,P(20,:),'g');
plot(xdistvalues,P(100,:),'c');

% %-------------------------------------------------











