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
    epochfilter{i} = ['($exposureday == ',num2str(i),') & ($dailyexposure == 1)'];
    %epochfilter{i} = ['($experimentday == ',num2str(i),') & ($dailyexposure == 1) & isequal($environment, ''TrackA'')'];
end

ca1cellfilter = '(isequal($area, ''CA1'') && ($meanrate < 7))';
ca3cellfilter = '(isequal($area, ''CA3'') && ($meanrate < 7))';

timefilter = { {'getlinstate', '(($traj ~= -1) & (abs($velocity) >= 3))', 6} };
%timefilter = {};

ca1f = createfilter('animal',animals,'epochs',epochfilter,'cells',ca1cellfilter,'excludetime', timefilter);
ca3f = createfilter('animal',animals,'epochs',epochfilter,'cells',ca3cellfilter,'excludetime', timefilter);
%-----------------------------------------------------------


%run function- single cells
%--------------------------------------------
iterator = 'singlecellanal';

ca1f = setfilteriterator(ca1f,iterator);
ca3f = setfilteriterator(ca3f,iterator);

%calculation will be thrown away, we just want the cell indices
ca1f = setfilterfunction(ca1f, 'funcSwitchBox', {'spikes','linpos'},'calctotalmeanrate');
ca3f = setfilterfunction(ca3f, 'funcSwitchBox', {'spikes','linpos'},'calctotalmeanrate');

ca1f = runfilter(ca1f);
ca3f = runfilter(ca3f);

ca1groups = numericgroupcombine(ca1f,1);
ca3groups = numericgroupcombine(ca3f,1);

%-------------------------------------------------


%Filter creation (for pairs)
%--------------------------------------------------------
epochfilter = [];

epochfilter{1} = ['( (($dailyexposure == 1) & isequal($environment, ''TrackB'') & ($experimentday > 3)) | (($dailyexposure == 1) & isequal($environment, ''TrackA'') & ($experimentday <= 3)) )'];
epochfilter{2} = ['( (($dailyexposure == 2) & isequal($environment, ''TrackB'') & ($experimentday > 3)) | (($dailyexposure == 2) & isequal($environment, ''TrackA'') & ($experimentday <= 3)) )']; 

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

ca1f = setfilterfunction(ca1f, 'funcSwitchBox', {'spikes','linpos'},'calcpeakrate','calcTrackActive');
ca3f = setfilterfunction(ca3f, 'funcSwitchBox', {'spikes','linpos'},'calcpeakrate','calcTrackActive');

%ca1f = setfilterfunction(ca1f, 'funcSwitchBox', {'spikes','linpos'},'calctotalmeanrate','calcTrackActive');
%ca3f = setfilterfunction(ca3f, 'funcSwitchBox', {'spikes','linpos'},'calctotalmeanrate','calcTrackActive');

ca1f = runfilter(ca1f);
ca3f = runfilter(ca3f);

ca1pairs = numericgroupcombine(ca1f,1);
ca3pairs = numericgroupcombine(ca3f,1);

indexcolumns = [1 2 4 5];
ca1cmp = indexmatch(ca1pairs{1},ca1pairs{2},indexcolumns,'outindex',[1 2 3 4 5]);
ca3cmp = indexmatch(ca3pairs{1},ca3pairs{2},indexcolumns,'outindex',[1 2 3 4 5]);
%--------------------------------------------------------------


% Plot the within-day change as a function of peak rate
%-------------------------------------------------------------
ca1daycmp = [];
ca3daycmp = [];
for i = 1:length(ca1groups)
    tmpdata = ca1groups{i}(:,[1 2 3 4 5]);
    ca1daycmp{i} = indexmatch(ca1cmp,tmpdata,[1 2 3 4 5]);
    tmpdata = ca3groups{i}(:,[1 2 3 4 5]);
    ca3daycmp{i} = indexmatch(ca3cmp,tmpdata,[1 2 3 4 5]);
end

activethresh = [.1 .1];
peakratenorm = 1;

%CA1------------------
figure
ca1B = [];
ca1propactive = [];
projection = [];
projection2 = [];
directionvect = [];
for i = 1:3 %length(ca1daycmp)
    subplot(3,1,i);
    arrow([ca1daycmp{i}(:,6)/peakratenorm ca1daycmp{i}(:,7)],[ca1daycmp{i}(:,8)/peakratenorm ca1daycmp{i}(:,9)],'TipAngle',0);
    hold on
    plot(ca1daycmp{i}(:,8)/peakratenorm,ca1daycmp{i}(:,9),'.');
    axis([0 60 0 1]);
    
    %[B,Bint] = regress(ca1groups{i}(:,7),[ca1groups{i}(:,6) ones(length(ca1groups{i}(:,6)),1)]);
    [B,Bint] = regress(ca1groups{i}(:,7),[ca1groups{i}(:,6)/peakratenorm]);
    B = robustfit([ca1groups{i}(:,6)/peakratenorm],ca1groups{i}(:,7),'bisquare',1,'off');
    line([0 1],[0 B],'color',[1 0 0])
    
    %unitvect = [1/sqrt((B^2)+1) B/sqrt((B^2)+1)]';
    unitvect = [1 0]';
    %unitvect = [0 1]';
    projection{i} = [ca1daycmp{i}(:,6)/peakratenorm ca1daycmp{i}(:,7)]*unitvect;
    projection2{i} = [ca1daycmp{i}(:,8)/peakratenorm ca1daycmp{i}(:,9)]*unitvect;
    
    directionvect = [directionvect; [projection{i} projection2{i}-projection{i}]];
    
    ca1B(i) = B(1);
    ca1propactive(i) = length(find((ca1daycmp{i}(:,6)/peakratenorm >= activethresh(1)) & (ca1daycmp{i}(:,7) > activethresh(2))))/length(ca1daycmp{i}(:,6));
    %plot(ca1daycmp{i}(:,6),ca1daycmp{i}(:,7),'.');
    
end
directionvect = sortrows(directionvect);
smooth = smoothvect(directionvect(:,2),gaussian(20,50));
figure
plot(directionvect(:,1),smooth);

%CA3------------------------------
figure
ca3B = [];
ca3projection = [];
ca3projection2 = [];
ca3propactive = [];
ca3directionvect = [];
for i = 1:3 %length(ca3daycmp)
    subplot(3,1,i);
    arrow([ca3daycmp{i}(:,6)/peakratenorm ca3daycmp{i}(:,7)],[ca3daycmp{i}(:,8)/peakratenorm ca3daycmp{i}(:,9)],'TipAngle',0);
    hold on
    plot(ca3daycmp{i}(:,8)/peakratenorm,ca3daycmp{i}(:,9),'r.');
    %[B,Bint] = regress(ca3groups{i}(:,7),[ca3groups{i}(:,6) ones(length(ca3groups{i}(:,6)),1)]);
    [B,Bint] = regress(ca3groups{i}(:,7),[ca3groups{i}(:,6)/peakratenorm]);
    %B = robustfit([ca3groups{i}(:,6)/peakratenorm],ca3groups{i}(:,7),'bisquare',1,'off');
    ca3B(i) = B(1);
    %unitvect = [1/sqrt((B^2)+1) B/sqrt((B^2)+1)]';
    unitvect = [1 0]';
    %unitvect = [0 1]';
    ca3projection{i} = [ca3daycmp{i}(:,6)/peakratenorm ca3daycmp{i}(:,7)]*unitvect;
    ca3projection2{i} = [ca3daycmp{i}(:,8)/peakratenorm ca3daycmp{i}(:,9)]*unitvect;
    
    ca3directionvect = [ca3directionvect; [ca3projection{i} ca3projection2{i}-ca3projection{i}]];
    
    ca3propactive(i) = length(find((ca3daycmp{i}(:,6) >= activethresh(1)) & (ca3daycmp{i}(:,7) > activethresh(2))))/length(ca3daycmp{i}(:,6));
    %plot(ca1daycmp{i}(:,6),ca1daycmp{i}(:,7),'.');
    axis([0 60 0 1]);
end
ca3directionvect = sortrows(ca3directionvect);
smooth = smoothvect(ca3directionvect(:,2),gaussian(20,50));
figure
plot(ca3directionvect(:,1),smooth,'r');

%----------------------------------------------










