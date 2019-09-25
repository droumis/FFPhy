%Animal selection
%-----------------------------------------------------
animals = {'Dudley','Miles','Conley','Bond','Frank','Nine','Ten'};
%animals = {'Miles','Conley','Bond','Frank','Nine','Ten'};
%animals = {'Dudley','Miles','Conley','Alex','Bond','Frank','Nine','Ten'};
%animals = {'Ten'};
%animals = {'Dudley','Miles','Conley', 'Ten'};
%animals = {'Bond','Frank','Nine'};
%animals = {'Bond'};
%animals = {'Nine'};
%-----------------------------------------------------


%Filter creation
%--------------------------------------------------------
epochfilter = [];
for i = 1:14
    
    %epochfilter{i} = ['(isequal($type,''sleep'')) & ($sleepnum == 2) & ($runbefore.exposureday == ',num2str(i),') & ($runbefore.dailyexposure == 1)'];  
    epochfilter{i} = ['($exposureday == ',num2str(i),') & ($dailyexposure == 1) & isequal($environment, ''TrackA'')'];
    %epochfilter{i} = ['($exposureday == ',num2str(i),') & ($dailyexposure == 1)'];
    %epochfilter{i} = ['(($exposureday == ',num2str(i),') & ($dailyexposure == 1) & isequal($environment, ''TrackA''))|(($exposureday == ',num2str(i),') & ($exposureday < 9) & ($dailyexposure == 1) & isequal($environment, ''TrackB''))'];
end

ca1cellfilter = '(isequal($area, ''CA1'') && ($meanrate < 7))';
ca3cellfilter = '(isequal($area, ''CA3'') && ($meanrate < 7))';

timefilter = { {'getlinstate', '(($traj ~= -1) & (abs($velocity) >= 5))', 6} };
%timefilter = { {'getlinstate', '(($traj ~= -1) & (abs($velocity) >= 3))', 6},{'getriptimes','($nripples == 0)', [], 'cellfilter', '((isequal($area, ''CA3''))|(isequal($area, ''CA1'')))','minstd',3} };
%timefilter = { {'getlinstate', '(($traj ~= -1) & (abs($velocity) < 3))', 6}, {'getriptimes','($nripples > 2)', [], 'cellfilter', '((isequal($area, ''CA3''))|(isequal($area, ''CA1'')))'}};
%timefilter = {{'getriptimes', '($nripples >= 1)', [], 'cellfilter', '(isequal($area, ''CA1''))|(isequal($area, ''CA3''))','minstd',3}};

%timefilter = {};

ca1f = createfilter('animal',animals,'epochs',epochfilter,'cells',ca1cellfilter,'excludetime', timefilter);
ca3f = createfilter('animal',animals,'epochs',epochfilter,'cells',ca3cellfilter,'excludetime', timefilter);
%-----------------------------------------------------------




%run function- single cells
%--------------------------------------------
iterator = 'singlecellanal';

ca1f = setfilteriterator(ca1f,iterator);
ca3f = setfilteriterator(ca3f,iterator);

%non occ-norm mean rate
%ca1f = setfilterfunction(ca1f, 'calclinadapt', {'spikes','linpos'},'modeltype','SpatialTemporal');
%ca3f = setfilterfunction(ca3f, 'calclinadapt', {'spikes','linpos'},'modeltype','SpatialTemporal');

%firing rate vs. velocity

%non occ-norm mean rate
%ca1f = setfilterfunction(ca1f, 'calctotalmeanrate', {'spikes'});
%ca3f = setfilterfunction(ca3f, 'calctotalmeanrate', {'spikes'});

%occ-norm mean rate
%ca1f = setfilterfunction(ca1f, 'calcoccnormmeanrate', {'spikes', 'linpos'});
%ca3f = setfilterfunction(ca3f, 'calcoccnormmeanrate', {'spikes', 'linpos'});

%peak rate
%ca1f = setfilterfunction(ca1f, 'calcpeakrate', {'spikes', 'linpos'});
%ca3f = setfilterfunction(ca3f, 'calcpeakrate', {'spikes', 'linpos'});

%out of field firing
%ca1f = setfilterfunction(ca1f, 'calcoutfieldfiring', {'spikes', 'linpos'},'thresh',.1);
%ca3f = setfilterfunction(ca3f, 'calcoutfieldfiring', {'spikes', 'linpos'},'thresh',.1);

%velocity vs. rate
ca1f = setfilterfunction(ca1f, 'calclinvelratecorr', {'spikes', 'linpos'});
ca3f = setfilterfunction(ca3f, 'calclinvelratecorr', {'spikes', 'linpos'});

%test
%ca1f = setfilterfunction(ca1f, 'calclinadapt', {'spikes', 'linpos'});
%ca3f = setfilterfunction(ca3f, 'calclinadapt', {'spikes', 'linpos'});

%ca1f = setfilterfunction(ca1f, 'funcSwitchBox', {'spikes','linpos'},'calcpeakrate','calcTrackActive');
%ca3f = setfilterfunction(ca3f, 'funcSwitchBox', {'spikes','linpos'},'calcpeakrate','calcTrackActive');

%ca1f = setfilterfunction(ca1f, 'funcSwitchBox', {'spikes','linpos'},'getmeanrate','getspikewidth');
%ca3f = setfilterfunction(ca3f, 'funcSwitchBox',{'spikes','linpos'},'getmeanrate','getspikewidth');

ca1f = runfilter(ca1f);
ca3f = runfilter(ca3f);




ca1groups = numericgroupcombine(ca1f,0);
ca3groups = numericgroupcombine(ca3f,0);

groups = {[1 2],[3 4 5 6],[7 8 9 10 11 12 13 14]};
%groups = {[1:7]};
ca1perioddata = [];
ca3perioddata = [];
for groupnum = 1:length(groups)
    ca1perioddata{groupnum} = [];
    ca3perioddata{groupnum} = [];
    for n = groups{groupnum} %which days to combine for end data
        tmpdata = [ca1groups{n}(:,:)]; 
        tmpdata(:,end+1) = n;
        ca1perioddata{groupnum} = [ca1perioddata{groupnum}; tmpdata];
        tmpdata = [ca3groups{n}(:,:)];
        tmpdata(:,end+1) = n;
        ca3perioddata{groupnum} = [ca3perioddata{groupnum}; tmpdata];
    end
    % temporary fix for linear velocity firing rate only
%    valid = find(ca1perioddata{groupnum}(:,1) ~= -10);
%    ca1perioddata{groupnum} = ca1perioddata{groupnum}(valid,:);
%    valid = find(ca3perioddata{groupnum}(:,1) ~= -10);
%    ca3perioddata{groupnum} = ca3perioddata{groupnum}(valid,:);
end


%ca1perioddata{1}(:,3) = 1;
% [B,Bint,R,Rint,stats] = regress(ca1perioddata{1}(:,1),ca1perioddata{1}(:,2:3));

'pause'
pause
%figure
%plotgroups(ca1groups);
%hold on
%plotgroups(ca3groups,'r');
%-------------------------------------------------



%population mean rate
%---------------------------------------------------

iterator = 'multicellanal';

ca1f = setfilteriterator(ca1f,iterator);
ca3f = setfilteriterator(ca3f,iterator);

figure
ca1f = setfilterfunction(ca1f, 'calcpopulationrate', {'spikes'},60);
ca3f = setfilterfunction(ca3f, 'calcpopulationrate', {'spikes'},60);
disp('ca1');
ca1f = runfilter(ca1f);
disp('ca3');
ca3f = runfilter(ca3f);

ca1groups = numericgroupcombine(ca1f);
ca3groups = numericgroupcombine(ca3f);

plotgroups(ca1groups);
hold on
plotgroups(ca3groups,'r');
%----------------------------------------------------
