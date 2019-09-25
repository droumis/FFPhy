
%Animal selection
%-----------------------------------------------------
animals = {'Dudley','Miles','Conley','Bond','Frank','Nine','Ten'};
%animals = {'Dudley','Miles','Conley','Bond','Frank','Nine'};
%animals = {'Ten'};
%animals = {'Dudley','Miles','Conley', 'Ten'};
%animals = {'Frank'};
%animals = {'Bond'};
%-----------------------------------------------------


%Filter creation
%--------------------------------------------------------
epochfilter = [];
epochfilter = [];
for i = 1:14
    epochfilter{i} = ['($exposureday == ',num2str(i),') & ($dailyexposure == 1)'];
    %epochfilter{i} = ['(($exposureday == ',num2str(i),') & ($dailyexposure == 1) & isequal($environment, ''TrackA''))|(($exposureday == ',num2str(i),') & ($exposureday < 9) & ($dailyexposure == 1) & isequal($environment, ''TrackB''))'];
end

%epochfilter{1} = ['($exposure < 2)'];
%epochfilter{2} = ['($exposure == 6)'];

cellpairfilter = {'allcomb', '(isequal($area, ''CA3'') && ($meanrate < 7) && ($meanrate > .2))', '(isequal($area, ''CA3'') && ($meanrate < 7)) && ($meanrate > .2)'};

timefilter = {{'getriptimes', '($nripples > 2)', [], 'cellfilter', '(isequal($area, ''CA1''))','minenergy',2}};


iterator = 'singlecellanal';

f = createfilter('animal',animals,'epochs',epochfilter,'cellpairs',cellpairfilter,'excludetimefilter', timefilter, 'iterator', iterator);

f = setfilterfunction(f, 'calcxcorrmeasures', {'spikes', 'linpos'}, 'edgespikes', 1);
f = runfilter(f);

groups = numericgroupcombine(f,0);
figure
plotgroups(groups,'',1);