%Animal selection
%-----------------------------------------------------
animals = {'Dudley','Miles','Conley','Bond','Frank','Nine','Ten'};
%animals = {'Dudley','Miles','Conley','Alex','Bond','Frank','Nine','Ten'};
%animals = {'Ten'};
%animals = {'Dudley','Miles','Conley', 'Ten'};
%animals = {'Bond','Frank','Nine'};
%animals = {'Bond'};
%-----------------------------------------------------


%Filter creation
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

ca1f = setfilterfunction(ca1f, 'funcSwitchBox', {'spikes','linpos'},'calcpeakrate','calcTrackActive');
ca3f = setfilterfunction(ca3f, 'funcSwitchBox', {'spikes','linpos'},'calcpeakrate','calcTrackActive');

ca1f = runfilter(ca1f);
ca3f = runfilter(ca3f);

ca1groups = numericgroupcombine(ca1f,1);
ca3groups = numericgroupcombine(ca3f,1);

%-------------------------------------------------


%plot peak vs prop active
%-------------------------------------------
%groups = {[1 2 3],[4 5 6 7 8 9],[10 11 12 13 14]};
groups = {[1],[2],[3],[4],[5],[6],[7]};
ca1col =  pink(length(groups)+2);
ca3col = pink(length(groups)+2);
ca1slope = [];
ca3slope = [];
ca1plot = [];
ca3plot = [];
ca1Bint = [];
ca3Bint = [];
for i = 1:length(groups)
    ca1plot{i} = [];
    ca3plot{i} = [];
    for j = 1:length(groups{i})
        ca1plot{i} = [ca1plot{i}; ca1groups{groups{i}(j)}(:,[6:7])];
        ca3plot{i} = [ca3plot{i}; ca3groups{groups{i}(j)}(:,[6:7])];
    end
end
figure
for i = 1:length(ca1plot)
    plot(ca1plot{i}(:,1),ca1plot{i}(:,2),'.','color',ca1col(i,:));
    hold on
    [B,ca1Bint(i,1:2)] = regress(ca1plot{i}(:,2),[ca1plot{i}(:,1)]);   
    %B = robustfit([ca1plot{i}(:,1)],ca1plot{i}(:,2),'fair',1,'off');
    line([0 60],[0 B*60],'color',ca1col(i,:))
    ca1slope(i) = B;
end
figure
for i = 1:length(ca3plot)
    plot(ca3plot{i}(:,1),ca3plot{i}(:,2),'.','color',ca3col(i,:));
    hold on
    [B,ca3Bint(i,1:2)] = regress(ca3plot{i}(:,2),[ca3plot{i}(:,1)]);   
    %B = robustfit([ca3plot{i}(:,1)],ca3plot{i}(:,2),'fair',1,'off');
    line([0 60],[0 B*60],'color',ca3col(i,:))
    ca3slope(i) = B;
end
%------------------------------------------------------


%plot probability distribution plots
%-------------------------------------------
groups = {[1 2 3], [4 5 6 7 8],[9 10 11 12 13]};
%groups = {[1],[2],[3],[4],[5],[6],[7]};
ca1col =  pink(length(groups)+2);
ca3col = pink(length(groups)+2);
ca1plot = [];
ca3plot = [];

for i = 1:length(groups)
    ca1plot{i} = [];
    ca3plot{i} = [];
    for j = 1:length(groups{i})
        ca1plot{i} = [ca1plot{i}; ca1groups{groups{i}(j)}(:,[6:7])];
        ca3plot{i} = [ca3plot{i}; ca3groups{groups{i}(j)}(:,[6:7])];
    end
end

figure
ca1axeshandle = axes;
for i = 1:length(ca1plot)   
    ca1plot{i} = flipud(sortrows(ca1plot{i},1));
    numcells = size(ca1plot{i},1);
    ca1plot{i}(:,3) = 1;
    ca1plot{i}(:,3) = cumsum(ca1plot{i}(:,3))/numcells;
    [ca1X, Xindex] = unique(ca1plot{i}(:,1));
    ca1Y = ca1plot{i}(Xindex,3);
    
    ca1Xi = [min(ca1X):.01:max(ca1X)]';
    ca1Yi = interp1(ca1X,ca1Y,ca1Xi);
    
    ca1Yslope = diff(ca1Yi);
    ca1Xslope = ca1Xi(2:end);
    ca1Yslope = -smoothvect(ca1Yslope,gaussian(150,400));
   
    plot(ca1Xslope, ca1Yslope, 'color', ca1col(i,:));
    %plot(ca1Xi,ca1Yi,'color',ca1col(i,:));
    
    hold on
end

figure
ca3axeshandle = axes;
for i = 1:length(ca3plot)   
    ca3plot{i} = flipud(sortrows(ca3plot{i},1));
    numcells = size(ca3plot{i},1);
    ca3plot{i}(:,3) = 1;
    ca3plot{i}(:,3) = cumsum(ca3plot{i}(:,3))/numcells;
    [ca3X, Xindex] = unique(ca3plot{i}(:,1));
    ca3Y = ca3plot{i}(Xindex,3);
    
    ca3Xi = [min(ca3X):.01:max(ca3X)]';
    ca3Yi = interp1(ca3X,ca3Y,ca3Xi);
   
    ca3Yslope = diff(ca3Yi);
    ca3Xslope = ca3Xi(2:end);
    ca3Yslope = -smoothvect(ca3Yslope,gaussian(150,400));
   
    plot(ca3Xslope, ca3Yslope, 'color', ca3col(i,:));
    %plot(ca3Xi,ca3Yi,'color',ca3col(i,:));
    
    hold on
end

%------------------------------------------------------


        









