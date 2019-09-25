%Animal selection
%-----------------------------------------------------
%animals = {'Dudley','Miles','Conley','Ten','Bond','Frank','Nine'};
%animals = {'Dudley','Miles','Conley','Alex','Bond','Frank','Nine','Ten'};
animals = {'Conley','Bond','Frank'};
%animals = {'Ten'};
%animals = {'Dudley','Miles','Conley', 'Ten'};
%animals = {'Bond','Frank','Nine'};
%animals = {'Bond'};
%-----------------------------------------------------


%Filter creation for training data
%--------------------------------------------------------
epochfilter = [];
epochfilter{1} = ['($dailyexposure == 1) & isequal($environment, ''TrackB'')'];  
%cellfilter = '(isequal($area, ''CA1'') && ($meanrate < 7))';
cellfilter = '($meanrate < 7)';
%timefilter = { {'getlinstate', '(($traj ~= -1) & (abs($velocity) >= 3))', 6},{'getriptimes','($nripples == 0)', [], 'cellfilter', '((isequal($area, ''CA3''))|(isequal($area, ''CA1'')))'} };
timefilter = { {'getriptimes','($nripples == 0)', [], 'cellfilter', '((isequal($area, ''CA3''))|(isequal($area, ''CA1'')))'} };
trainingfilter = createfilter('animal',animals,'epochs',epochfilter,'cells',cellfilter,'excludetime', timefilter);
%-----------------------------------------------------------

%create training data by calulating the linearized rates of all cells 
%--------------------------------------------
iterator = 'multicellanal';
trainingfilter = setfilteriterator(trainingfilter,iterator);
trainingfilter = setfilterfunction(trainingfilter, 'calcpopulationlinfields', {'spikes','linpos'});
trainingfilter = runfilter(trainingfilter);
%-------------------------------------------------


%Filter creation for position decoding
%--------------------------------------------------------
epochfilter = [];
epochfilter{1} = ['($dailyexposure == 1) & isequal($environment, ''TrackA'')'];  
%epochfilter{1} = ['($dailyexposure == 2) & isequal($environment, ''TrackB'')'];  
%timefilter = { {'getlinstate', '(($traj ~= -1) & (abs($velocity) >= 3))', 6},{'getriptimes','($nripples == 0)', [], 'cellfilter', '((isequal($area, ''CA3''))|(isequal($area, ''CA1'')))'} };
timefilter = { {'getlinstate', '($traj ~= -1)', 6} };

%timefilter = { {'getriptimes','($nripples >= 2)', [], 'cellfilter', '((isequal($area, ''CA3''))|(isequal($area, ''CA1'')))','minenergy',2} };
decodefilter = createfilter('animal',animals,'epochs',epochfilter,'cells',cellfilter,'excludetime', timefilter);
%-----------------------------------------------------------

%Get the binned spike counts for all cells
%-----------------------------------------------------------
timebin = .05;
iterator = 'multicellanal';
decodefilter = setfilteriterator(decodefilter,iterator);
decodefilter = setfilterfunction(decodefilter, 'getpopulationrates', {'spikes','linpos'},timebin);
decodefilter = runfilter(decodefilter);
%----------------------------------------------------------

%Decode position
%----------------------------------------------------------
animalnum = 2;
epochnum = 2;

matches = rowfind(trainingfilter(animalnum).output{1}(epochnum).index(:,[1 3 4]),decodefilter(animalnum).output{1}(epochnum).index(:,[1 3 4])); %find the matching cell indices
trainingdata = [];
spikedata = [];
decodedata = [];
indexlist = [];

%pick out all the matching cells from the training data and the
%decoding data
%traindata contains linear rates, and is n by x, where n is the
%number of cells and x is the number of spatial bins
%spikedata contains spikecounts, and is n by t, where t is the
%number of temporal bins in the data to be decoded.
for trainingcell = 1:length(matches)
    if (matches(trainingcell) > 0) %we have a match
        indexlist = [indexlist; trainingfilter(animalnum).output{1}(epochnum).index(trainingcell,:)];
        trainingdata(trainingcell,1:size(trainingfilter(animalnum).output{1}(epochnum).rates,2)) = trainingfilter(animalnum).output{1}(epochnum).rates(trainingcell,:);
        spikedata(trainingcell,1:size(decodefilter(animalnum).output{1}(epochnum).spikecounts,2)) = decodefilter(animalnum).output{1}(epochnum).spikecounts(matches(trainingcell),:);
    end
end
spikedata = double(spikedata);
trainingdata = trainingdata*timebin; %transform rates to expected number of spikes

%the decoded data contains spatial probabilities, and is x by t
decodedata = zeros(size(trainingdata,2),size(spikedata,2));
naninds = find(isnan(trainingdata(1,:)));
for t = 1:size(spikedata,2) %time  
    Tspikecount = repmat(spikedata(:,t),1,size(trainingdata,2));
    %calculate P(spikecount|x) for this timebin across all cells and all x 
    spatialprob = prod(((trainingdata.^Tspikecount)./factorial(Tspikecount)).*exp(-trainingdata),1)'; 
    spatialprob(naninds) = 0;
    spatialprob = spatialprob/sum(spatialprob);  %normalize across space
    decodedata(:,t) = spatialprob;  
end
%----------------------------------------------------------

%Play movie
%---------------------------------------------------------
figure
trajmapping = [1 1 2 2];
%maxprob = max(decodedata(:));
maxprob = .4;
xdata = {[],[]};
ydata = {[],[]};
for traj = [1:4];
    
    %subplot(2,1,trajmapping(traj));
    trajindex = find(trainingfilter(animalnum).output{1}(epochnum).traj == traj);
    xdata{trajmapping(traj)} = trainingfilter(animalnum).output{1}(epochnum).dist(trajindex);
    ydata{trajmapping(traj)} = [ydata{trajmapping(traj)}; decodedata(trajindex,1)'];
    %p(traj) = plot(trainingfilter(animalnum).output{1}(epochnum).dist(trajindex), decodedata(trajindex,1));
    %set(p(traj),'LineWidth',4);
    %hold on
    
end
for plotnum = 1:2
    subplot(2,1,plotnum);
    p(plotnum) = plot(xdata{plotnum}, sum(ydata{plotnum}));
    set(p(plotnum),'LineWidth',4);
    hold on
         
    locator{plotnum} = line([0 0],[0 maxprob],'Color',[1 0 0]);
    axis([0 180 0 maxprob]);
    hold on
end

for t = 1:size(decodedata,2)
    ydata = {[],[]};
    for traj = 1:4
        trajindex = find(trainingfilter(animalnum).output{1}(epochnum).traj == traj);
        ydata{trajmapping(traj)} = [ydata{trajmapping(traj)}; decodedata(trajindex,t)'];
        %set(p(traj),'YData',decodedata(trajindex,t));
    end
    for plotnum = 1:2
        set(p(plotnum),'YData',sum(ydata{plotnum}));
        animallocation = decodefilter(animalnum).output{1}(epochnum).dist(t);
        animaltraj = decodefilter(animalnum).output{1}(epochnum).traj(t);      
        if (animaltraj > 0)
            if (trajmapping(animaltraj) == plotnum)
                set(locator{plotnum}, 'XData', [animallocation animallocation]);
                set(locator{plotnum},'Visible','on');
            else
                set(locator{plotnum},'Visible','off');
            end
        else
            set(locator{plotnum},'Visible','off');
        end
    end
    drawnow
    pause(.2)
end
%--------------------------------------------------------









