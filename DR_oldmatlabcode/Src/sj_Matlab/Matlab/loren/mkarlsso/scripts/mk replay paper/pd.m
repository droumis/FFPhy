%Animal selection
%-----------------------------------------------------
%animals = {'Dudley','Miles','Conley','Ten','Bond','Frank','Nine'};
%animals = {'Dudley','Miles','Conley','Alex','Bond','Frank','Nine','Ten'};


%animals = {'Conley','Bond','Frank','Miles','Nine', 'Ten'};
%animals = {'Miles','Nine', 'Ten'};
animals = {'Frank'};
%animals = {'Ten'};
%animals = {'Dudley','Miles','Conley', 'Ten'};
%animals = {'Bond','Frank','Nine'};
%animals = {'Bond'};
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

%Filter creation for place field data
%placefieldfilter = createfilter('animal',animals,'epochs',epochfilter,'cells',cellfilter,'excludetime', timefilter);
%-----------------------------------------------------------

%placefieldfilter = setfilteriterator(trainingfilter,iterator);
% create a list of all of the place fields without a peak rate threshold
%placefieldfilter = setfilterfunction(trainingfilter, 'calcpopulationlinfields', {'spikes','linpos'},2, 0);

%placefieldfilter = runfilter(placefieldfilter);
%Filter creation for position decoding
%--------------------------------------------------------
epochfilter = [];
cellfilter = '(($meanrate < 7))';
%epochfilter{1} = ['($dailyexposure == 1) & isequal($environment, ''TrackB'')'];  
%epochfilter{1} = ['($dailyexposure == 2) & isequal($environment, ''TrackB'')'];  
epochfilter{1} = ['($dailyexposure == 1) & isequal($environment, ''TrackA'') & ($exposureday > 3)'];  
%epochfilter{1} = ['(isequal($type,''sleep'')) & isequal($runbefore.environment, ''TrackA'') & ($runbefore.exposureday > 3)'];  
%epochfilter{1} = ['(isequal($type,''sleep'')) & isequal($runbefore.environment, ''TrackB'') & isequal($runafter.environment, ''TrackB'')'];  
%epochfilter{1} = ['(isequal($type,''sleep'')) & isequal($runbefore.environment, ''TrackB'') & isequal($runafter.environment, ''TrackA'')'];  


%timefilter = { {'getlinstate', '(($traj ~= -1) & (abs($velocity) >= 3))', 6},{'getriptimes','($nripples == 0)', [], 'cellfilter', '((isequal($area, ''CA3''))|(isequal($area, ''CA1'')))'} };
%timefilter = { {'getlinstate', '(abs($velocity) < 3)', 6},{'getriptimes','($nripples >= 2)', [], 'cellfilter', '((isequal($area, ''CA3''))|(isequal($area, ''CA1'')))'} };
%timefilter = { {'getriptimes','($nripples >= 2)', [], 'cellfilter', '((isequal($area, ''CA3''))|(isequal($area, ''CA1'')))','minenergy',2} };
%timefilter = { {'getlinstate', '(abs($velocity) < 3)', 6} };
%timefilter = {{'get2dstate', '(($immobilitytime > 10)'}};
%timefilter = {{'get2dstate', '((abs($velocity) < 3))'}};

timefilter = {};

decodefilter = createfilter('animal',animals,'epochs',epochfilter,'cells',cellfilter,'excludetime', timefilter);
%-----------------------------------------------------------

%Get the binned spike counts for all cells
%-----------------------------------------------------------

cellcountthresh = 5;
iterator = 'multicellanal';
decodefilter = setfilteriterator(decodefilter,iterator);
decodefilter = setfilterfunction(decodefilter, 'getpopulationevents2', {'spikes','linpos','pos','ripples','cellinfo'},cellcountthresh);
%decodefilter = setfilterfunction(decodefilter, 'getpopulationevents_nolinpos', {'spikes','pos', 'ripples','cellinfo'},cellcountthresh);
%timebin = .2;
%decodefilter = setfilterfunction(decodefilter, 'getpopulationrates', {'spikes','linpos'},timebin);
decodefilter = runfilter(decodefilter);
%----------------------------------------------------------



%Decode position
%----------------------------------------------------------
% animalnum = 2;
% epochnum = 2;
% 
% matches = rowfind(trainingfilter(animalnum).output{1}(epochnum).index(:,[1 3 4]),decodefilter(animalnum).output{1}(epochnum).index(:,[1 3 4])); %find the matching cell indices
% trainingdata = [];
% spikedata = [];
% decodedata = [];
% indexlist = [];
% 
% %pick out all the matching cells from the training data and the
% %decoding data
% %traindata contains linear rates, and is n by x, where n is the
% %number of cells and x is the number of spatial bins
% %spikedata contains spikecounts, and is n by t, where t is the
% %number of temporal bins in the data to be decoded.
% for trainingcell = 1:length(matches)
%     if (matches(trainingcell) > 0) %we have a match
%         indexlist = [indexlist; trainingfilter(animalnum).output{1}(epochnum).index(trainingcell,:)];
%         trainingdata(trainingcell,1:size(trainingfilter(animalnum).output{1}(epochnum).rates,2)) = trainingfilter(animalnum).output{1}(epochnum).rates(trainingcell,:);
%         spikedata(trainingcell,1:size(decodefilter(animalnum).output{1}(epochnum).spikecounts,2)) = decodefilter(animalnum).output{1}(epochnum).spikecounts(matches(trainingcell),:);
%     end
% end
% spikedata = double(spikedata);
% trainingdata = trainingdata*timebin; %transform rates to expected number of spikes
% 
% %the decoded data contains spatial probabilities, and is x by t
% decodedata = zeros(size(trainingdata,2),size(spikedata,2));
% naninds = find(isnan(trainingdata(1,:)));
% for t = 1:size(spikedata,2) %time  
%     Tspikecount = repmat(spikedata(:,t),1,size(trainingdata,2));
%     %calculate P(spikecount|x) for this timebin across all cells and all x 
%     spatialprob = prod(((trainingdata.^Tspikecount)./factorial(Tspikecount)).*exp(-trainingdata),1)'; 
%     spatialprob(naninds) = 0;
%     spatialprob = spatialprob/sum(spatialprob);  %normalize across space
%     decodedata(:,t) = spatialprob;  
% end
% %----------------------------------------------------------

%save /data/loren/mkarlsso/Replay/trackBreplaydata.mat
%keyboard
viewdecodeevent([1 7], trainingfilter, decodefilter);  %For loren: this is where you view the events

'pause'
pause

%list = [1 1;1 2;3 2;3 3;3 4;3 5;3 6;3 7;3 8;3 9;3 10;3 11];
list = [1 1;1 2;2 1;2 2;2 3;2 4;2 5;2 6;2 7;2 8;3 2;3 3;3 4;3 5;3 6;3 7;3 8;3 9;3 10;3 11];
%list = [1 1;1 2;3 3;3 4;3 5;3 6;3 7;3 8;3 9;3 10;3 11];
%list = [1 1;1 2;3 2;3 3;3 4];
%list = [4 1;5 1;6 1;6 2;6 3;6 4];
out = [];
% allow for updates of the decode filter to add ripple stats
global decodefilter
for i = 1:size(list,1)
    % get the place field rates for each activated cell in each event
%    out{i} =  calcepochreplayspikestats(list(i,:), 2, placefieldfilter, decodefilter, .05);
    out = [out; calcepochreplaystats(list(i,:), 1, trainingfilter, decodefilter)];
    %out = [out; calcepochreplaystatsbothtrainers(list(i,:), trainingfilter, decodefilter)];
    %out = [out; calcepochripplestats(list(i,:), 1, trainingfilter, decodefilter)];
    %out = [out; calcepochactivationstats(list(i,:), 1, trainingfilter, decodefilter)]; %[numcellsactivated percentofactivecells totalspikesacrosscells]
%     [tmpcounts, totals] = calcepochactivation(list(i,:), trainingfilter, decodefilter);
%     if ~isempty(tmpcounts)
%         tmpcounts(:,1) = tmpcounts(:,1)/totals(1);
%         tmpcounts(:,2) = tmpcounts(:,2)/totals(2);
%         out = [out;tmpcounts];
%     end
end
%---------------------------------------------------------

'pause'
pause


% if we just got the list of place field rates, we need to calculate the
% interpotaled rate from 0 to 1 representing the beginning and end of each
% event
nbins = 100;
peakrate = [];
flpeak = [];
localrate = [];
fllocal = [];
relloc = [];
localr = [];
for i = 1:length(out)
    for e = 1:length(out{i})
	len = length(out{i}{e}(:,1));
	xi = linspace(1, len, nbins);
        ptmp = out{i}{e}(:,1); 
        ltmp = out{i}{e}(:,2); 
	relloc = [relloc ; [1:len]' ./ len];
	localr = [localr ; ltmp];
	peakrate = [peakrate ; interp1(1:len, ptmp', xi, 'nearest')];
	localrate = [localrate ; interp1(1:len, ltmp', xi, 'nearest')];
	if (isfinite(ptmp(1)) & isfinite(ptmp(end)))
	    flpeak = [flpeak; ptmp(1) ptmp(end)];
	end
	if (isfinite(ltmp(1)) & isfinite(ltmp(end)))
	    fllocal = [fllocal; ltmp(1) ltmp(end)];
	end
    end
end



immobile = out((out(:,4)>1),:);
%immobile = out((out(:,4)>2.4),:);
mobile = out((out(:,4)<1),:);
%immobile = out(((out(:,4)>0)&(out(:,4)<10)),:);
%[P,z] = zproptest([sum(immobile(:,3) < .05) sum(mobile(:,3) < .05)],[length(immobile) length(mobile)])
[P,z] = ztestprop2([sum(immobile(:,3) < .05) length(immobile)],[sum(mobile(:,3) < .05) length(mobile)])



% 
% %Play movie
% %---------------------------------------------------------
% figure
% trajmapping = [1 1 2 2];
% %maxprob = max(decodedata(:));
% maxprob = .4;
% xdata = {[],[]};
% ydata = {[],[]};
% for traj = [1:4];
%     
%     %subplot(2,1,trajmapping(traj));
%     trajindex = find(trainingfilter(animalnum).output{1}(epochnum).traj == traj);
%     xdata{trajmapping(traj)} = trainingfilter(animalnum).output{1}(epochnum).dist(trajindex);
%     ydata{trajmapping(traj)} = [ydata{trajmapping(traj)}; decodedata(trajindex,1)'];
%     %p(traj) = plot(trainingfilter(animalnum).output{1}(epochnum).dist(trajindex), decodedata(trajindex,1));
%     %set(p(traj),'LineWidth',4);
%     %hold on
%     
% end
% for plotnum = 1:2
%     subplot(2,1,plotnum);
%     p(plotnum) = plot(xdata{plotnum}, sum(ydata{plotnum}));
%     set(p(plotnum),'LineWidth',4);
%     hold on
%          
%     locator{plotnum} = line([0 0],[0 maxprob],'Color',[1 0 0]);
%     axis([0 180 0 maxprob]);
%     hold on
% end
% 
% for t = 1:size(decodedata,2)
%     ydata = {[],[]};
%     for traj = 1:4
%         trajindex = find(trainingfilter(animalnum).output{1}(epochnum).traj == traj);
%         ydata{trajmapping(traj)} = [ydata{trajmapping(traj)}; decodedata(trajindex,t)'];
%         %set(p(traj),'YData',decodedata(trajindex,t));
%     end
%     for plotnum = 1:2
%         set(p(plotnum),'YData',sum(ydata{plotnum}));
%         animallocation = decodefilter(animalnum).output{1}(epochnum).dist(t);
%         animaltraj = decodefilter(animalnum).output{1}(epochnum).traj(t);      
%         if (animaltraj > 0)
%             if (trajmapping(animaltraj) == plotnum)
%                 set(locator{plotnum}, 'XData', [animallocation animallocation]);
%                 set(locator{plotnum},'Visible','on');
%             else
%                 set(locator{plotnum},'Visible','off');
%             end
%         else
%             set(locator{plotnum},'Visible','off');
%         end
%     end
%     drawnow
%     pause(.2)
% end
% %--------------------------------------------------------
