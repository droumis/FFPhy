function h = viewdecodeevent2(varargin)
%H = viewdecodeevent2(index, trainingfilter,decodefilter)
%  
%This function creates a playback gui for viewing the output of the
%adaptive filter. It assumes a very specific structure to the adaptest 
%and lindistpos data structures.
%
%index- [day epoch tetrode cellnum]
%datapath- example '/data13/mkarlsso/Con/'
%fileprefix- any animal specific prefix to the data files, example 'con'.
%                Enter '' if none.
%lowercase- any animal specific prefix to the variable names, example
%                'con'. Enter '' if none.
if (~ischar(varargin{1}))  %figure switchbox
    %create the figure
    cellindex = varargin{1};
    h = figure('MenuBar','none','Color',get(0,'DefaultUicontrolBackgroundColor'), ... 
        'Tag','viewdecodeevent2','NumberTitle','off','Resize','on','Name',['Place cell movie viewer: cell ',num2str(cellindex)], ... 
        'CloseRequestFcn','viewdecodeevent2(''viewdecodeevent2_CloseRequestFnc'',gcbo,guidata(gcbo));', ...
        'WindowStyle','normal');
     %try
         %fill the figure
         fill_figure(varargin{1},varargin{2},varargin{3},h);
%      catch
%          %something crashed during the fill. Display error and delete
%          %figure. Output an empty handle.
%          disp(lasterr)
%          delete(h);
%          h = [];
%      end
else
    %call an internal function
    feval(varargin{:});
    h = [];
end
%---------------------------------------------------------------
function fill_figure(index, trainingfilter, decodefilter, fighandle)
% this is the fill function.  It sets up the figure and all of the objects
% within

frameperiod = .1;   %this is the default run speed
binsize = .015; %default temporal bin
currentbin = 1;
eventindex = 1;
trainingindex = 1;
animalnum = index(1);
epochnum = index(2);
eventtime = decodefilter(animalnum).output{1}(epochnum).eventtime(eventindex,1);
eventtimeend = decodefilter(animalnum).output{1}(epochnum).eventtime(eventindex,2);
homesegmentlength = trainingfilter(animalnum).output{1}(epochnum).homesegmentlength;
maxheight = .4;

figpos = get(fighandle,'Position');   

%create the uicontrols
b(1) = uicontrol('Tag','playbutton','Style','pushbutton','Units','pixel','Position',[5 105 50 20],'String','Play','Callback','viewdecodeevent2(''playbutton_Callback'',gcbo,guidata(gcbo));');
b(2) = uicontrol('Tag','downspeedbutton','Style','pushbutton','Units','pixel','Position',[165 105 20 20],'String','<','Callback','viewdecodeevent2(''downspeedbutton_Callback'',gcbo,guidata(gcbo));');
b(3) = uicontrol('Tag','upspeedbutton','Style','pushbutton','Units','pixel','Position',[185 105 20 20],'String','>','Callback','viewdecodeevent2(''upspeedbutton_Callback'',gcbo,guidata(gcbo));');
b(4) = uicontrol('Tag','downeventbutton','Style','pushbutton','Units','pixel','Position',[165 70 20 20],'String','<','Callback','viewdecodeevent2(''downeventbutton_Callback'',gcbo,guidata(gcbo));');
b(5) = uicontrol('Tag','upeventbutton','Style','pushbutton','Units','pixel','Position',[185 70 20 20],'String','>','Callback','viewdecodeevent2(''upeventbutton_Callback'',gcbo,guidata(gcbo));');
b(6) = uicontrol('Tag','downtrainingbutton','Style','pushbutton','Units','pixel','Position',[55 70 20 20],'String','<','Callback','viewdecodeevent2(''downtrainingbutton_Callback'',gcbo,guidata(gcbo));');
b(7) = uicontrol('Tag','uptrainingbutton','Style','pushbutton','Units','pixel','Position',[75 70 20 20],'String','>','Callback','viewdecodeevent2(''uptrainingbutton_Callback'',gcbo,guidata(gcbo));');
b(8) = uicontrol('Tag','makefigurebutton','Style','pushbutton','Units','pixel','Position',[5 40 50 20],'String','Figure','Callback','viewdecodeevent2(''makefigurebutton_Callback'',gcbo,guidata(gcbo));');
b(9) = uicontrol('Tag','makeseriesbutton','Style','pushbutton','Units','pixel','Position',[57 40 50 20],'String','Series','Callback','viewdecodeevent2(''makeseriesbutton_Callback'',gcbo,guidata(gcbo));');


e(1) = uicontrol('Tag','timedisplay','Style','edit','Units','pixel','Position',[55 105 60 20],'String',num2str(eventtime),'HorizontalAlignment','left');
e(2) = uicontrol('Tag','speeddisplay','Style','edit','Units','pixel','Position',[115 105 50 20],'String',num2str(frameperiod),'HorizontalAlignment','left');
e(3) = uicontrol('Tag','eventdisplay','Style','edit','Units','pixel','Position',[115 70 50 20],'String',num2str(eventindex),'HorizontalAlignment','left');
e(4) = uicontrol('Tag','trainingdisplay','Style','edit','Units','pixel','Position',[5 70 50 20],'String',num2str(trainingindex),'HorizontalAlignment','left');


t(1) = uicontrol('Tag','timetext','Style','text','Units','pixel','Position',[55 125 50 12],'String','Time');
t(2) = uicontrol('Tag','speedtext','Style','text','Units','pixel','Position',[115 125 80 12],'String','Frame period');
t(3) = uicontrol('Tag','eventtext','Style','text','Units','pixel','Position',[115 90 80 12],'String','Event');
t(4) = uicontrol('Tag','trainingtext','Style','text','Units','pixel','Position',[5 90 80 12],'String','Training data');

a(1) = axes('Tag','locationaxes1','Position',[.04 .55 .40 .245]);
a(2) = axes('Tag','locationaxes2','Position',[.50 .74 .49 .245]);
a(3) = axes('Tag','locationaxes3','Position',[.50 .4 .49 .245]);
a(4) = axes('Tag','rasteraxes','Units','pixel','Position',[310 5 figpos(3)-320 120]);

displayeventcells = {'None'};
l(1) = uicontrol('Style','listbox','Tag','listbox1','Units','pixel','FontUnits','pixels','FontSize',11,'BackgroundColor',[1 1 1],'ForegroundColor',[0 0 0], ...
'Position',[210 5 70 120],'String',displayeventcells,'Max',1,'Min',1,'Callback','viewdecodeevent2(''list1_Callback'',gcbo,guidata(gcbo));', 'Value',[1]);

numplots = 2;
matches = rowfind(trainingfilter(animalnum).output{trainingindex}(epochnum).index(:,[1 3 4]),decodefilter(animalnum).output{1}(epochnum).index(:,[1 3 4])); %find the matching cell indices
eventcellsactive = [];
eventcellpeaks = [];
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

startevent = decodefilter(animalnum).output{1}(epochnum).eventtime(eventindex,1);
endevent = decodefilter(animalnum).output{1}(epochnum).eventtime(eventindex,2);
timebins = startevent:binsize:endevent;
for trainingcell = 1:length(matches)
    if (matches(trainingcell) > 0) %we have a match
        indexlist = [indexlist; trainingfilter(animalnum).output{trainingindex}(epochnum).index(trainingcell,:)];
        trainingdata = [trainingdata; trainingfilter(animalnum).output{trainingindex}(epochnum).rates(trainingcell,:)];      
        [tmppeak, tmppeakind] = max(trainingfilter(animalnum).output{trainingindex}(epochnum).rates(trainingcell,:));
        tmppeakdist = trainingfilter(animalnum).output{trainingindex}(epochnum).dist(tmppeakind);
        tmpspiketimes = decodefilter(animalnum).output{1}(epochnum).eventdata(eventindex).spiketimes(find(decodefilter(animalnum).output{1}(epochnum).eventdata(eventindex).cellindex == matches(trainingcell)));
        if ~isempty(tmpspiketimes)
            eventcellsactive = [eventcellsactive matches(trainingcell)];
            eventcellpeaks = [eventcellpeaks tmppeakdist];
        end
        spikebins = lookup(tmpspiketimes,timebins);
        spikecount = zeros(1,length(timebins));
        for i = 1:length(spikebins)
            spikecount(spikebins(i)) = spikecount(spikebins(i))+1;
        end
        spikedata = [spikedata; spikecount];
    end
end
trainingdata = trainingdata*binsize; %transform rates to expected number of spikes

if ~isempty(eventcellsactive)
    eventcells = [eventcellsactive' eventcellpeaks'];
    eventcells = sortrows(eventcells,2);
    eventcellsactive = eventcells(:,1);
    eventcellpeaks = eventcells(:,2);
else
    eventcells = [];
end
axes(a(4));
tickplot = [nan nan];
for i = 1:size(eventcells,1)
    displayeventcells{i+1} = ['Cell ',num2str(eventcells(i))];
    tmpspiketimes = decodefilter(animalnum).output{1}(epochnum).eventdata(eventindex).spiketimes(find(decodefilter(animalnum).output{1}(epochnum).eventdata(eventindex).cellindex == eventcells(i)));
    tmpspiketimes(:,2) = i;
    tickplot = [tickplot; tmpspiketimes];
end
rasterplot = plot(tickplot(:,1),tickplot(:,2),'k+');
hold on;
set(a(4),'YDir','reverse');
%set(a(4), 'XLim', [median(tickplot(2:end,1))-.5 median(tickplot(2:end,1))+.5]);
etime = decodefilter(animalnum).output{1}(epochnum).eventtime(eventindex,:);
set(a(4), 'XLim', [etime(1) - .2 etime(2) + .2]);
set(l(1), 'String',displayeventcells);

%the decoded data contains spatial probabilities, and is x by t
decodedata = zeros(size(trainingdata,2),size(spikedata,2));
naninds = find(isnan(trainingdata(1,:)));
for t = 1:size(spikedata,2) %time  
    Tspikecount = repmat(spikedata(:,t),1,size(trainingdata,2));
    %calculate P(spikecount|x) for this timebin across all cells and all x 
    spatialprob = prod(((trainingdata.^Tspikecount)./factorial(Tspikecount)).*exp(-trainingdata),1)'; 
    spatialprob(naninds) = 0;
    if (sum(spikedata(:,t)) ~= 0)
	spatialprob = spatialprob/sum(spatialprob);  %normalize across space
    else
	spatialprob(:) = 0;
    end
    decodedata(:,t) = spatialprob;  
end

trajmapping = [1 1 2 2];
xdata = {[],[]};
ydata = {[],[]};
plotxdata = {[],[],[]};
plotydata = {[],[],[]};

%combine the outbound and inbound trajectories
for traj = [1:4];    
    trajindex = find(trainingfilter(animalnum).output{trainingindex}(epochnum).traj == traj);
    xdata{trajmapping(traj)} = trainingfilter(animalnum).output{trainingindex}(epochnum).dist(trajindex);
    ydata{trajmapping(traj)} = stack(ydata{trajmapping(traj)}, decodedata(trajindex,1)');    
end

%combine data for the home arm
for traj = 1:2
    ydata{traj} = sum(ydata{traj});
    homeindex = find(xdata{traj} < homesegmentlength);
    outerindex = find(xdata{traj} >= homesegmentlength);
    plotxdata{1} = xdata{traj}(homeindex);
    plotxdata{traj+1} = xdata{traj}(outerindex);
    plotydata{1} = [plotydata{1}; ydata{traj}(homeindex)];
    plotydata{traj+1} = ydata{traj}(outerindex);
end
plotydata{1} = sum(plotydata{1});

%plot the three plots
for plotnum = 1:3
    axes(a(plotnum))
    p(plotnum) = plot(plotxdata{plotnum}, plotydata{plotnum});
    set(p(plotnum),'LineWidth',4);
    hold on
         
    locator{plotnum} = line([0 0],[0 maxheight],'Color',[1 0 0]);
    try
        animallocation = decodefilter(animalnum).output{1}(epochnum).eventdist(eventindex);
        animaltraj = decodefilter(animalnum).output{1}(epochnum).eventtraj(eventindex);
    catch
        animallocation = 0;
        animaltraj = 0;
    end
    if (animallocation < homesegmentlength)
        set(locator{1}, 'XData', [animallocation animallocation]);
        set(locator{1},'Visible','on');       
    else
        if (animaltraj > 0)
            if (trajmapping(animaltraj)+1 == plotnum)
                set(locator{plotnum}, 'XData', [animallocation animallocation]);
                set(locator{plotnum},'Visible','on');
            else
                set(locator{plotnum},'Visible','off');
            end
        else
            set(locator{plotnum},'Visible','off');
        end
    end
    if (plotnum == 1)
        axis([min(plotxdata{plotnum}) max(plotxdata{plotnum}) 0 maxheight]);
    else
        axis([min(plotxdata{plotnum}) 185 0 maxheight]);
    end
    hold on
end

%set up the place field plots
for plotnum = 1:3
    axes(a(plotnum))
    tmpratedata = zeros(1,length(plotxdata{plotnum}));
    placeplot(plotnum) = plot(plotxdata{plotnum},tmpratedata,'g');
    set(placeplot(plotnum),'Visible','off');
end




%creat the timer object that is called during playback
plottimer = timer('Tag',['plottimer',num2str(fighandle)],'Period',frameperiod,'TimerFcn',['viewdecodeevent2(''plotTimerFcn'',timerfind(''Tag'',''plottimer',num2str(fighandle),'''));'],'TasksToExecute',1000000,'ExecutionMode','fixedDelay','UserData',[fighandle 1]);

%save important data to the figure 
fighandles = guihandles(fighandle);
fighandles.fighandle = fighandle;
fighandles.placeplot = placeplot;
fighandles.a = a;
fighandles.p = p;
fighandles.rasterplot = rasterplot;
fighandles.trajmapping = trajmapping;
fighandles.timebins = timebins;
fighandles.plottimer = plottimer;
fighandles.binsize = binsize;
fighandles.currentbin = currentbin;
fighandles.displayeventcells = displayeventcells;
fighandles.eventcellsactive = eventcellsactive;
fighandles.eventcellpeaks = eventcellpeaks;
fighandles.index = index;
fighandles.eventindex = eventindex;
fighandles.trainingindex = trainingindex;
fighandles.trainingfilter = trainingfilter;
fighandles.decodefilter = decodefilter;
fighandles.decodedata = decodedata;
fighandles.matches = matches;
fighandles.locator = locator;
fighandles.maxheight = maxheight;
fighandles.homesegmentlength = homesegmentlength;
fighandles.frameperiod = frameperiod;
fighandles.makeseries = 0;
fighandles.direction = 0;

guidata(fighandle,fighandles);
%-------------------------------------------------------------
function downtrainingbutton_Callback(hObject,fighandles)
% callback for the decrease speed button
stop(fighandles.plottimer);
currenttraining = fighandles.trainingindex;
newtrain = max([1 currenttraining-1]);
fighandles.trainingindex = newtrain;
set(fighandles.trainingdisplay,'String',num2str(newtrain));
guidata(fighandles.fighandle,fighandles);
loadevent(fighandles);
%-----------------------------------------------------------
function uptrainingbutton_Callback(hObject,fighandles)
% callback for the decrease speed button
stop(fighandles.plottimer);
currenttraining = fighandles.trainingindex;
newtrain = currenttraining+1;
fighandles.trainingindex = newtrain;
set(fighandles.trainingdisplay,'String',num2str(newtrain));
guidata(fighandles.fighandle,fighandles);
loadevent(fighandles);
%------------------------------------------------------
function downspeedbutton_Callback(hObject,fighandles)
% callback for the decrease speed button
currentspeed = fighandles.frameperiod;
newspeed = max([.01 currentspeed-.01]);
fighandles.frameperiod = newspeed;
set(fighandles.speeddisplay,'String',num2str(newspeed));
set(fighandles.plottimer,'Period',fighandles.frameperiod);
guidata(fighandles.fighandle,fighandles);
%--------------------------------------------------------------
function upspeedbutton_Callback(hObject,fighandles)
%callback for the increase speed button
currentspeed = fighandles.frameperiod;
newspeed = min([1 currentspeed+.01]);
fighandles.frameperiod = newspeed;
set(fighandles.speeddisplay,'String',num2str(newspeed));
set(fighandles.plottimer,'Period',fighandles.frameperiod);
guidata(fighandles.fighandle,fighandles);
%-------------------------------------------------------------
function downeventbutton_Callback(hObject,fighandles)
% callback for the decrease speed button

stop(fighandles.plottimer);
currentevent = fighandles.eventindex;
newevent = max([1 currentevent-1]);
fighandles.eventindex = newevent;
set(fighandles.eventdisplay,'String',num2str(newevent));
guidata(fighandles.fighandle,fighandles);
loadevent(fighandles);
%--------------------------------------------------------------
function upeventbutton_Callback(hObject,fighandles)
%callback for the increase speed button

stop(fighandles.plottimer);
numevents = size(fighandles.decodefilter(fighandles.index(1)).output{1}(fighandles.index(2)).eventtime,1);
currentevent = fighandles.eventindex;
newevent = min([numevents currentevent+1]);
fighandles.eventindex = newevent;
set(fighandles.eventdisplay,'String',num2str(newevent));
guidata(fighandles.fighandle,fighandles);
loadevent(fighandles);
%--------------------------------------------------------
function list1_Callback(hObject,fighandles)

newvalue = get(hObject,'Value');
if (newvalue > 1)
    cellnum = fighandles.eventcellsactive(newvalue-1);
    trainingcell = find(fighandles.matches == cellnum);
    placefieldplot(trainingcell,fighandles);
    
else
    for plotnum = 1:3
        set(fighandles.placeplot(plotnum),'Visible','off');
    end
end
%---------------------------------------------------------
function makefigurebutton_Callback(hObject,fighandles)

newfig = figure;
copyobj(fighandles.a(4),newfig);
%----------------------------------------------------------
function makeseriesbutton_Callback(hObject,fighandles)

fighandles.makeseries = 1;
playbutton_Callback(hObject,fighandles);

%------------------------------------------------------------
function playbutton_Callback(hObject,fighandles)
%callback for the play button

fighandles.direction = 1;
r = get(fighandles.plottimer,'Running');
if (strcmp(r,'off'))
    start(fighandles.plottimer);
end
guidata(fighandles.fighandle,fighandles);
%------------------------------------------------------------
function viewdecodeevent2_CloseRequestFnc(hObject,fighandles)
%this function is called when the window is closed

r = get(fighandles.plottimer,'Running');
fighandle = fighandles.fighandle;
if (strcmp(r,'on'))
    stop(fighandles.plottimer);
end
delete(timerfind('Tag',['plottimer',num2str(fighandle)]));
delete(hObject);
%----------------------------------------------------------
function loadevent(fighandles)

animalnum = fighandles.index(1);
epochnum = fighandles.index(2);
trajmapping = fighandles.trajmapping;
binsize = fighandles.binsize;
eventindex = fighandles.eventindex;
trainingindex = fighandles.trainingindex;
a = fighandles.a;
p = fighandles.p;
matches = fighandles.matches;

trainingdata = [];
spikedata = [];
decodedata = [];
indexlist = [];
activeindexlist = [];
activespiketimes = [];
activerates = [];

%pick out all the matching cells from the training data and the
%decoding data
%traindata contains linear rates, and is n by x, where n is the
%number of cells and x is the number of spatial bins
%spikedata contains spikecounts, and is n by t, where t is the
%number of temporal bins in the data to be decoded.

matches = rowfind(fighandles.trainingfilter(animalnum).output{trainingindex}(epochnum).index(:,[1 3 4]),fighandles.decodefilter(animalnum).output{1}(epochnum).index(:,[1 3 4])); %find the matching cell indices


% matches1 = rowfind(fighandles.trainingfilter(animalnum).output{1}(epochnum).index(:,[1 3 4]),fighandles.decodefilter(animalnum).output{1}(epochnum).index(:,[1 3 4])); %find the matching cell indices
% matches2 = rowfind(fighandles.trainingfilter(animalnum).output{2}(epochnum).index(:,[1 3 4]),fighandles.decodefilter(animalnum).output{1}(epochnum).index(:,[1 3 4])); %find the matching cell indices
% 
% if (trainingindex == 1)
%     matches = setdiff(matches1,matches2);
% else
%     matches = setdiff(matches2,matches1);
% end

fighandles.matches = matches;
startevent = fighandles.decodefilter(animalnum).output{1}(epochnum).eventtime(eventindex,1);
endevent = fighandles.decodefilter(animalnum).output{1}(epochnum).eventtime(eventindex,2);
if (endevent-startevent) > 2
    disp('Event too long');
    return
end
fighandles.timebins = startevent:binsize:endevent;
fighandles.currentbin = 1;
fighandles.eventcellsactive = [];
fighandles.eventcellpeaks = [];
activecount = 0;
for trainingcell = 1:length(matches)
    if (matches(trainingcell) > 0) %we have a match
        indexlist = [indexlist; fighandles.trainingfilter(animalnum).output{trainingindex}(epochnum).index(trainingcell,:)];
        trainingdata = [trainingdata; fighandles.trainingfilter(animalnum).output{trainingindex}(epochnum).rates(trainingcell,:)];
        [tmppeak, tmppeakind] = max(fighandles.trainingfilter(animalnum).output{trainingindex}(epochnum).rates(trainingcell,:));
        tmppeakdist = fighandles.trainingfilter(animalnum).output{trainingindex}(epochnum).dist(tmppeakind);
        tmpspiketimes = fighandles.decodefilter(animalnum).output{1}(epochnum).eventdata(eventindex).spiketimes(find(fighandles.decodefilter(animalnum).output{1}(epochnum).eventdata(eventindex).cellindex == matches(trainingcell)));
        %save all the info for the active cells
        if ~isempty(tmpspiketimes)
            activecount = activecount+1;
            activespiketimes{activecount} = tmpspiketimes;
            activerates = [activerates; fighandles.trainingfilter(animalnum).output{trainingindex}(epochnum).rates(trainingcell,:)];
            fighandles.eventcellsactive = [fighandles.eventcellsactive matches(trainingcell)];
            activeindexlist = [activeindexlist; fighandles.trainingfilter(animalnum).output{trainingindex}(epochnum).index(trainingcell,:)];
            fighandles.eventcellpeaks = [fighandles.eventcellpeaks tmppeakdist];
        end
        spikebins = lookup(tmpspiketimes,fighandles.timebins);
        spikecount = zeros(1,length(fighandles.timebins));
        for i = 1:length(spikebins)
            spikecount(spikebins(i)) = spikecount(spikebins(i))+1;
        end
        spikedata = [spikedata; spikecount];
    end
end
trainingdata = trainingdata*binsize; %transform rates to expected number of spikes
activerates = activerates*binsize;

if ~isempty(fighandles.eventcellsactive)
    eventcells = [fighandles.eventcellsactive' fighandles.eventcellpeaks'];
    [eventcells, sorti] = sortrows(eventcells,2);
    fighandles.eventcellsactive = eventcells(:,1);
    fighandles.eventcellpeaks = eventcells(:,2);
    fighandles.activeindexlist = activeindexlist(sorti,:);
else
    eventcells = [];    
end
fighandles.displayeventcells = {'None'};
tickplot = [];
for i = 1:length(fighandles.eventcellsactive)
    fighandles.displayeventcells{i+1} = ['Cell ',num2str(fighandles.eventcellsactive(i))];
    tmpspiketimes = fighandles.decodefilter(animalnum).output{1}(epochnum).eventdata(eventindex).spiketimes(find(fighandles.decodefilter(animalnum).output{1}(epochnum).eventdata(eventindex).cellindex == eventcells(i)));
    tmpspiketimes(:,2) = i;
    tickplot = [tickplot; tmpspiketimes];
end
set(fighandles.rasterplot,'XData',tickplot(:,1));
set(fighandles.rasterplot,'YData',tickplot(:,2));
ptmp = get(fighandles.rasterplot, 'Parent');
%set(ptmp, 'XLim', [min(tickplot(:,1))-.1 max(tickplot(:,1))+.1]);
etime = fighandles.decodefilter(animalnum).output{1}(epochnum).eventtime(eventindex,:);
set(ptmp, 'XLim', [etime(1) - .1 etime(2) + .1]);
set(fighandles.listbox1, 'String',fighandles.displayeventcells);




%the decoded data contains spatial probabilities, and is x by t
decodedata = zeros(size(trainingdata,2),size(spikedata,2));
naninds = find(isnan(trainingdata(1,:)));
for t = 1:size(spikedata,2) %time  
    Tspikecount = repmat(spikedata(:,t),1,size(trainingdata,2));
    %calculate P(spikecount|x) for this timebin across all cells and all x 
    spatialprob = prod(((trainingdata.^Tspikecount)./factorial(Tspikecount)).*exp(-trainingdata),1)'; 
    spatialprob(naninds) = 0;
    if (sum(spikedata(:,t)) ~= 0)
	spatialprob = spatialprob/sum(spatialprob);  %normalize across space
    else
	spatialprob(:) = 0;
    end
    decodedata(:,t) = spatialprob;  
end

% if (length(activespiketimes) > 4)
%     stats = calcReplayStats(activespiketimes,activerates,fighandles.timebins,fighandles.trainingfilter(animalnum).output{trainingindex}(epochnum).dist);
% else
%     stats = nan;
% end
% disp(stats)
xdata = {[],[]};
ydata = {[],[]};
plotydata = {[],[],[]};
plotxdata = {[],[],[]};


%combine the outbound and inbound trajectories
for traj = [1:4];    
    trajindex = find(fighandles.trainingfilter(animalnum).output{trainingindex}(epochnum).traj == traj);
    xdata{trajmapping(traj)} = fighandles.trainingfilter(animalnum).output{trainingindex}(epochnum).dist(trajindex);
    ydata{trajmapping(traj)} = [ydata{trajmapping(traj)}; decodedata(trajindex,fighandles.currentbin)'];    
end

%combine data for the home arm
for traj = 1:2
    ydata{traj} = sum(ydata{traj});
    homeindex = find(xdata{traj} < fighandles.homesegmentlength);
    outerindex = find(xdata{traj} >= fighandles.homesegmentlength);
    plotxdata{1} = xdata{traj}(homeindex);
    plotxdata{traj+1} = xdata{traj}(outerindex);
    plotydata{1} = [plotydata{1}; ydata{traj}(homeindex)];
    plotydata{traj+1} = ydata{traj}(outerindex);
end
plotydata{1} = sum(plotydata{1});

%plot the three plots
for plotnum = 1:3
    
    set(p(plotnum),'YData',plotydata{plotnum});
    set(p(plotnum),'XData',plotxdata{plotnum});
    try
        animallocation = fighandles.decodefilter(animalnum).output{1}(epochnum).eventdist(fighandles.eventindex);
        animaltraj = fighandles.decodefilter(animalnum).output{1}(epochnum).eventtraj(fighandles.eventindex);
    catch
        animallocation = 0;
        animaltraj = 0;
    end
    if (animallocation < fighandles.homesegmentlength)
        set(fighandles.locator{1}, 'XData', [animallocation animallocation]);
        set(fighandles.locator{1},'Visible','on');
        set(fighandles.locator{2},'Visible','off');
        set(fighandles.locator{3},'Visible','off');
    else
        if (animaltraj > 0)
            if (trajmapping(animaltraj)+1 == plotnum)
                set(fighandles.locator{plotnum}, 'XData', [animallocation animallocation]);
                set(fighandles.locator{plotnum},'Visible','on');
            else
                set(fighandles.locator{plotnum},'Visible','off');
            end
        else
            set(fighandles.locator{plotnum},'Visible','off');
        end
    end
end

for plotnum = 1:3
    set(fighandles.placeplot(plotnum),'Visible','off');
end
set(fighandles.listbox1,'Value',1);

fighandles.decodedata = decodedata;
set(fighandles.timedisplay,'String',num2str(fighandles.timebins(fighandles.currentbin)));
drawnow
guidata(fighandles.fighandle,fighandles);
%------------------------------------------------------------
function plotTimerFcn(hobject)
% this is called by the plottimer object when its 'Running' property is set
% to 'on'

userdata = get(hobject,'UserData');
fighandle = userdata(1);
fighandles = guidata(fighandle);
p = fighandles.p;

if (fighandles.makeseries)
    makereplayfig(fighandles);
end


%  newfig = figure(101);
%  figuredata = [];
%  currentfiguredata = get(newfig,'UserData');
%  if isempty(currentfiguredata)
%     figuredata(1,1) = copyobj(fighandles.a(1),newfig);
%     cla(figuredata(1,1));
%     figuredata(1,2) = copyobj(fighandles.a(2),newfig);
%     cla(figuredata(1,2));
%     figuredata(1,3) = copyobj(fighandles.a(3),newfig);
%     cla(figuredata(1,3));
%     figuredata(1,4) = 1;
%  else
%      figuredata = currentfiguredata;
%      colors = jet(size(fighandles.decodedata,2));
%      tmpplot1 = copyobj(p(1),figuredata(1));
%      tmpplot2 = copyobj(p(2),figuredata(2));
%      tmpplot3 = copyobj(p(3),figuredata(3));
%      figuredata(4) = figuredata(4)+1;
%      
%      set(tmpplot1,'Color',colors(figuredata(4)-1,:));
%      set(tmpplot1,'LineWidth',.5);
%      set(tmpplot2,'Color',colors(figuredata(4)-1,:));
%      set(tmpplot2,'LineWidth',.5);
%      set(tmpplot3,'Color',colors(figuredata(4)-1,:));
%      set(tmpplot3,'LineWidth',.5);
%  end
%  set(newfig,'UserData',figuredata);
%  
    

animalnum = fighandles.index(1);
epochnum = fighandles.index(2);
trajmapping = fighandles.trajmapping;
trainingindex = fighandles.trainingindex;



if (fighandles.currentbin < size(fighandles.decodedata,2))
    fighandles.currentbin = fighandles.currentbin+1;
else
    fighandles.currentbin = 1;
    stop(fighandles.plottimer);
    if (fighandles.makeseries)
        fighandles.makeseries = 0;
    end
end

%timelength = fighandles.timelength;

% calculate the time elapsed since the last timer call
% and the best time point within the data
instantperiod = get(hobject,'InstantPeriod');


xdata = {[],[]};
ydata = {[],[]};
plotydata = {[],[],[]};
plotxdata = {[],[],[]};

%combine the outbound and inbound trajectories
for traj = [1:4];    
    trajindex = find(fighandles.trainingfilter(animalnum).output{trainingindex}(epochnum).traj == traj);
    xdata{trajmapping(traj)} = fighandles.trainingfilter(animalnum).output{trainingindex}(epochnum).dist(trajindex);
    ydata{trajmapping(traj)} = [ydata{trajmapping(traj)}; fighandles.decodedata(trajindex,fighandles.currentbin)'];    
end

%combine data for the home arm
for traj = 1:2
    ydata{traj} = sum(ydata{traj});
    homeindex = find(xdata{traj} < fighandles.homesegmentlength);
    outerindex = find(xdata{traj} >= fighandles.homesegmentlength);
    plotxdata{1} = [xdata{traj}(homeindex)];
    plotxdata{traj+1} = xdata{traj}(outerindex);
    plotydata{1} = [plotydata{1}; ydata{traj}(homeindex)];
    plotydata{traj+1} = ydata{traj}(outerindex);
end
plotydata{1} = sum(plotydata{1});

%plot the three plots
for plotnum = 1:3
    set(p(plotnum),'XData',plotxdata{plotnum});
    set(p(plotnum),'YData',plotydata{plotnum});
    try
        animallocation = fighandles.decodefilter(animalnum).output{1}(epochnum).eventdist(fighandles.eventindex);
        animaltraj = fighandles.decodefilter(animalnum).output{1}(epochnum).eventtraj(fighandles.eventindex);
    catch
        animallocation = 0;
        animaltraj = 0;
    end
    if (animallocation < fighandles.homesegmentlength)
        set(fighandles.locator{1}, 'XData', [animallocation animallocation]);
        set(fighandles.locator{1},'Visible','on');
        set(fighandles.locator{2},'Visible','off');
        set(fighandles.locator{3},'Visible','off');
    else
        if (animaltraj > 0)
            if (trajmapping(animaltraj)+1 == plotnum)
                set(fighandles.locator{plotnum}, 'XData', [animallocation animallocation]);
                set(fighandles.locator{plotnum},'Visible','on');
            else
                set(fighandles.locator{plotnum},'Visible','off');
            end
        else
            set(fighandles.locator{plotnum},'Visible','off');
        end
    end
end

set(fighandles.timedisplay,'String',num2str(fighandles.timebins(fighandles.currentbin)));
drawnow;

%fighandles.timeindex = timeindex;
guidata(fighandles.fighandle,fighandles);
%-------------------------------------------------------------
function makereplayfig(fighandles)


p = fighandles.p;


 newfig = figure(101);
 figuredata = [];
 currentfiguredata = get(newfig,'UserData');
 if isempty(currentfiguredata)


     figuredata(1,1) = copyobj(fighandles.a(1),newfig);
     cla(figuredata(1,1));
     figuredata(1,2) = copyobj(fighandles.a(2),newfig);
     cla(figuredata(1,2));
     figuredata(1,3) = copyobj(fighandles.a(3),newfig);
     cla(figuredata(1,3));
     figuredata(1,4) = 1;
     tmpdata = copyobj(fighandles.a(4),newfig);
%    figuredata(1,5) = copyobj(fighandles.a(4),newfig);
     spikeaxis = axis(tmpdata);

     % look up the cell indeces for the cells active in this event
     animalnum = fighandles.index(1);
     epochnum = fighandles.index(2);
     adir = fighandles.decodefilter(animalnum).animal{2};
     aname = fighandles.decodefilter(animalnum).animal{3};
     eventcellsactive = fighandles.eventcellsactive;

     % HACK: Reorganize the cell list for Frank event #37
     %eventcellsactive = [1 38 24 54 61 70 66 48 49 56 4];

%     clist = fighandles.decodefilter(animalnum).output{1}(epochnum).index(eventcellsactive,:)
    % for Bond :
    clist = [ 4     3    11     4 ;
	      4     3    12     1 ;
	      4     3    14     5 ;
	      4     3    13     1 ;
	      4     3    14     3 ;
	      4     3    10     1 ;
	      4     3    12     3 ;
	      4     3     5     1 ;
	      4     3    18     1 ;
	      4     3     1     1 ;
	      4     3     5     2 ;
	      4     3    19     2 ;
	      4     3     1    10 ;
	      4     3    17     1 ;
	      4     3    19     1 ;
	      4     3    29     4 ;
	      4     3     2     4 ;
	      4     3    11     5 ]


    epoch = fighandles.decodefilter(fighandles.index(1)).epochs{1}(fighandles.index(2),:);

    if (epoch(1) < 10)
         daystr = ['0',num2str(epoch(1))];
     else
         daystr = num2str(epoch(1));
     end

     figuredata(1,5) = axes('position',[.07 .07 .4 .2]);
     datarange = [spikeaxis(1) spikeaxis(2)];
     % plot the raster of spikes
     plotrasters(adir, aname, clist, spikeaxis(1:2));
     set(gca, 'YLim', [0 18]);

     figuredata(1,7) = axes('position',[.07 .27 .4 .15]);
     loadfile = [adir,'EEG/',aname,'ripple',daystr,'-',num2str(epoch(2)),'-05.mat'];
%     loadfile = [adir,'EEG/',aname,'ripple',daystr,'-',num2str(epoch(2)),'-01.mat'];
     load(loadfile);
     ripple = ripple{epoch(1)}{epoch(2)}{5};
%     ripple = ripple{epoch(1)}{epoch(2)}{1};
     times = geteegtimes(ripple);
     goodind = find((times >= spikeaxis(1))&(times <= spikeaxis(2)));
     plot(times(goodind),ripple.data(goodind,1), 'k');
     axis([spikeaxis(1) spikeaxis(2) -150 150]);
     axis off;
     hold on;


     %axis([spikeaxis(1) spikeaxis(2) -150 150]);

     loadfile = [adir, aname,'pos',daystr,'.mat'];
     load(loadfile);
     pos = pos{epoch(1)}{epoch(2)};
     goodind = find((pos.data(:,1) >= spikeaxis(1)-5)&(pos.data(:,1) <= spikeaxis(2)+5));

     figuredata(1,9) = axes('position',[.1 .85 .15 .15]);
     plot(pos.data(:,2),pos.data(:,3),'.', 'Color', [.5 .5 .5], 'MarkerSize', 3);
     hold on;
     %plot(pos.data(goodind,2),pos.data(goodind,3),'y.');
     plot(pos.data(goodind(1:round(length(goodind)/2)),2),pos.data(goodind(1:round(length(goodind)/2)),3),'r.','MarkerSize',6);
     plot(pos.data(goodind(round(length(goodind)/2:end)),2),pos.data(goodind(round(length(goodind)/2):end),3),'g.','MarkerSize',6);
     plot(pos.data(goodind(round(length(goodind)/2)),2),pos.data(goodind(round(length(goodind)/2)),3),'k.','MarkerSize',6);

     axis([min(pos.data(find(~isnan(pos.data(:,2))),2)) max(pos.data(find(~isnan(pos.data(:,2))),2)) min(pos.data(find(~isnan(pos.data(:,2))),3)) max(pos.data(find(~isnan(pos.data(:,2))),3))]);

      figure
      orient tall
      geom = [8 ceil(size(clist,1) / 2)];
      pind = [fighandles.index(1) fighandles.trainingindex fighandles.index(2)];
      plotfilterfields2d(fighandles.trainingfilter, pind, clist, geom);


 %    tmptickplot = [];
 %    load('/data/mkarlsso/Bon/bonspikes04');
 %    load('/data/mkarlsso/Bon/bonlinpos04');
 %    for i = 1:length(fighandles.eventcellsactive)
 %        tmpspiketimes = spikes{4}{2}{fighandles.activeindexlist(i,3)}{fighandles.activeindexlist(i,4)}.data(:,1);
 %        tmpspiketimes(:,2) = i;
 %        tmptickplot = [tmptickplot; tmpspiketimes];
 %    end
 %    plot(tmptickplot(:,1),tmptickplot(:,2),'k+');
 %    hold on
 %    plot(linpos{4}{2}.statematrix.time,16*(linpos{4}{2}.statematrix.lindist/150));
 else
     figuredata = currentfiguredata;
     nonzero = length(find(sum(fighandles.decodedata)));
     %colors = jet(size(fighandles.decodedata,2));
     colors = jet(nonzero);
     %colors = [repmat([0 0 0],7,1);colors];
     tmpplot1 = copyobj(p(1),figuredata(1));
     tmpplot2 = copyobj(p(2),figuredata(2));
     tmpplot3 = copyobj(p(3),figuredata(3));
     if (sum(get(tmpplot1, 'YData')) | ...
         sum(get(tmpplot2, 'YData')) | ...
         sum(get(tmpplot3, 'YData')))
	 figuredata(4) = figuredata(4)+1;
	 set(tmpplot1,'Color',colors(figuredata(4)-1,:));
	 set(tmpplot1,'LineWidth',.5);
	 set(tmpplot2,'Color',colors(figuredata(4)-1,:));
	 set(tmpplot2,'LineWidth',.5);
	 set(tmpplot3,'Color',colors(figuredata(4)-1,:));
	 set(tmpplot3,'LineWidth',.5);
	 % add a patch to the ripple plot
     	 axes(figuredata(1,7));
	 x = fighandles.timebins(fighandles.currentbin);
	 y = 100;
	 patch([x-0.015 x-0.015 x x], [y y+30 y+30 y], ...
		 colors(figuredata(4)-1, :), ...
		 'EdgeColor', colors(figuredata(4)-1, :));
     else
	 set(tmpplot1,'LineWidth',0.001);
	 set(tmpplot2,'LineWidth',0.001);
	 set(tmpplot3,'LineWidth',0.001);
     end

     %set(tmpplot1,'Color',colors(figuredata(4)-1,:));
     %set(tmpplot1,'LineWidth',.5);
     %set(tmpplot2,'Color',colors(figuredata(4)-1,:));
     %set(tmpplot2,'LineWidth',.5);
     %set(tmpplot3,'Color',colors(figuredata(4)-1,:));
     %set(tmpplot3,'LineWidth',.5);
 end
 set(newfig,'UserData',figuredata);
 
 
%------------------------------------------------------------
function placefieldplot(trainingcell,fighandles)

animalnum = fighandles.index(1);
epochnum = fighandles.index(2);
trajmapping = fighandles.trajmapping;
binsize = fighandles.binsize;
eventindex = fighandles.eventindex;
trainingindex = fighandles.trainingindex;
a = fighandles.a;
p = fighandles.p;
matches = fighandles.matches;

cellrates = fighandles.trainingfilter(animalnum).output{trainingindex}(epochnum).rates(trainingcell,:);
maxcellrate = max(cellrates);
cellrates = (cellrates/maxcellrate)*fighandles.maxheight;

xdata = {[],[]};
ydata = {[],[]};
plotydata = {[],[],[]};
plotxdata = {[],[],[]};

%combine the outbound and inbound trajectories
for traj = [1:4];    
    trajindex = find(fighandles.trainingfilter(animalnum).output{trainingindex}(epochnum).traj == traj);
    xdata{trajmapping(traj)} = fighandles.trainingfilter(animalnum).output{trainingindex}(epochnum).dist(trajindex);
    ydata{trajmapping(traj)} = [ydata{trajmapping(traj)}; cellrates(trajindex)];    
end

%combine data for the home arm
for traj = 1:2
    ydata{traj} = max(ydata{traj});
    homeindex = find(xdata{traj} < fighandles.homesegmentlength);
    outerindex = find(xdata{traj} >= fighandles.homesegmentlength);
    plotxdata{1} = [xdata{traj}(homeindex)];
    plotxdata{traj+1} = xdata{traj}(outerindex);
    plotydata{1} = [plotydata{1}; ydata{traj}(homeindex)];
    plotydata{traj+1} = ydata{traj}(outerindex);
end
plotydata{1} = max(plotydata{1});

%plot the three plots
for plotnum = 1:3
    set(fighandles.placeplot(plotnum),'XData',plotxdata{plotnum});
    set(fighandles.placeplot(plotnum),'YData',plotydata{plotnum});
    set(fighandles.placeplot(plotnum),'Visible','on');
end





%---------------------------------------------------------------
