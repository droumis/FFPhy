function h = viewdecodeevent(varargin)
%H = viewdecodeevent(index, trainingfilter,decodefilter)
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
        'Tag','viewdecodeevent','NumberTitle','off','Resize','on','Name',['Place cell movie viewer: cell ',num2str(cellindex)], ... 
        'CloseRequestFcn','viewdecodeevent(''viewdecodeevent_CloseRequestFnc'',gcbo,guidata(gcbo));', ...
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
trainingindex = 2;
animalnum = index(1);
epochnum = index(2);
eventtime = decodefilter(animalnum).output{1}(epochnum).eventtime(eventindex);
homesegmentlength = trainingfilter(animalnum).output{1}(epochnum).homesegmentlength;
maxheight = .4;

figpos = get(fighandle,'Position');   

%create the uicontrols
b(1) = uicontrol('Tag','playbutton','Style','pushbutton','Units','pixel','Position',[5 105 50 20],'String','Play','Callback','viewdecodeevent(''playbutton_Callback'',gcbo,guidata(gcbo));');
b(2) = uicontrol('Tag','downspeedbutton','Style','pushbutton','Units','pixel','Position',[165 105 20 20],'String','<','Callback','viewdecodeevent(''downspeedbutton_Callback'',gcbo,guidata(gcbo));');
b(3) = uicontrol('Tag','upspeedbutton','Style','pushbutton','Units','pixel','Position',[185 105 20 20],'String','>','Callback','viewdecodeevent(''upspeedbutton_Callback'',gcbo,guidata(gcbo));');
b(4) = uicontrol('Tag','downeventbutton','Style','pushbutton','Units','pixel','Position',[165 70 20 20],'String','<','Callback','viewdecodeevent(''downeventbutton_Callback'',gcbo,guidata(gcbo));');
b(5) = uicontrol('Tag','upeventbutton','Style','pushbutton','Units','pixel','Position',[185 70 20 20],'String','>','Callback','viewdecodeevent(''upeventbutton_Callback'',gcbo,guidata(gcbo));');
b(6) = uicontrol('Tag','downtrainingbutton','Style','pushbutton','Units','pixel','Position',[55 70 20 20],'String','<','Callback','viewdecodeevent(''downtrainingbutton_Callback'',gcbo,guidata(gcbo));');
b(7) = uicontrol('Tag','uptrainingbutton','Style','pushbutton','Units','pixel','Position',[75 70 20 20],'String','>','Callback','viewdecodeevent(''uptrainingbutton_Callback'',gcbo,guidata(gcbo));');
b(8) = uicontrol('Tag','makefigurebutton','Style','pushbutton','Units','pixel','Position',[5 40 50 20],'String','Figure','Callback','viewdecodeevent(''makefigurebutton_Callback'',gcbo,guidata(gcbo));');


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
'Position',[210 5 70 120],'String',displayeventcells,'Max',1,'Min',1,'Callback','viewdecodeevent(''list1_Callback'',gcbo,guidata(gcbo));', 'Value',[1]);

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
rasterplot = plot(tickplot(:,1),tickplot(:,2),'+');
hold on;
set(a(4),'YDir','reverse');
set(l(1), 'String',displayeventcells);

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
    animallocation = decodefilter(animalnum).output{1}(epochnum).eventdist(eventindex);
    animaltraj = decodefilter(animalnum).output{1}(epochnum).eventtraj(eventindex);
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
plottimer = timer('Tag',['plottimer',num2str(fighandle)],'Period',frameperiod,'TimerFcn',['viewdecodeevent(''plotTimerFcn'',timerfind(''Tag'',''plottimer',num2str(fighandle),'''));'],'TasksToExecute',1000000,'ExecutionMode','fixedDelay','UserData',[fighandle 1]);

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
function viewdecodeevent_CloseRequestFnc(hObject,fighandles)
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
activespiketimes = [];
activerates = [];

%pick out all the matching cells from the training data and the
%decoding data
%traindata contains linear rates, and is n by x, where n is the
%number of cells and x is the number of spatial bins
%spikedata contains spikecounts, and is n by t, where t is the
%number of temporal bins in the data to be decoded.
matches = rowfind(fighandles.trainingfilter(animalnum).output{trainingindex}(epochnum).index(:,[1 3 4]),fighandles.decodefilter(animalnum).output{1}(epochnum).index(:,[1 3 4])); %find the matching cell indices
fighandles.matches = matches;
startevent = fighandles.decodefilter(animalnum).output{1}(epochnum).eventtime(eventindex,1);
endevent = fighandles.decodefilter(animalnum).output{1}(epochnum).eventtime(eventindex,2);
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
    eventcells = sortrows(eventcells,2);
    fighandles.eventcellsactive = eventcells(:,1);
    fighandles.eventcellpeaks = eventcells(:,2);
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
set(fighandles.listbox1, 'String',fighandles.displayeventcells);

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
    animallocation = fighandles.decodefilter(animalnum).output{1}(epochnum).eventdist(fighandles.eventindex);
    animaltraj = fighandles.decodefilter(animalnum).output{1}(epochnum).eventtraj(fighandles.eventindex);
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

% newfig = figure;
% copyobj(fighandles.a(1),newfig);
% copyobj(fighandles.a(2),newfig);
% copyobj(fighandles.a(3),newfig);

animalnum = fighandles.index(1);
epochnum = fighandles.index(2);
trajmapping = fighandles.trajmapping;
trainingindex = fighandles.trainingindex;

p = fighandles.p;

if (fighandles.currentbin < size(fighandles.decodedata,2))
    fighandles.currentbin = fighandles.currentbin+1;
else
    fighandles.currentbin = 1;
    stop(fighandles.plottimer);
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
    animallocation = fighandles.decodefilter(animalnum).output{1}(epochnum).eventdist(fighandles.eventindex);
    animaltraj = fighandles.decodefilter(animalnum).output{1}(epochnum).eventtraj(fighandles.eventindex);
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
