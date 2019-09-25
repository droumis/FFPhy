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
binsize = .02; %default temporal bin
currentbin = 1;
eventindex = 1;
animalnum = index(1);
epochnum = index(2);
eventtime = decodefilter(animalnum).output{1}(epochnum).eventtime(eventindex);
maxheight = .4;




b = get(fighandle,'Position');   

%create the uicontrols
b(1) = uicontrol('Tag','playbutton','Style','pushbutton','Units','pixel','Position',[5 40 50 20],'String','Play','Callback','viewdecodeevent(''playbutton_Callback'',gcbo,guidata(gcbo));');
b(2) = uicontrol('Tag','downspeedbutton','Style','pushbutton','Units','pixel','Position',[165 40 20 20],'String','<','Callback','viewdecodeevent(''downspeedbutton_Callback'',gcbo,guidata(gcbo));');
b(3) = uicontrol('Tag','upspeedbutton','Style','pushbutton','Units','pixel','Position',[185 40 20 20],'String','>','Callback','viewdecodeevent(''upspeedbutton_Callback'',gcbo,guidata(gcbo));');
b(4) = uicontrol('Tag','downeventbutton','Style','pushbutton','Units','pixel','Position',[165 5 20 20],'String','<','Callback','viewdecodeevent(''downeventbutton_Callback'',gcbo,guidata(gcbo));');
b(5) = uicontrol('Tag','upeventbutton','Style','pushbutton','Units','pixel','Position',[185 5 20 20],'String','>','Callback','viewdecodeevent(''upeventbutton_Callback'',gcbo,guidata(gcbo));');



e(1) = uicontrol('Tag','timedisplay','Style','edit','Units','pixel','Position',[55 40 60 20],'String',num2str(eventtime),'HorizontalAlignment','left');
e(2) = uicontrol('Tag','speeddisplay','Style','edit','Units','pixel','Position',[115 40 50 20],'String',num2str(frameperiod),'HorizontalAlignment','left');
e(3) = uicontrol('Tag','eventdisplay','Style','edit','Units','pixel','Position',[115 5 50 20],'String',num2str(eventindex),'HorizontalAlignment','left');

t(1) = uicontrol('Tag','timetext','Style','text','Units','pixel','Position',[55 60 50 12],'String','Time');
t(2) = uicontrol('Tag','speedtext','Style','text','Units','pixel','Position',[115 60 80 12],'String','Frame period');
t(3) = uicontrol('Tag','eventtext','Style','text','Units','pixel','Position',[115 25 80 12],'String','Event');

a(1) = axes('Tag','locationaxes1','Position',[.001 .74 .99 .245]);
a(2) = axes('Tag','locationaxes2','Position',[.001 .4 .99 .245]);

 
numplots = 2;



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

startevent = decodefilter(animalnum).output{1}(epochnum).eventtime(eventindex)-.5;
timebins = startevent:binsize:(startevent+1);
for trainingcell = 1:length(matches)
    if (matches(trainingcell) > 0) %we have a match
        indexlist = [indexlist; trainingfilter(animalnum).output{1}(epochnum).index(trainingcell,:)];
        trainingdata = [trainingdata; trainingfilter(animalnum).output{1}(epochnum).rates(trainingcell,:)];
        
        tmpspiketimes = decodefilter(animalnum).output{1}(epochnum).eventdata(eventindex).spiketimes(find(decodefilter(animalnum).output{1}(epochnum).eventdata(eventindex).cellindex == matches(trainingcell)));
        spikebins = lookup(tmpspiketimes,timebins);
        spikecount = zeros(1,length(timebins));
        for i = 1:length(spikebins)
            spikecount(spikebins(i)) = spikecount(spikebins(i))+1;
        end
        spikedata = [spikedata; spikecount];
    end
end
trainingdata = trainingdata*binsize; %transform rates to expected number of spikes

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
    axes(a(plotnum))
    %subplot(2,1,plotnum);
    p(plotnum) = plot(xdata{plotnum}, sum(ydata{plotnum}));
    set(p(plotnum),'LineWidth',4);
    hold on
         
    locator{plotnum} = line([0 0],[0 maxheight],'Color',[1 0 0]);
    animallocation = decodefilter(animalnum).output{1}(epochnum).eventdist(eventindex);
    animaltraj = decodefilter(animalnum).output{1}(epochnum).eventtraj(eventindex);
    if (animaltraj > 0)
        if (trajmapping(animaltraj) == plotnum)
            set(locator{plotnum}, 'XData', [animallocation animallocation]);
            set(locator{plotnum},'Visible','on');
        else
            set(locator{plotnum},'Visible','off');
        end
    else
        set(locator{plotnum}, 'XData', [animallocation animallocation]);
        set(locator{plotnum},'Visible','on');
    end
    
    axis([0 180 0 maxheight]);
    hold on
end




%creat the timer object that is called during playback
plottimer = timer('Tag',['plottimer',num2str(fighandle)],'Period',frameperiod,'TimerFcn',['viewdecodeevent(''plotTimerFcn'',timerfind(''Tag'',''plottimer',num2str(fighandle),'''));'],'TasksToExecute',1000000,'ExecutionMode','fixedDelay','UserData',[fighandle 1]);

%save important data to the figure 
fighandles = guihandles(fighandle);
fighandles.fighandle = fighandle;
fighandles.a = a;
fighandles.p = p;
fighandles.trajmapping = trajmapping;
fighandles.timebins = timebins;
fighandles.plottimer = plottimer;
fighandles.binsize = binsize;
fighandles.currentbin = currentbin;
fighandles.index = index;
fighandles.eventindex = eventindex;
fighandles.trainingfilter = trainingfilter;
fighandles.decodefilter = decodefilter;
fighandles.decodedata = decodedata;
fighandles.matches = matches;
fighandles.locator = locator;
fighandles.maxheight = maxheight;

fighandles.frameperiod = frameperiod;
fighandles.direction = 0;

%fighandles.ticplot = ticplot;

guidata(fighandle,fighandles);
%-------------------------------------------------------------
function plotTimerFcn(hobject)
% this is called by the plottimer object when its 'Running' property is set
% to 'on'

userdata = get(hobject,'UserData');
fighandle = userdata(1);
fighandles = guidata(fighandle);

animalnum = fighandles.index(1);
epochnum = fighandles.index(2);
trajmapping = fighandles.trajmapping;

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


ydata = {[],[]};
for traj = 1:4
    trajindex = find(fighandles.trainingfilter(animalnum).output{1}(epochnum).traj == traj);
    ydata{trajmapping(traj)} = [ydata{trajmapping(traj)}; fighandles.decodedata(trajindex,fighandles.currentbin)'];
    %set(p(traj),'YData',decodedata(trajindex,t));
end
for plotnum = 1:2
    set(p(plotnum),'YData',sum(ydata{plotnum}));
    animallocation = fighandles.decodefilter(animalnum).output{1}(epochnum).eventdist(fighandles.eventindex);
    animaltraj = fighandles.decodefilter(animalnum).output{1}(epochnum).eventtraj(fighandles.eventindex);
    if (animaltraj > 0)
        if (trajmapping(animaltraj) == plotnum)
            set(fighandles.locator{plotnum}, 'XData', [animallocation animallocation]);
            set(fighandles.locator{plotnum},'Visible','on');
        else
            set(fighandles.locator{plotnum},'Visible','off');
        end
    else
        set(fighandles.locator{plotnum}, 'XData', [animallocation animallocation]);
        set(fighandles.locator{plotnum},'Visible','on');
    end
end
set(fighandles.timedisplay,'String',num2str(fighandles.timebins(fighandles.currentbin)));
drawnow






%update the screen to the current time
%timeplot(timeindex,fighandles);


%fighandles.timeindex = timeindex;
guidata(fighandles.fighandle,fighandles);

%-------------------------------------------------------------
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
numevents = length(fighandles.decodefilter(fighandles.index(1)).output{1}(fighandles.index(2)).eventtime);
currentevent = fighandles.eventindex;
newevent = min([numevents currentevent+1]);
fighandles.eventindex = newevent;
set(fighandles.eventdisplay,'String',num2str(newevent));
guidata(fighandles.fighandle,fighandles);
loadevent(fighandles);
%------------------------------------------------------------
function loadevent(fighandles)


animalnum = fighandles.index(1);
epochnum = fighandles.index(2);
trajmapping = fighandles.trajmapping;
binsize = fighandles.binsize;
eventindex = fighandles.eventindex;
a = fighandles.a;
p = fighandles.p;
matches = fighandles.matches;

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

startevent = fighandles.decodefilter(animalnum).output{1}(epochnum).eventtime(eventindex)-.5;
fighandles.timebins = startevent:binsize:(startevent+1);
fighandles.currentbin = 1;
for trainingcell = 1:length(matches)
    if (matches(trainingcell) > 0) %we have a match
        indexlist = [indexlist; fighandles.trainingfilter(animalnum).output{1}(epochnum).index(trainingcell,:)];
        trainingdata = [trainingdata; fighandles.trainingfilter(animalnum).output{1}(epochnum).rates(trainingcell,:)];
        
        tmpspiketimes = fighandles.decodefilter(animalnum).output{1}(epochnum).eventdata(eventindex).spiketimes(find(fighandles.decodefilter(animalnum).output{1}(epochnum).eventdata(eventindex).cellindex == matches(trainingcell)));
        spikebins = lookup(tmpspiketimes,fighandles.timebins);
        spikecount = zeros(1,length(fighandles.timebins));
        for i = 1:length(spikebins)
            spikecount(spikebins(i)) = spikecount(spikebins(i))+1;
        end
        spikedata = [spikedata; spikecount];
    end
end
trainingdata = trainingdata*binsize; %transform rates to expected number of spikes

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

%xdata = {[],[]};
ydata = {[],[]};
for traj = [1:4];
    
    %subplot(2,1,trajmapping(traj));
    trajindex = find(fighandles.trainingfilter(animalnum).output{1}(epochnum).traj == traj);
    %xdata{trajmapping(traj)} = trainingfilter(animalnum).output{1}(epochnum).dist(trajindex);
    ydata{trajmapping(traj)} = [ydata{trajmapping(traj)}; decodedata(trajindex,1)'];
    %p(traj) = plot(trainingfilter(animalnum).output{1}(epochnum).dist(trajindex), decodedata(trajindex,1));
    %set(p(traj),'LineWidth',4);
    %hold on
    
end
for plotnum = 1:2
    set(p(plotnum),'YData',sum(ydata{plotnum}));
    animallocation = fighandles.decodefilter(animalnum).output{1}(epochnum).eventdist(fighandles.eventindex);
    animaltraj = fighandles.decodefilter(animalnum).output{1}(epochnum).eventtraj(fighandles.eventindex);
    if (animaltraj > 0)
        if (trajmapping(animaltraj) == plotnum)
            set(fighandles.locator{plotnum}, 'XData', [animallocation animallocation]);
            set(fighandles.locator{plotnum},'Visible','on');
        else
            set(fighandles.locator{plotnum},'Visible','off');
        end
    else
        set(fighandles.locator{plotnum}, 'XData', [animallocation animallocation]);
        set(fighandles.locator{plotnum},'Visible','on');
    end
end
fighandles.decodedata = decodedata;
set(fighandles.timedisplay,'String',num2str(fighandles.timebins(fighandles.currentbin)));
drawnow
guidata(fighandles.fighandle,fighandles);
%---------------------------------------------------------
function tobeginningbutton_Callback(hObject,fighandles)
%callback for the beginning button

%set the current time index for the timer object to 1
plottimer = fighandles.plottimer;
userdata = get(plottimer,'UserData');
userdata(2) = 1;
set(plottimer,'UserData',userdata);

%set the time display
currtime = timetrans(fighandles.times(1)+fighandles.basetime-fighandles.offset,1,1);
set(fighandles.timedisplay,'String',currtime);   

%direction is 0
fighandles.direction = 0;
fighandles.timeindex = 1;
%stop the timer if it is running
r = get(fighandles.plottimer,'Running');
if (strcmp(r,'on'))
    stop(fighandles.plottimer);
end
nearestindex = findnearesttime(fighandles.times(1), fighandles.basetime,fighandles.pos.data(:,1),fighandles.offset);
set(fighandles.locplot(2),'XData',fighandles.pos.data(nearestindex,2));
set(fighandles.locplot(2),'YData',fighandles.pos.data(nearestindex,3));
set(fighandles.locplot(3),'XData',fighandles.pos.data(nearestindex,4));
set(fighandles.locplot(3),'YData',fighandles.pos.data(nearestindex,5));


for i = 1:length(fighandles.ticplot)  
    set(fighandles.ticplot(i),'XData',[]);
    set(fighandles.ticplot(i),'YData',[]);   
end
guidata(fighandles.fighandle,fighandles);
timeplot(1,fighandles);
%------------------------------------------------------------
function toendbutton_Callback(hObject,fighandles)
%callback for the end button

%set the current time index for the timer object to end
plottimer = fighandles.plottimer;
userdata = get(plottimer,'UserData');
userdata(2) = fighandles.timelength;
set(plottimer,'UserData',userdata);

%set the time display
currtime = timetrans(fighandles.times(end)+fighandles.basetime-fighandles.offset,1,1);
set(fighandles.timedisplay,'String',currtime);   

%direction is 0
fighandles.direction = 0;
fighandles.timeindex = fighandles.timelength;

%stop the timer if it is running
r = get(fighandles.plottimer,'Running');
if (strcmp(r,'on'))
    stop(fighandles.plottimer);
end
nearestindex = findnearesttime(fighandles.times(end), fighandles.basetime,fighandles.pos.data(:,1),fighandles.offset);
set(fighandles.locplot(2),'XData',fighandles.pos.data(nearestindex,2));
set(fighandles.locplot(2),'YData',fighandles.pos.data(nearestindex,3));
set(fighandles.locplot(3),'XData',fighandles.pos.data(nearestindex,4));
set(fighandles.locplot(3),'YData',fighandles.pos.data(nearestindex,5));

for i = 1:length(fighandles.ticplot)  
    set(fighandles.ticplot(i),'XData',[]);
    set(fighandles.ticplot(i),'YData',[]);   
end
timeplot(fighandles.timelength,fighandles);
%------------------------------------------------------------
function playbutton_Callback(hObject,fighandles)
%callback for the play button

fighandles.direction = 1;
r = get(fighandles.plottimer,'Running');
if (strcmp(r,'off'))
    start(fighandles.plottimer);
end
guidata(fighandles.fighandle,fighandles);
%--------------------------------------------------------------
function reversebutton_Callback(hObject,fighandles)
%callback for the reverse button

fighandles.direction = -1;
r = get(fighandles.plottimer,'Running');
if (strcmp(r,'off'))
  start(fighandles.plottimer);
end
guidata(fighandles.fighandle,fighandles);
%---------------------------------------------------------------
function pausebutton_Callback(hObject,fighandles)
%callback for the pause button

fighandles.direction = 0;
r = get(fighandles.plottimer,'Running');
if (strcmp(r,'on'))
    stop(fighandles.plottimer);
end
guidata(fighandles.fighandle,fighandles);
%-------------------------------------------------------------------
function quitbutton_Callback(hObject,fighandles)
%callback for the quit button

r = get(fighandles.plottimer,'Running');
fighandle = fighandles.fighandle;
if (strcmp(r,'on'))
    stop(fighandles.plottimer);
end
delete(timerfind('Tag',['plottimer',num2str(fighandle)]));
delete(fighandles.viewdecodeevent);
%-------------------------------------------------------------
function makestill_Callback(hObject,fighandles)

pausebutton_Callback(hObject,fighandles);
newfig = figure;
for i = 1:length(fighandles.a)
    newhandle(i) = copyobj(fighandles.a(i),newfig);
end
%newaxes = axes;
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
function out = timetrans(times,UnitsPerSec,dir)
% used to translate numerical time to the ??:??:?? format

if (dir==1)
	for i = 1:length(times)
      if (times(i) > 0)
            hours = floor(times(i)/(60*60*UnitsPerSec));
	        minutes = floor(times(i)/(60*UnitsPerSec))-(hours*60);
	        seconds = floor(times(i)/(UnitsPerSec))-(hours*60*60)-(minutes*60);
      else
            hours = abs(ceil(times(i)/(60*60*UnitsPerSec)));
	        minutes = abs(ceil(times(i)/(60*UnitsPerSec))+(hours*60));
	        seconds = abs(ceil(times(i)/(UnitsPerSec))+(hours*60*60)+(minutes*60));
      end
      if (minutes<10)
          tempmin = ['0',num2str(minutes)];
      else
          tempmin = [num2str(minutes)];
      end
      if (seconds<10)
          tempseconds = ['0',num2str(seconds)];
      else
          tempseconds = [num2str(seconds)];
      end
      out{i,1} = [num2str(hours),':',tempmin,':',tempseconds];
      if (times(i)<0)
          out{i,1} = ['-',out{i,1}];
      end
    end
    
elseif (dir==2)
    for i = 1:length(times)
        t = [0 0 0 0 0 0];
        temptime = times{i};
        colons = findstr(temptime,':');
        startind = length(temptime);
        count = 6;
        for j = length(colons):-1:1
            t(count) = str2num(temptime((colons(j)+1):startind));
            startind = colons(j)-1;
            count = count-1;    
        end
        t(count) = str2num(temptime(1:colons(1)-1));
        out(i,1) = etime(t,[0 0 0 0 0 0])*UnitsPerSec;
    end
end
%-----------------------------------------------------------
function out = findnearesttime(time, basetime,postimes,offset)


time = time+basetime+offset;
[value, out] = min(abs(postimes-time));