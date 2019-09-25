function out = sj_plottrajdatafind4eeg(animdirect,prefix,index, varargin)
% sj_plottrajdatafind4eeg('/data25/sjadhav/RippleInterruption/REf_direct/','REf',[7 4]);
% Plots only EEG and Speed

% 4: Only look for eeg examples - forget spikes

% Do not use subplots. All on one figure and put linearized position on top
% of rasters. 1 eeg and ripple band above. And Speed above that


%Shantanu - CHange version 1.
% Look for only 1 kind of traj. Pick cells that you want. and change
% subplot dimensions

% Shantanu
% Load linpos data for day and epoch. Get traj start sequences and plot a
% giant figure for each trajectory / pair of trajectories in a loop so that
% you can pick example figures. You want to plot linearized position,
% linearvelocity/abs veloctiy (need pos field for that), place cell spike
% rasters (spike field), and LFP - either EEG or ripple band (eeg/ripple file),
% and DIO times marked (DIO structure - jbfill it).

% Can also load ripples structure to superimpose. Ideally, I should have a
% file which has looked for ripples combined across tetrodes

% After each figure, stop in a loop

day=index(1,1);
epoch = index(1,2);
out = 1;

Fspos = 30; %Hz
respos = 1/30; % sec
Fseeg = 1500; %Hz
reseeg = 1/1500; % sec
Fsspikes = 10000; %Hz
resspikes = 1/10000; %sec



% Use subset of eegtets for eeg and ripple: 2 each
tets = [1,5,9,10,11,12];
eegtets = [5, 9, 10, 12];
maineegtets = [9, 5];
maineegidxs = [2, 1];


% Cell indices within spikes structure - can also get by doing
% setfiltercellswith the cellinfo file

% Trajectories:
% 1 = Out Left    well 1to2
% 2 = In Left     well 2to1
% 3 = Out Right   well 1to3
% 4 = In Right    well 3to1


linposfile = sprintf('%s/%slinpos%02d.mat', animdirect, prefix, day);
load(linposfile);
statematrix = linpos{day}{epoch}.statematrix;
linpostime = statematrix.time;

posfile = sprintf('%s/%spos%02d.mat', animdirect, prefix, day);
load(posfile);
absvel = abs(pos{day}{epoch}.data(:,5)); % Can also use field 9
posn = pos{day}{epoch}.data(:,2:3); % Can also use fields 6:7
postime = pos{day}{epoch}.data(:,1); % same as linpostime

DIOfile = sprintf('%s/%sDIO%02d.mat', animdirect, prefix, day);
load(DIOfile);
% Get DIO times
stim = DIO{day}{epoch}{16};
if isempty(stim)
    stim = DIO{day}{epoch}{15};
end
stim_starttime = stim.pulsetimes(:,1)./10000; %sec
stim_endtime = stim.pulsetimes(:,2)./10000; %sec


% EEG and if desired ripple band LFP data + theta band LFP data
cnt=0;
for tet=eegtets
    cnt=cnt+1;
    eegfile = sprintf('%s/EEG/%seeg%02d-%01d-%02d.mat', animdirect, prefix, day,epoch,tet);
    load(eegfile);
    eval(['eegdata',num2str(cnt),'= eeg{day}{epoch}{tet}.data;'])
    ripplefile = sprintf('%s/EEG/%sripple%02d-%01d-%02d.mat', animdirect, prefix, day,epoch,tet);
    load(ripplefile);
    eval(['rippledata',num2str(cnt),'= ripple{day}{epoch}{tet}.data;'])
    thetafile = sprintf('%s/EEG/%stheta%02d-%01d-%02d.mat', animdirect, prefix, day,epoch,tet);
    load(thetafile);
    eval(['thetadata',num2str(cnt),'= theta{day}{epoch}{tet}.data;'])
    if cnt==1
        teeg = geteegtimes(eeg{day}{epoch}{tet});
        eegstart = eeg{day}{epoch}{tet}.starttime;
        eegsamprate = round(eeg{day}{epoch}{tet}.samprate);
        % Theta is downsampled 10 times
        ttheta = geteegtimes(theta{day}{epoch}{tet});
        thetastart = theta{day}{epoch}{tet}.starttime;
        thetasamprate = round(theta{day}{epoch}{tet}.samprate);
    end
end

%--------------------------------------------------------------------------

% Get track, wells and traj info
trajwells = linpos{day}{epoch}.trajwells;
wellCoord = linpos{day}{epoch}.wellSegmentInfo.wellCoord;

%--------------------------------------------------------------------------

% Get start index/time of each trajectory in sequence:
% 0s are skipped automatically
% You can either 1) skip -1s, or 2) include them in the previous trajectory

statematrix = linpos{day}{epoch}.statematrix;
traj = statematrix.traj;
start = find(diff(traj)~=0)+1;
start = [1; start];

% 1) skip -1s
% Get trajseq and startidxs by skipping
start(find(traj(start)==-1))=[]; % For skipping -1s
trajseq = traj(start); % Traj Identitites

%If skipped -1s, then get rid of start idxs where the ttraj just continues - unlikely
trajseq(find(diff(trajseq)==0))=[]; start(find(diff(trajseq)==0))=[];


%--------------------------------------------------------------------------

%% EEG and DIO time
%pt = stim.pulsetimes ./ 10000; % in s
%eind = lookup(pt(:,1), t);

speedthrs1 = 5; %cm/sec;



start = [start; length(traj)];
ntraj = length(start);
% Set fig prop
set(0,'defaultaxesfontsize',12);set(0,'defaultaxesfontweight','normal');
set(0,'defaultaxeslinewidth',1);
titlestr = ['Out Left';'In  Left';'OutRight';'In Right'];

% If Looking for traj 3: Out Right
% titlestr = ['OutRight'];

figno = 4;

% Plot each trajectory in a loop
for n=1:ntraj-1
    
    currstart = start(n);
    currend = start(n+1);
    currtraj = trajseq(n);
    
    if currtraj>=3,
        currwells = trajwells(2,:);
    else
        currwells = trajwells(1,:);
    end
    well1 = wellCoord(currwells(1),:); %xy position of wells
    well2 = wellCoord(currwells(2),:);
    
    % Only do if currtraj=3  
    %if currtraj == 3
        oritrajsttime = statematrix.time(currstart);
        oritrajendtime = statematrix.time(currend);
        oritimeintraj = oritrajendtime - oritrajsttime; % in sec
        trajstartidx = currstart;
        trajendidx = currend;
        
        
        trajsttime = statematrix.time(currstart)-3; % Add time to both
        trajendtime = statematrix.time(currend)+3;
        timeintraj = trajendtime - trajsttime; % in sec
        
        % Update start and end indices after adding secs to both ends
        currstart = lookup(trajsttime, statematrix.time);
        currend = lookup(trajendtime, statematrix.time);
        plotsttime = statematrix.time(currstart);
        plotendtime = statematrix.time(currend);
        
        if oritimeintraj > 2 % more than 2 secs
            
            figure(figno); hold on;
            redimscreen_widehor;
            
            % Speed
            %--------
            
            % Speed ylimits in plot: 0-5
            
            figure(figno); %subplot(nplots,1,1); hold on;
            hold on;
            currspeed = absvel(currstart:currend);
            currpostime = statematrix.time(currstart:currend);
            currpostime = currpostime - statematrix.time(currstart);
            
            % Change scale
            speedscale = max(currspeed)-min(currspeed);
            plotscale = 3-0;
            currspeed = currspeed.*(plotscale/speedscale);
            
            plot(currpostime, currspeed,'Linewidth',2);
            set(gca,'XLim',[0 max(currpostime)]);
            
            % Line for speedthrs
            yaxis = speedthrs1*(plotscale/speedscale)*1.5*ones(size([0:length(currspeed)]'));
            xaxis = [0:length(currspeed)]';
            plot(xaxis,yaxis,'r');
            speedthrscale = speedthrs1*(plotscale/speedscale);
            
            y1=min(currspeed); y2=max(currspeed);
            
            % Plot Traj Start and End points with circles
            plott = oritrajsttime-plotsttime;
            speedt = currspeed(lookup(plott, currpostime));
            plot(plott,speedt,'ro','MarkerSize',8,'LineWidth',2);
            plott = oritrajendtime-plotsttime;
            speedt = currspeed(lookup(plott, currpostime));
            plot(plott,speedt,'ro','MarkerSize',8,'LineWidth',2);
            
            % Plot 1 sec scale bar
            xaxis = 0:0.1:1;
            yaxis=(y2)*0.9*ones(size(xaxis));
            plot(xaxis,yaxis,'k-','Linewidth',4);
           % title(['Traj',num2str(currtraj),' ',titlestr, ': ',num2str(round(trajsttime)),'-',num2str(round(trajendtime))],'FontSize',14,'Fontweight','bold')
   
           
            % 2) 2 EEGs
            baseline = 4; %
            
            % Get EEG Ind in current traj
            eind1 = lookup(trajsttime, teeg);
            eind2 = lookup(trajendtime, teeg);
            
            eegtimeaxis = teeg(eind1:eind2);
            eegtimeaxis = eegtimeaxis - teeg(eind1);               
                
            tind1 = lookup(trajsttime, ttheta);
            tind2 = lookup(trajendtime, ttheta);
            thetatimeaxis = ttheta(tind1:tind2);
            thetatimeaxis = thetatimeaxis - ttheta(tind1);
                    
            figure(figno); hold on;
            for i=1:length(maineegtets)              
%                 % Get theta and ripple
%                 if i==1
%                     eval(['thetadata = thetadata',num2str(maineegidxs(i)),';']);
%                     curreegdata = double(thetadata(tind1:tind2,1));
%                     eegscale = max(curreegdata)-min(curreegdata);
%                     downeeg = baseline+1;
%                     upeeg = downeeg+5;
%                     plotscale = 5;
%                     curreegdata = downeeg + (plotscale/2) + curreegdata.*(plotscale/eegscale);
%                     plot(thetatimeaxis,curreegdata,'Linewidth',1);
%                     
%                     % Update baseline
%                     baseline = upeeg;                                     
%                     eval(['rippledata = rippledata',num2str(maineegidxs(i)),';']);
%                     curreegdata = double(rippledata(eind1:eind2,1));
%                     eegscale = max(curreegdata)-min(curreegdata);
%                     downeeg = baseline+1;
%                     upeeg = downeeg+5;
%                     plotscale = 5;
%                     curreegdata = downeeg + (plotscale/2) + curreegdata.*(plotscale/eegscale);
%                     plot(eegtimeaxis,curreegdata,'Linewidth',1);
%                     
%                     % Update baseline
%                     baseline = upeeg;
%                 end
                eval(['eegdata = eegdata',num2str(maineegidxs(i)),';']);
                curreegdata = eegdata(eind1:eind2);
                eegscale = max(curreegdata)-min(curreegdata);
                downeeg = baseline+1;
                upeeg = downeeg+5;
                plotscale = 5;
                curreegdata = downeeg + (plotscale/2) + curreegdata.*(plotscale/eegscale);              
                plot(eegtimeaxis,curreegdata,'Linewidth',1);
                
                % Update baseline
                baseline = upeeg;
                
               
                
            end
            
            % DIOs within eeg
%            currDIOtime = stim_starttime(find((stim_starttime>starttime) & (stim_starttime<endtime)));
%            y1=min(eegdata(eind1:eind2)); y2=max(eegdata(eind1:eind2));
%             for s=1:length(currDIOtime)
%                 
%                 yaxis = 0:25;
%                 DIOt = currDIOtime(s)-teeg(eind1);
%                 speedt = currspeed(lookup(DIOt, currpostime))
%                 if speedt < 3*speedthrscale
%                     %plot(DIOt*ones(size(yaxis)),yaxis,'m--','LineWidth',1.5);
%                     xaxis = DIOt:0.01:DIOt+0.1;
%                     jbfill(xaxis,25*ones(size(xaxis)),0*ones(size(xaxis)),'m','m',1,0.2);
%                 end
%                 
%             end
            
            
            
            % Ripple
            %-------
            %         figure(figno); hold on; %subplot(nplots,1,i+1); hold on;
            %         eval(['rippledata = rippledata',num2str(maineegidx),';']);
            %         curreegdata = double(rippledata(eind1:eind2,1));
            %         eegscale = max(curreegdata)-min(curreegdata);
            %         plotscale = 29-24;
            %         curreegdata =  24 + curreegdata.*(plotscale/eegscale);
            %
            %         plot(eegtimeaxis,curreegdata,'Linewidth',1);
            
            
            %set(gca,'YLim',[0 18]);
            %set(gca,'YTick',[]); set(gca,'YTickLabel',[]);
            
            
            % Return to user
            keyboard;
            close(figno);
            %        close(figno+100);
            
        end  % end timeintraj
        
    %end  % End if current traj
    
    
    
    
    
end







