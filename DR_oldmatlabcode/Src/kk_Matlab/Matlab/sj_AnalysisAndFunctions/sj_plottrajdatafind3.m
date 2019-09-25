function out = sj_plottrajdatafind3(animdirect,prefix,index, varargin)
% Same as sj_plottrajdatafind4egg
% sj_plottrajdatafind3('/data25/sjadhav/RippleInterruption/REf_direct/','REf',[7 4]);
% 1 EEG Spikes from 4 cells - ordered and speed with DIO on top

% REf - Day 7 Ep 4

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
maineegtet = 9;
maineegidx = 2;


% Cell indices within spikes structure - can also get by doing
% setfiltercellswith the cellinfo file

% Trajectories:
% 1 = Out Left    well 1to2
% 2 = In Left     well 2to1
% 3 = Out Right   well 1to3
% 4 = In Right    well 3to1


%day epoch tet cell
cellsi(1,:) = [day epoch 5 2];  % Traj 2 -In Left - Near Center well at trajend
cellsi(2,:) = [day epoch 9 1];  % Traj 3 -Out Right - Near Choice Point at trajmid
cellsi(3,:) = [day epoch 9 4];  % Traj 3 -Out Right - Near Right EndArm at traj9/10th
cellsi(4,:) = [day epoch 9 5];  % Traj 2 -In Left - Near Center well at trajend, some traj 3&4
cellsi(5,:) = [day epoch 12 1]; % Traj 3&1 -Out Right and Out Left - At Center well at trajstart
cellsi(6,:) = [day epoch 12 2]; % Traj 3&1 -Out Right and Out Left - At Center well at trajstart
cellsi(7,:) = [day epoch 12 3]; % Traj 2 -In Left - At Left arm Corner at trajmid
cellsi(8,:) = [day epoch 12 4]; % Traj 3&4 -Out Right and In Right - At Right arm center at traj 3/4th
cellsi(9,:) = [day epoch 12 5]; % Traj 1&2 -Out Left and In Left - Near Left Arm Corner at trajmid
cellsi(10,:) = [day epoch 12 6];% Traj 1 -Out Left - At Left arm center at traj3/4th

usecells = [2 3 5 6 8]; % 5 or 4 cells
% Reorder them
usecells = [6,2,8,3];


linposfile = sprintf('%s/%slinpos%02d.mat', animdirect, prefix, day);
load(linposfile);
statematrix = linpos{day}{epoch}.statematrix;
linpostime = statematrix.time;

posfile = sprintf('%s/%spos%02d.mat', animdirect, prefix, day);
load(posfile);
absvel = abs(pos{day}{epoch}.data(:,5)); % Can also use field 9
posn = pos{day}{epoch}.data(:,2:3); % Can also use fields 6:7
postime = pos{day}{epoch}.data(:,1); % same as linpostime

spikefile = sprintf('%s/%sspikes%02d.mat', animdirect, prefix, day);
load(spikefile);

DIOfile = sprintf('%s/%sDIO%02d.mat', animdirect, prefix, day);
load(DIOfile);
% Get DIO times
stim = DIO{day}{epoch}{16};
if isempty(stim)
    stim = DIO{day}{epoch}{15};
end
stim_starttime = stim.pulsetimes(:,1)./10000; %sec
stim_endtime = stim.pulsetimes(:,2)./10000; %sec


% EEG and if desired ripple band LFP data
cnt=0;
for tet=eegtets
    cnt=cnt+1;
    eegfile = sprintf('%s/EEG/%seeg%02d-%01d-%02d.mat', animdirect, prefix, day,epoch,tet);
    load(eegfile);
    eval(['eegdata',num2str(cnt),'= eeg{day}{epoch}{tet}.data;'])
    ripplefile = sprintf('%s/EEG/%sripple%02d-%01d-%02d.mat', animdirect, prefix, day,epoch,tet);
    load(ripplefile);
    eval(['rippledata',num2str(cnt),'= ripple{day}{epoch}{tet}.data;'])
    if cnt==1
        teeg = geteegtimes(eeg{day}{epoch}{tet});
        eegstart = eeg{day}{epoch}{tet}.starttime;
        eegsamprate = round(eeg{day}{epoch}{tet}.samprate);
    end
end

% Spike data
for i=1:size(cellsi,1)
    eval(['spiketime{',num2str(i),'}= spikes{cellsi(',num2str(i),',1)}{cellsi(',num2str(i),',2)}'...
        '{cellsi(',num2str(i),',3)}{cellsi(',num2str(i),',4)}.data;']);
    eval(['spiketime{',num2str(i),'}= spikes{cellsi(',num2str(i),',1)}{cellsi(',num2str(i),',2)}'...
        '{cellsi(',num2str(i),',3)}{cellsi(',num2str(i),',4)}.data(:,1);']);
    eval(['spikepos{',num2str(i),'}= spikes{cellsi(',num2str(i),',1)}{cellsi(',num2str(i),',2)}'...
        '{cellsi(',num2str(i),',3)}{cellsi(',num2str(i),',4)}.data(:,2:3);']);
    eval(['spikeposidx{',num2str(i),'}= spikes{cellsi(',num2str(i),',1)}{cellsi(',num2str(i),',2)}'...
        '{cellsi(',num2str(i),',3)}{cellsi(',num2str(i),',4)}.data(:,7);']);
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

% Look for traj 3: Out Right
% trajseq = trajseq(find(trajseq==3));
% start = start(find(trajseq==3));
% ntraj = length(start);
titlestr = ['OutRight'];

% If using subplots
fight = 0.09; % Height of each subplot as fraction
nplots = 9;

figno = 3;

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
    
    if currtraj == 3
        oritrajsttime = statematrix.time(currstart);
        oritrajendtime = statematrix.time(currend);
        oritimeintraj = oritrajendtime - oritrajsttime; % in sec
        trajstartidx = currstart;
        trajendidx = currend;
        
        
        trajsttime = statematrix.time(currstart)-5; % Add time to both
        trajendtime = statematrix.time(currend)+5;
        timeintraj = trajendtime - trajsttime; % in sec
        
        % Update start and end indices after adding secs to both ends
        currstart = lookup(trajsttime, statematrix.time);
        currend = lookup(trajendtime, statematrix.time);
        plotsttime = statematrix.time(currstart);
        plotendtime = statematrix.time(currend);
        
        if oritimeintraj > 2 % more than 2 secs
            
            % 1) Make a current position + spkposition figure
            %             figure(figno+100); hold on;
            %             redimscreen_figforppt1;
            %             plot(well1(1),well1(2),'ro','MarkerSize',16,'LineWidth',2);
            %             plot(well2(1),well2(2),'ro','MarkerSize',16,'LineWidth',2);
            %             currpos = posn(currstart:currend,:);
            %             plot(currpos(:,1),currpos(:,2),'b.-');
            
            figure(figno); hold on;
            redimscreen;
            
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
            title(['Traj',num2str(currtraj),' ',titlestr, ': ',num2str(round(trajsttime)),'-',num2str(round(trajendtime))],'FontSize',14,'Fontweight','bold')
            
            
            
            % 2) Then rasters Spikes
            
            % Get EEG Ind in current traj
            eind1 = lookup(trajsttime, teeg);
            eind2 = lookup(trajendtime, teeg);
            
            baseline = 4; % Max yplot for currspeed is 5. So go above that
  
            % --------
            clr = {'m','c','g','r','b','k','r','g','c','m'};
            starttime = teeg(eind1); % More accurate then position time due to better resln closer to spike resln
            endtime = teeg(eind2);
            cnt=0;
            for c=usecells
                eval(['currspkt = spiketime{',num2str(c),'};']);
                eval(['currspkpos = spikepos{',num2str(c),'};']);
                eval(['currspkposidx = spikeposidx{',num2str(c),'};']);
                currspkt = currspkt(find(currspkt>=teeg(eind1) & currspkt<=teeg(eind2)));
                % By time currstart:currend,:
                currspkpos1 = currspkpos(find(currspkt>=teeg(eind1) & currspkt<=teeg(eind2)),:);
                % Or by posidx
                currspkpos = currspkpos(find(currspkposidx>=currstart & currspkposidx<=currend),:);
                
                % If spikes, subtract from subtract from start time and bin
                if ~isempty(currspkt)
                    currspkt = currspkt - starttime;
                    %raster = starttime:0.001:endtime;
                end
                
                cnt=cnt+1;
                figure(figno); hold on; %subplot(nplots,1,cnt+5); hold on;
                if ~isempty(currspkt)
                    if size(currspkt,2)~=1, currspkt=currspkt'; end
                    % Use plotraster or spikeTrain
                    plotraster(currspkt,(baseline+2*(cnt-1))*ones(size(currspkt)),1.8,[],'Color',[clr{cnt}],'LineWidth',2);
                    %spikeTrain(currspkt,(baseline+(c-1))*ones(size(currspkt)),0.8,[],'Color',[clr{cnt}]);
                else
                    plot(0,0,'k.');
                end
                set(gca,'XLim',[0 endtime-starttime]);
                
                
                % Plot Spikes positions on top of all positions in figure 2
                %                 figure(figno+100); hold on;
                %                 plot(currspkpos(:,1),currspkpos(:,2),[clr{cnt} 'o'],'MarkerSize',8,'MarkerFaceColor',[clr{cnt}]);
                %                 plot(currspkpos1(:,1),currspkpos1(:,2),[clr{cnt} 'o'],'MarkerSize',8,'MarkerFaceColor',[clr{cnt}]);
                
            end
            
            % 3) Then Linearized Position on top of Spike Rasters
            % How far up on Yaxis do spike rasters go
            up = baseline+2*(cnt-1)+2;
            plotscale = up-baseline;
            % Plot a line at baseline and up for linearized position to end if you want
            plot(currpostime,baseline*ones(size(currspeed)),'k--','LineWidth',1);
            plot(currpostime,up*ones(size(currspeed)),'k--','LineWidth',1);
            % Get current linear distance
            lindist1 = statematrix.linearDistanceToWells(currstart:currend,currwells(1));
            lindist2 = statematrix.linearDistanceToWells(currstart:currend,currwells(2));
            distscale = max(lindist1)-min(lindist1);
            lindist1 = baseline + lindist1.*(plotscale/distscale);
            plot(currpostime,lindist1,'k-','LineWidth',2);
           
            
            % EEG
            %-------
            
            currDIOtime = stim_starttime(find((stim_starttime>starttime) & (stim_starttime<endtime)));
            
            %for i=1:length(eegtets)
            figure(figno); hold on; %subplot(nplots,1,i+1); hold on;
            eval(['eegdata = eegdata',num2str(maineegidx),';']);
            curreegdata = eegdata(eind1:eind2);
            eegscale = max(curreegdata)-min(curreegdata);
            downeeg = up+1;
            upeeg = downeeg+5;
            plotscale = 5;
            curreegdata = downeeg + (plotscale/2) + curreegdata.*(plotscale/eegscale);
            eegtimeaxis = teeg(eind1:eind2);
            eegtimeaxis = eegtimeaxis - teeg(eind1);
            
            plot(eegtimeaxis,curreegdata,'Linewidth',1);
            
            % DIOs within eeg
            y1=min(eegdata(eind1:eind2)); y2=max(eegdata(eind1:eind2));
            for s=1:length(currDIOtime)
                
                yaxis = 0:25;
                DIOt = currDIOtime(s)-teeg(eind1);
                speedt = currspeed(lookup(DIOt, currpostime));
                if speedt < 3*speedthrscale
                    %plot(DIOt*ones(size(yaxis)),yaxis,'m--','LineWidth',1.5);
                    xaxis = DIOt:0.01:DIOt+0.1;
                    jbfill(xaxis,25*ones(size(xaxis)),0*ones(size(xaxis)),'m','m',1,0.2);
                end
                
            end
            
            
            
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
            
            
            set(gca,'YLim',[0 18]);
            %set(gca,'YTick',[]); set(gca,'YTickLabel',[]);
            
            
            % Return to user
            keyboard;
            close(figno);
            %        close(figno+100);
            
        end  % end timeintraj
        
    end  % End current traj
    
    
    
    
    
end







