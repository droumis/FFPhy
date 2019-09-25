function out = sj_plottrajdatafind2(animdirect,prefix,index, varargin)
% sj_plottrajdatafind2('/data25/sjadhav/RippleInterruption/REf_direct/','REf',[7 4]);
% 4 EEGs and 4 Spikes with speed. Also - Spikes on Top of Trajectory


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


linposfile = sprintf('%s/%slinpos%02d.mat', animdirect, prefix, day);
load(linposfile);
statematrix = linpos{day}{epoch}.statematrix;

posfile = sprintf('%s/%spos%02d.mat', animdirect, prefix, day);
load(posfile);
absvel = abs(pos{day}{epoch}.data(:,5)); % Can also use field 9
posn = pos{day}{epoch}.data(:,2:3); % Can also use fields 6:7

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
usecells = [2 3 6 8]; % 5 cells

fight = 0.09; % Height of each subplot as fraction
nplots = 9;

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
        trajsttime = statematrix.time(currstart);
        trajendtime = statematrix.time(currend);
        oritimeintraj = trajendtime - trajsttime; % in sec 
        trajstartidx = currstart;
        trajendidx = currend;
        
        
        trajsttime = statematrix.time(currstart)-2; % Add time to both
        trajendtime = statematrix.time(currend)+2;
        timeintraj = trajendtime - trajsttime; % in sec 
        
        % Update start and end indices after adding secs to both ends
        currstart = lookup(trajsttime, statematrix.time);
        currend = lookup(trajendtime, statematrix.time);
        
        
        if oritimeintraj > 2 % more than 2 secs
            
            % Get DIO in current traj
            % ------------------------
            currDIOtime = stim_starttime(find((stim_starttime>trajsttime) & (stim_starttime<trajendtime)));
            
            % Make a current position + spkposition figure
            figure(1); hold on;
            redimscreen_figforppt1;
            %plot(well1(1),well1(2),'ro','MarkerSize',16,'MarkerFaceColor','r');
            %plot(well2(1),well2(2),'ro','MarkerSize',16,'MarkerFaceColor','r');
            plot(well1(1),well1(2),'ro','MarkerSize',16,'LineWidth',2);
            plot(well2(1),well2(2),'ro','MarkerSize',16,'LineWidth',2);
            currpos = posn(currstart:currend,:);
            plot(currpos(:,1),currpos(:,2),'b.-');          
            
            figure(2); hold on;
            redimscreen;
            %title(['Traj' num2str(currtraj) ': '
            %titlestr],'FontSize',12,'Fontweight','normal')
    
            % Speed
            %--------
            
            figure(2); subplot(nplots,1,1); hold on;
            currspeed = absvel(currstart:currend);
            plot(currspeed,'Linewidth',2);
            set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
            set(gca,'XTick',[]); set(gca,'YTick',[]);
            set(gca,'YLim',[min(currspeed) max(currspeed)]);
            set(gca,'XLim',[0 length(currspeed)]);
            %axis off
            
            yaxis = speedthrs1*ones(size(currspeed));
            xaxis = [1:length(currspeed)]';
            plot(xaxis,yaxis,'r');
            
            y1=min(currspeed); y2=max(currspeed);
            
            % Plot Traj Start and End points with circles
            plot(trajstartidx-currstart,y2*0.9,'ro','MarkerSize',8,'LineWidth',2);
            plot(trajendidx-currstart,y2*0.9,'ro','MarkerSize',8,'LineWidth',2);
            
            % DIOs within speed
            
            for s=1:length(currDIOtime)
                sind1 = lookup(currDIOtime(s), statematrix.time);
                sind2 = lookup(currDIOtime(s)+0.1, statematrix.time); % 100ms
                sind1 = sind1 - currstart;
                sind2 = sind2 - currstart;
                %jbfill
                xaxis=sind1:sind2;
                jbfill(xaxis,y2*ones(size(xaxis)),y1*ones(size(xaxis)),'m','m',1,0.2);
                %       % or dotted box
                %         % vertical lines
                %       yaxis=[y1:0.1:y2];
                %       xaxis=sind1*ones(size(y1));
                %       plot(xaxis,yaxis,'m-')
                %       xaxis=sind2*ones(size(y1));
                %       plot(xaxis,yaxis,'m-')
                %         % hor lines
                %       xaxis=sind1:sind2;
                %       yaxis=y1*ones(size(xaxis));
                %       plot(xaxis,yaxis,'m-')
                %       yaxis=y2*ones(size(xaxis));
                %       plot(xaxis,yaxis,'m-')
            end
            
            % Plot 1 sec scale bar
            ind1 = currstart-currstart;
            ind2 = lookup(trajsttime+1, statematrix.time);
            ind2 = ind2 - currstart;
            xaxis=ind1:ind2;
            yaxis=(y2)*ones(size(xaxis));
            plot(xaxis,yaxis,'k-','Linewidth',4);
            
            set(gca,'Position',[0.13 0.8646 0.775 fight]);
            
            title(['Traj',num2str(currtraj),' ',titlestr, ': ',num2str(round(trajsttime)),'-',num2str(round(trajendtime))],'FontSize',14,'Fontweight','bold')
            
            
            % Figure 1  - All positions in current traj
            
          
            % EEG
            %-------
            
            % Get EEG Ind in current traj
            eind1 = lookup(trajsttime, teeg);
            eind2 = lookup(trajendtime, teeg);
            
            for i=1:length(eegtets)
                figure(2); subplot(nplots,1,i+1); hold on;
                eval(['eegdata = eegdata',num2str(i),';']);
                plot(eegdata(eind1:eind2),'Linewidth',1);
                set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
                set(gca,'XTick',[]); set(gca,'YTick',[]);
                set(gca,'YLim',[min(eegdata(eind1:eind2)) max(eegdata(eind1:eind2))]);
                set(gca,'XLim',[0 eind2-eind1+20]);
                
                % DIOs within eeg
                y1=min(eegdata(eind1:eind2)); y2=max(eegdata(eind1:eind2));
                for s=1:length(currDIOtime)
                    sind1 = lookup(currDIOtime(s), teeg);
                    sind2 = lookup(currDIOtime(s)+0.1, teeg); % 100ms
                    sind1 = sind1 - eind1;
                    sind2 = sind2 - eind1;
                    %jbfill
                    xaxis=sind1:sind2;
                    jbfill(xaxis,y2*ones(size(xaxis)),y1*ones(size(xaxis)),'m','m',1,0.2);
                end
                set(gca,'Position',[0.13 0.8646-(fight+0.004)*(i) 0.775 fight]);
                %             if i==1
                %                 set(gca,'Position',[0.13 0.7840 0.775 0.0594]);
                %             end
                %             if i==2
                %                 set(gca,'Position',[0.13 0.72 0.775 0.0594])
                %             end
            end
            
            
            % Ripple
            %-------
%             for i=1:length(eegtets)
%                 figure(2); subplot(nplots,1,i+3); hold on;
%                 eval(['rippledata = rippledata',num2str(i),';']);
%                 plot(rippledata(eind1:eind2,1),'Linewidth',1);
%                 set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
%                 set(gca,'XTick',[]); set(gca,'YTick',[]);
%                 set(gca,'YLim',[min(rippledata(eind1:eind2,1)) max(rippledata(eind1:eind2,1))]);
%                 set(gca,'XLim',[0 eind2-eind1+1]);
%                 % DIOs within ripple
%                 y1=double(min(rippledata(eind1:eind2,1))); y2=double(max(rippledata(eind1:eind2,1)));
%                 for s=1:length(currDIOtime)
%                     sind1 = lookup(currDIOtime(s), teeg);
%                     sind2 = lookup(currDIOtime(s)+0.1, teeg); % 100ms
%                     sind1 = sind1 - eind1;
%                     sind2 = sind2 - eind1;
%                     %jbfill
%                     xaxis=sind1:sind2;
%                     jbfill(xaxis,y2*ones(size(xaxis)),y1*ones(size(xaxis)),'m','m',1,0.2);
%                 end
%                 set(gca,'Position',[0.13 0.8646-(fight+0.004)*(i+2) 0.775 fight]);
%                 %             if i==1
%                 %                 set(gca,'Position',[0.13 0.656 0.775 0.0594]);
%                 %             end
%                 %             if i==2
%                 %                 set(gca,'Position',[0.13 0.592 0.775 0.0594])
%                 %             end
%             end
            
            % Spikes
            % --------
            clr = {'k','m','g','c','r','k','r','g','c','m'};
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
                    raster = starttime:0.001:endtime;
                end
                
                %             if i<6
                %                 subplot(10,1,6); hold on;
                %             else
                %                 subplot(8,1,8); hold on;
                %             end
                cnt=cnt+1;
                figure(2); subplot(nplots,1,cnt+5); hold on;
                if ~isempty(currspkt)
                    if size(currspkt,2)~=1, currspkt=currspkt'; end
                    plotraster(currspkt,ones(size(currspkt)),0.8,[],'Color',[clr{cnt}]);
                else
                    plot(0,0,'k.');
                end
                set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
                set(gca,'XTick',[]); set(gca,'YTick',[]);
                set(gca,'YLim',[1 1.8]);
                set(gca,'XLim',[0 endtime-starttime]);
                set(gca,'Position',[0.13 0.8646-(fight+0.004)*(cnt+4) 0.775 fight]);
                
                
                % Plot Spikes positions on top of all positions
                figure(1); hold on;
                plot(currspkpos(:,1),currspkpos(:,2),[clr{cnt} 'o'],'MarkerSize',8,'MarkerFaceColor',[clr{cnt}]); 
                %plot(currspkpos1(:,1),currspkpos1(:,2),[clr{cnt} 'o'],'MarkerSize',8,'MarkerFaceColor',[clr{cnt}]); 
                
            end
            
           
            
            
            
            % Return to user
            keyboard;
            close(1);
            close(2);
            
        end  % end timeintraj
        
    end  % End current traj
    
    
    
    
    
end







