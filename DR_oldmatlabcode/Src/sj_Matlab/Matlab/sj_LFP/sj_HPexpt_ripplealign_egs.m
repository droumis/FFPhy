function out = sj_HPexpt_ripplealign_egs(prefix, day, epoch, varargin)
% Shantanu - Nov 2012
% Align and plot LFP and spiking to single ripples. From sj_plottrajdatafind3.m
% sj_HPexpt_ripplealign_egs('HPb', 6, 2,'docells',0,'dopos',0,'plotrippleband',1);

if nargin<1,
    keyboard
    error('Please enter Expt Prefix and Day No!');
end
if nargin<2,
    keyboard
    error('Please enter Day No!');
end
if nargin<3,
    epoch=1; %% Epoch - Usually 1 for calibration
end

% Some parameters
Fspos = 30; %Hz
respos = 1/30; % sec
Fseeg = 1500; %Hz
reseeg = 1/1500; % sec
Fsspikes = 10000; %Hz
resspikes = 1/10000; %sec
speedthrs1 = 3; %cm/sec;

docells=0; dopos=0;plotrippleband=0;

pret=100; postt=200; %in ms

for option = 1:2:length(varargin)-1
    switch varargin{option}
        case 'docells'
            docells = varargin{option+1};
        case 'dopos'
            dopos = varargin{option+1};
        case 'plotrippleband'
            plotrippleband = varargin{option+1};
    end
end



set(0,'defaultaxesfontweight','normal'); set(0,'defaultaxeslinewidth',2);
forppr=0;
if forppr==1
    set(0,'defaultaxesfontsize',16);
    tfont = 18; % title font
    xfont = 16;
    yfont = 16;
else
    set(0,'defaultaxesfontsize',24);
    tfont = 28;
    xfont = 20;
    yfont = 20;
end
clr = {'b','r','g','c','m','y','k','r'};

%% -----------------------------------------
% SET DATA
% -------------------------------------------

switch prefix
    case 'HPa'
        rawdir = '/data25/sjadhav/HPExpt/HPa/';
        directoryname = '/data25/sjadhav/HPExpt/HPa_direct/';
    case 'HPb'
        rawdir = '/data25/sjadhav/HPExpt/HPb/';
        directoryname = '/data25/sjadhav/HPExpt/HPb_direct/';
        
        switch day
            case 1
                riptet = 4; % 18 for iHp
                dHp_eegtet = 4;
                iHp_eegtet = 6;
                PFC_eegtet = 9;
                reftet = 7;
                eegtets = [dHp_eegtet, iHp_eegtet, PFC_eegtet, reftet];
                
                % For cells
                dHptets = [1,2,3,4,5,6];
                iHptets = 18; % All other cells are silent
                PFCtets = [8,9,10,12,14];
                reftag = 'Ref7Hp'; % for PFC wrt to Hipp Ref
                
             case 6
                riptet = 4; % 18 for iHp
                dHp_eegtet = 4;
                iHp_eegtet = 6;
                PFC_eegtet = 9;
                reftet = 7;
                eegtets = [dHp_eegtet, iHp_eegtet, PFC_eegtet, reftet];
                reftag = 'Ref7Hp'; % for PFC wrt to Hipp Ref
                
                
        end
        
end

dir2 = directoryname;
animdirect=dir2;
clr = {'b','g','c','m','y','k','r'};

% --------------------------
% Load extracted ripple file -
% --------------------------
% SHOULD SWITCH THIS TO ALLTET. If alltet does not exist, then use given tet/tets.
ripfile = sprintf('%s/%sripples%02d.mat', directoryname, prefix, day);
load(ripfile);
rip_starttime=[]; rip_sizes=[];
for i=1:length(riptet)
    currriptet=riptet(i);
    rip_starttime = [rip_starttime; 1000* ripples{day}{epoch}{currriptet}.starttime];   % in msec
    rip_sizes = [rip_sizes; ripples{day}{epoch}{currriptet}.maxthresh];   % in units of std dev
end
%rem = find(rip_sizes<3 & rip_sizes>4);
rem = find(rip_sizes>5);
rip_starttime(rem) = [];
rip_sizes(rem) = [];
[rip_starttime,sortidx] = sort(rip_starttime);
rip_sizes = rip_sizes(sortidx);
% Define triggering events as the start of each ripple
triggers = rip_starttime;

% --------------------------
% Get Other Data -
% --------------------------

if dopos==1
    
    % linposfile = sprintf('%s/%slinpos%02d.mat', animdirect, prefix, day);
    % load(linposfile);
    % statematrix = linpos{day}{epoch}.statematrix;
    % linpostime = statematrix.time;
    
    posfile = sprintf('%s/%spos%02d.mat', animdirect, prefix, day);
    load(posfile);
    absvel = abs(pos{day}{epoch}.data(:,5)); % Can also use field 9
    posn = pos{day}{epoch}.data(:,2:3); % Can also use fields 6:7
    postime = pos{day}{epoch}.data(:,1); % same as linpostime
    postime = postime*1000; % in ms
    
end



%DIOfile = sprintf('%s/%sDIO%02d.mat', animdirect, prefix, day);
%load(DIOfile);
% Get DIO times - eg. for reward
%stim = DIO{day}{epoch}{6};
% stim_starttime = stim.pulsetimes(:,1)./10000; %sec
% stim_endtime = stim.pulsetimes(:,2)./10000; %sec

if (day<10)
    daystring = ['0',num2str(day)];
else
    daystring = num2str(day);
end

% ---------------------------------------
% EEG / Gnd EEG / and if desired ripple band LFP data
% ----------------------------------------
cnt=0;
for tet=eegtets
    cnt=cnt+1;
     if (tet<10)
          tetstring = ['0',num2str(tet)];
     else
          tetstring = num2str(tet);
     end
        
    % EEG
    % ----
    
    % Check if REFTAG file exists - this is especially for PFC tetrodes
    eegreffile = [dir2,'/EEG/',prefix,'eeg',reftag,daystring,'-',num2str(epoch),'-',tetstring];
    if (exist([eegreffile,'.mat'],'file'))==2
        disp(['REFTAG reference file exist for Tet ',num2str(tet),'; so using that.'])
        load(eegreffile);
        eval(['eegdata',num2str(cnt),'= eegref{day}{epoch}{tet}.data;'])
    else
        eegfile = sprintf('%sEEG/%seeg%02d-%01d-%02d.mat', animdirect, prefix, day,epoch,tet);
        load(eegfile);
        eval(['eegdata',num2str(cnt),'= eeg{day}{epoch}{tet}.data;'])
    end
    % EEG GND
    % --------
    eeggndfile = sprintf('%sEEG/%seeggnd%02d-%01d-%02d.mat', animdirect, prefix, day,epoch,tet);
    if (exist([eeggndfile],'file'))==2
        load(eeggndfile);
        eval(['eeggnddata',num2str(cnt),'= eeggnd{day}{epoch}{tet}.data;'])
    else
         disp(['EEG Gnd File does not exist for Tet ',num2str(tet),'; Getting regular file']);
         eval(['eeggnddata',num2str(cnt),'= eeg{day}{epoch}{tet}.data;'])
    end
    
    % Ripple band
    % -----------
    ripplefile = sprintf('%sEEG/%sripple%02d-%01d-%02d.mat', animdirect, prefix, day,epoch,tet);
    if (exist([ripplefile],'file'))==2
        load(ripplefile);
        eval(['rippledata',num2str(cnt),'= ripple{day}{epoch}{tet}.data;'])
    else
        disp(['Ripple File does not exist for Tet ',num2str(tet)]);
    end
    
    % Ripple Band Gnd
    % --------------
    ripplegndfile = sprintf('%sEEG/%sripplegnd%02d-%01d-%02d.mat', animdirect, prefix, day,epoch,tet);
    if (exist([ripplegndfile],'file'))==2
        load(ripplegndfile);
        eval(['ripplegnddata',num2str(cnt),'= ripplegnd{day}{epoch}{tet}.data;'])
    else
        disp(['Ripple GND File does not exist for Tet ',num2str(tet),'; Getting regular file']);
        load(ripplefile);
        eval(['ripplegnddata',num2str(cnt),'= ripple{day}{epoch}{tet}.data;'])
    end
        
    if cnt==1
        teeg = geteegtimes(eeg{day}{epoch}{tet})*1000; % in msec
        eegstart = eeg{day}{epoch}{tet}.starttime;
        eegsamprate = round(eeg{day}{epoch}{tet}.samprate);
    end
end

if docells==1
    
    spikefile = sprintf('%s/%sspikes%02d.mat', animdirect, prefix, day);
    load(spikefile);
    
    % Get cell indices
    n_dHpcells = 0; n_iHpcells = 0; n_PFCcells = 0;
    dHp_idxs = []; iHp_idxs = []; PFC_idxs = [];
    
    for i=1:length(dHptets),
        currtet = dHptets(i);
        ncurrcells = length(spikes{day}{epoch}{currtet});
        for c=1:ncurrcells
            if ~isempty(spikes{day}{epoch}{currtet}{c})
                if ~isempty(spikes{day}{epoch}{currtet}{c}.data)
                    if strcmp(spikes{day}{epoch}{currtet}{c}.tag,'MU')==0
                        n_dHpcells = n_dHpcells+1;
                        dHp_idxs = [dHp_idxs; [currtet, c]];
                    end
                end
            end
        end
    end
    
    for i=1:length(iHptets),
        currtet = iHptets(i);
        ncurrcells = length(spikes{day}{epoch}{currtet});
        for c=1:ncurrcells
            if ~isempty(spikes{day}{epoch}{currtet}{c})
                if ~isempty(spikes{day}{epoch}{currtet}{c}.data)
                    if strcmp(spikes{day}{epoch}{currtet}{c}.tag,'MU')==0
                        n_iHpcells = n_iHpcells+1;
                        iHp_idxs = [iHp_idxs; [currtet, c]];
                    end
                end
            end
        end
    end
    
    for i=1:length(PFCtets),
        currtet = PFCtets(i);
        ncurrcells = length(spikes{day}{epoch}{currtet});
        for c=1:ncurrcells
            if ~isempty(spikes{day}{epoch}{currtet}{c})
                if ~isempty(spikes{day}{epoch}{currtet}{c}.data)
                    if strcmp(spikes{day}{epoch}{currtet}{c}.tag,'MU')==0
                        n_PFCcells = n_PFCcells+1;
                        PFC_idxs = [PFC_idxs; [currtet, c]];
                    end
                end
            end
        end
    end
    
    allspikeidxs = [dHp_idxs; iHp_idxs; PFC_idxs];

end


% Set fig prop
set(0,'defaultaxesfontsize',12);set(0,'defaultaxesfontweight','normal');
set(0,'defaultaxeslinewidth',1);

% If using subplots
fight = 0.09; % Height of each subplot as fraction
nplots = 9; figno = 3;


for nr = 1:length(triggers),
    
    figure; hold on;
    %set(gcf,'Position',[120 150 900 930]);
    
    set(gcf,'Position',[45 652 900 430]);

    currtime = triggers(nr); % Rip starttime in ms
    sttime = currtime - pret;
    endtime = currtime + postt;
    
    % EEG
    % ----
    % Get EEG Ind for current ripple
    eind1 = lookup(sttime, teeg);
    eind2 = lookup(endtime, teeg);
    
    currdHpeeg = eeggnddata1(eind1:eind2);
    currPFCeeg = eeggnddata3(eind1:eind2);
    currdHpeegr = eegdata1(eind1:eind2);
    currPFCeegr = eegdata3(eind1:eind2);
    
    currrefeeg = eeggnddata4(eind1:eind2); % Ref is only wrt to gnd
    
    if plotrippleband == 1
        currdHpripg = double(ripplegnddata1(eind1:eind2,1)); % 1st column is filtered amplitude,2nd is phase, 3rd is envelope_magnitude
        currPFCripg = double(ripplegnddata3(eind1:eind2,1));
        currdHpripr = double(rippledata1(eind1:eind2,1)); % 1st column is filtered amplitude,2nd is phase, 3rd is envelope_magnitude
        currPFCripr = double(rippledata3(eind1:eind2,1));
        
        currrefrip = double(rippledata4(eind1:eind2,1));
        
    end
    
    subplot(3,1,1); hold on;
    eegtimeaxis = teeg(eind1:eind2); eegtimeaxis = eegtimeaxis - currtime;
    plot(eegtimeaxis,currdHpeeg,'k','Linewidth',2);
    plot(eegtimeaxis,currdHpeegr,'b','Linewidth',2);
    mina = min(min(currdHpeeg),min(currdHpeegr)); maxa = max(max(currdHpeeg),max(currdHpeegr));
    set(gca,'XLim',[min(eegtimeaxis) max(eegtimeaxis)]); set(gca,'YLim',1.1*[mina maxa]);
    ypts = mina:0.01:maxa; xpts = 0*ones(size(ypts));
    plot(xpts , ypts, 'k--','Linewidth',1);
    %axis off
    
    subplot(3,1,2); hold on;
    eegtimeaxis = teeg(eind1:eind2); eegtimeaxis = eegtimeaxis - currtime;
    plot(eegtimeaxis,currPFCeeg,'k','Linewidth',2);
    plot(eegtimeaxis,currPFCeegr,'r','Linewidth',2);
    mina = min(min(currPFCeeg),min(currPFCeegr)); maxa = max(max(currPFCeeg),max(currPFCeegr));
    set(gca,'XLim',[min(eegtimeaxis) max(eegtimeaxis)]); set(gca,'YLim',1.1*[mina maxa]);
    ypts = mina:0.01:maxa; xpts = 0*ones(size(ypts));
    plot(xpts , ypts, 'k--','Linewidth',1);
    %axis off
    
    subplot(3,1,3); hold on;
    eegtimeaxis = teeg(eind1:eind2); eegtimeaxis = eegtimeaxis - currtime;
    plot(eegtimeaxis,currrefeeg,'k','Linewidth',2);
    mina = min(currrefeeg); maxa = max(currrefeeg);
    set(gca,'XLim',[min(eegtimeaxis) max(eegtimeaxis)]); set(gca,'YLim',1.1*[mina maxa]);
    ypts = mina:0.01:maxa; xpts = 0*ones(size(ypts));
    plot(xpts , ypts, 'k--','Linewidth',1);
    %axis off
    
    if plotrippleband == 1
        figure; hold on;
        set(gcf,'Position',[44 124 900 430]);
        subplot(3,1,1); hold on;
        eegtimeaxis = teeg(eind1:eind2); eegtimeaxis = eegtimeaxis - currtime;
        plot(eegtimeaxis,currdHpripg,'k','Linewidth',2);
        plot(eegtimeaxis,currdHpripr,'b','Linewidth',2);
        mina = min(min(currdHpripg),min(currdHpripr)); maxa = max(max(currdHpripg),max(currdHpripr));
        set(gca,'XLim',[min(eegtimeaxis) max(eegtimeaxis)]); set(gca,'YLim',1.1*[mina maxa]);
        ypts = mina:0.01:maxa; xpts = 0*ones(size(ypts));
        plot(xpts , ypts, 'k--','Linewidth',1);
        %axis off
        
        subplot(3,1,2); hold on;
        eegtimeaxis = teeg(eind1:eind2); eegtimeaxis = eegtimeaxis - currtime;
        plot(eegtimeaxis,currPFCripg,'k','Linewidth',2);
        plot(eegtimeaxis,currPFCripr,'r','Linewidth',2);
        mina = min(min(currPFCripg),min(currPFCripr)); maxa = max(max(currPFCripg),max(currPFCripr));
        set(gca,'XLim',[min(eegtimeaxis) max(eegtimeaxis)]); set(gca,'YLim',1.1*[mina maxa]);
        ypts = mina:0.01:maxa; xpts = 0*ones(size(ypts));
        plot(xpts , ypts, 'k--','Linewidth',1);
        %axis off
        
        subplot(3,1,3); hold on;
        eegtimeaxis = teeg(eind1:eind2); eegtimeaxis = eegtimeaxis - currtime;
        plot(eegtimeaxis,currrefrip,'k','Linewidth',2);
        mina = min(currrefrip); maxa = max(currrefrip);
        set(gca,'XLim',[min(eegtimeaxis) max(eegtimeaxis)]); set(gca,'YLim',1.1*[mina maxa]);
        ypts = mina:0.01:maxa; xpts = 0*ones(size(ypts));
        plot(xpts , ypts, 'k--','Linewidth',1);
        %axis off
        
    end
    
%     % Spike Raster
%     subplot(7,1,[3 4 5 6]); hold on;
%     
%     for i=1:length(allspikeidxs)
%         spikeu = spikes{day}{epoch}{allspikeidxs(i,1)}{allspikeidxs(i,2)}.data*1000;
%         currspks =  spikeu(find( (spikeu>=sttime) & (spikeu<=endtime) ));
%         currspks = currspks - currtime;
%         if i <= (n_dHpcells + n_iHpcells)
%             clr = 'b';
%         else
%             clr = 'r';
%         end
%         plotraster(currspks,i*ones(size(currspks)), 0.8,[],'LineWidth',3,'Color',clr);
%     end
%     set(gca,'YLim',[0.5 i+1.5]);
%     ypts = 0:1:i+1; xpts = 0*ones(size(ypts));
%     plot(xpts , ypts, 'k--','Linewidth',1);
%     
%     % Speed
%     subplot(7,1,[7]); hold on;
%     
%     pind1 = lookup(sttime, postime);
%     pind2 = lookup(endtime, postime);
%     currspeed = absvel(pind1:pind2);
%     ptimeaxis = postime(pind1:pind2); ptimeaxis = ptimeaxis - currtime;
%     plot(ptimeaxis,currspeed,'k','Linewidth',2);
%     set(gca,'XLim',[min(ptimeaxis) max(ptimeaxis)]); set(gca,'YLim',[min(currspeed) max(currspeed)]);
%     ypts = min(currspeed):0.01:max(currspeed); xpts = 0*ones(size(ypts));
%     plot(xpts , ypts, 'k--','Linewidth',1);


    %axis off

    %   currpostime = statematrix.time(currstart:currend);
%   currpostime = currpostime - statematrix.time(currstart)
    
    
    keyboard;
    
    close all;
end








% % Plot each trajectory in a loop
% for n=1:ntraj-1
%     
%     currstart = start(n);
%     currend = start(n+1);
%     currtraj = trajseq(n);
%     
%     if currtraj>=3,
%         currwells = trajwells(2,:);
%     else
%         currwells = trajwells(1,:);
%     end
%     well1 = wellCoord(currwells(1),:); %xy position of wells
%     well2 = wellCoord(currwells(2),:);
%     
%     % Only do if currtraj=3
%     
%     if currtraj == 3
%         oritrajsttime = statematrix.time(currstart);
%         oritrajendtime = statematrix.time(currend);
%         oritimeintraj = oritrajendtime - oritrajsttime; % in sec
%         trajstartidx = currstart;
%         trajendidx = currend;
%         
%         
%         trajsttime = statematrix.time(currstart)-5; % Add time to both
%         trajendtime = statematrix.time(currend)+5;
%         timeintraj = trajendtime - trajsttime; % in sec
%         
%         % Update start and end indices after adding secs to both ends
%         currstart = lookup(trajsttime, statematrix.time);
%         currend = lookup(trajendtime, statematrix.time);
%         plotsttime = statematrix.time(currstart);
%         plotendtime = statematrix.time(currend);
%         
%         if oritimeintraj > 2 % more than 2 secs
%             
%             % 1) Make a current position + spkposition figure
%             %             figure(figno+100); hold on;
%             %             redimscreen_figforppt1;
%             %             plot(well1(1),well1(2),'ro','MarkerSize',16,'LineWidth',2);
%             %             plot(well2(1),well2(2),'ro','MarkerSize',16,'LineWidth',2);
%             %             currpos = posn(currstart:currend,:);
%             %             plot(currpos(:,1),currpos(:,2),'b.-');
%             
%             figure(figno); hold on;
%             redimscreen;
%             
%             % Speed
%             %--------
%             
%             % Speed ylimits in plot: 0-5
%             
%             figure(figno); %subplot(nplots,1,1); hold on;
%             hold on;
%             currspeed = absvel(currstart:currend);
%             currpostime = statematrix.time(currstart:currend);
%             currpostime = currpostime - statematrix.time(currstart);
%             
%             % Change scale
%             speedscale = max(currspeed)-min(currspeed);
%             plotscale = 3-0;
%             currspeed = currspeed.*(plotscale/speedscale);
%             
%             plot(currpostime, currspeed,'Linewidth',2);
%             set(gca,'XLim',[0 max(currpostime)]);
%             
%             % Line for speedthrs
%             yaxis = speedthrs1*(plotscale/speedscale)*1.5*ones(size([0:length(currspeed)]'));
%             xaxis = [0:length(currspeed)]';
%             plot(xaxis,yaxis,'r');
%             speedthrscale = speedthrs1*(plotscale/speedscale);
%             
%             y1=min(currspeed); y2=max(currspeed);
%             
%             % Plot Traj Start and End points with circles
%             plott = oritrajsttime-plotsttime;
%             speedt = currspeed(lookup(plott, currpostime));
%             plot(plott,speedt,'ro','MarkerSize',8,'LineWidth',2);
%             plott = oritrajendtime-plotsttime;
%             speedt = currspeed(lookup(plott, currpostime));
%             plot(plott,speedt,'ro','MarkerSize',8,'LineWidth',2);
%             
%             % Plot 1 sec scale bar
%             xaxis = 0:0.1:1;
%             yaxis=(y2)*0.9*ones(size(xaxis));
%             plot(xaxis,yaxis,'k-','Linewidth',4);
%             title(['Traj',num2str(currtraj),' ',titlestr, ': ',num2str(round(trajsttime)),'-',num2str(round(trajendtime))],'FontSize',14,'Fontweight','bold')
%             
%             
%             
%             % 2) Then rasters Spikes
%             
%             % Get EEG Ind in current traj
%             eind1 = lookup(trajsttime, teeg);
%             eind2 = lookup(trajendtime, teeg);
%             
%             baseline = 4; % Max yplot for currspeed is 5. So go above that
%             
%             % --------
%             clr = {'m','c','g','r','b','k','r','g','c','m'};
%             starttime = teeg(eind1); % More accurate then position time due to better resln closer to spike resln
%             endtime = teeg(eind2);
%             cnt=0;
%             for c=usecells
%                 eval(['currspkt = spiketime{',num2str(c),'};']);
%                 eval(['currspkpos = spikepos{',num2str(c),'};']);
%                 eval(['currspkposidx = spikeposidx{',num2str(c),'};']);
%                 currspkt = currspkt(find(currspkt>=teeg(eind1) & currspkt<=teeg(eind2)));
%                 % By time currstart:currend,:
%                 currspkpos1 = currspkpos(find(currspkt>=teeg(eind1) & currspkt<=teeg(eind2)),:);
%                 % Or by posidx
%                 currspkpos = currspkpos(find(currspkposidx>=currstart & currspkposidx<=currend),:);
%                 
%                 % If spikes, subtract from subtract from start time and bin
%                 if ~isempty(currspkt)
%                     currspkt = currspkt - starttime;
%                     %raster = starttime:0.001:endtime;
%                 end
%                 
%                 cnt=cnt+1;
%                 figure(figno); hold on; %subplot(nplots,1,cnt+5); hold on;
%                 if ~isempty(currspkt)
%                     if size(currspkt,2)~=1, currspkt=currspkt'; end
%                     % Use plotraster or spikeTrain
%                     plotraster(currspkt,(baseline+2*(cnt-1))*ones(size(currspkt)),1.8,[],'Color',[clr{cnt}],'LineWidth',2);
%                     %spikeTrain(currspkt,(baseline+(c-1))*ones(size(currspkt)),0.8,[],'Color',[clr{cnt}]);
%                 else
%                     plot(0,0,'k.');
%                 end
%                 set(gca,'XLim',[0 endtime-starttime]);
%                 
%                 
%                 % Plot Spikes positions on top of all positions in figure 2
%                 %                 figure(figno+100); hold on;
%                 %                 plot(currspkpos(:,1),currspkpos(:,2),[clr{cnt} 'o'],'MarkerSize',8,'MarkerFaceColor',[clr{cnt}]);
%                 %                 plot(currspkpos1(:,1),currspkpos1(:,2),[clr{cnt} 'o'],'MarkerSize',8,'MarkerFaceColor',[clr{cnt}]);
%                 
%             end
%             
%             % 3) Then Linearized Position on top of Spike Rasters
%             % How far up on Yaxis do spike rasters go
%             up = baseline+2*(cnt-1)+2;
%             plotscale = up-baseline;
%             % Plot a line at baseline and up for linearized position to end if you want
%             plot(currpostime,baseline*ones(size(currspeed)),'k--','LineWidth',1);
%             plot(currpostime,up*ones(size(currspeed)),'k--','LineWidth',1);
%             % Get current linear distance
%             lindist1 = statematrix.linearDistanceToWells(currstart:currend,currwells(1));
%             lindist2 = statematrix.linearDistanceToWells(currstart:currend,currwells(2));
%             distscale = max(lindist1)-min(lindist1);
%             lindist1 = baseline + lindist1.*(plotscale/distscale);
%             plot(currpostime,lindist1,'k-','LineWidth',2);
%             
%             
%             % EEG
%             %-------
%             
%             currDIOtime = stim_starttime(find((stim_starttime>starttime) & (stim_starttime<endtime)));
%             
%             %for i=1:length(eegtets)
%             figure(figno); hold on; %subplot(nplots,1,i+1); hold on;
%             eval(['eegdata = eegdata',num2str(maineegidx),';']);
%             curreegdata = eegdata(eind1:eind2);
%             eegscale = max(curreegdata)-min(curreegdata);
%             downeeg = up+1;
%             upeeg = downeeg+5;
%             plotscale = 5;
%             curreegdata = downeeg + (plotscale/2) + curreegdata.*(plotscale/eegscale);
%             eegtimeaxis = teeg(eind1:eind2);
%             eegtimeaxis = eegtimeaxis - teeg(eind1);
%             
%             plot(eegtimeaxis,curreegdata,'Linewidth',1);
%             
%             % DIOs within eeg
%             y1=min(eegdata(eind1:eind2)); y2=max(eegdata(eind1:eind2));
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
%             
%             
%             
%             % Ripple
%             %-------
%             %         figure(figno); hold on; %subplot(nplots,1,i+1); hold on;
%             %         eval(['rippledata = rippledata',num2str(maineegidx),';']);
%             %         curreegdata = double(rippledata(eind1:eind2,1));
%             %         eegscale = max(curreegdata)-min(curreegdata);
%             %         plotscale = 29-24;
%             %         curreegdata =  24 + curreegdata.*(plotscale/eegscale);
%             %
%             %         plot(eegtimeaxis,curreegdata,'Linewidth',1);
%             
%             
%             set(gca,'YLim',[0 18]);
%             %set(gca,'YTick',[]); set(gca,'YTickLabel',[]);
%             
%             
%             % Return to user
%             keyboard;
%             close(figno);
%             %        close(figno+100);
%             
%         end  % end timeintraj
%         
%     end  % End current traj
%     
% end
    
    
    
    
    
    
    
    
