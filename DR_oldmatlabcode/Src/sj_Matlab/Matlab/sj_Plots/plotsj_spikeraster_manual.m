
% SJ  From plotsj_rip_spikeraster_Rolling. FOr plotting final example
% SJ - 12/29/2012. From plotsj_rip_spikeraster_LongEEG
% Instead of aligning to an event and plotting spike raster + EEG +
% Position, do it in a long rolling window for the given epoch

%function [pret] = plotsj_rip_spikeraster(prefix, day, epoch, saveg1)
% Shantanu- Dec 2012
% Win - Adding position or velocity to graph to see trajectory clearly

% Adapted from sj_ripegs1 and sj_plotrajdatafind3
% Plot raster of spikes during ripple using given cells. Not calling in function command right now.
% Can also plot LFP and posn. For these, esp. posn, check sj_plottrajdatafind3
% Shantanu 08Jun2012

% if nargin<1,
%     keyboard
%     error('Please enter Expt Prefix and Day No!');
% end
% if nargin<2,
%     keyboard
%     error('Please enter Day No!');
% end
% if nargin<3,
%     epoch=2; %% Epoch - 2 or 4 for runs
% end
% if nargin<4
%     saveg1=0; % Save summary figure
% end

clear; %close all;
prefix='HPa';
day=2; epoch=4;

win = 60; % 60 sec window
overlap = 10;

%sd=3; %% SD threshold for ripples
% Pret and Post from start of ripple? Also plot start-middle and end of ripple
%pret=30000;
%postt=30000;

n=2;


switch prefix
    case 'HPa'
        directoryname = 'D:/data25/sjadhav/HPExpt/HPa_direct';
        dire = 'D:/data25/sjadhav/HPExpt/HPa';
        riptetlist=[1,4,5,6];
end

currdir = pwd;
if (directoryname(end) == '/')
    directoryname = directoryname(1:end-1);
end
if (dire(end) == '/')
    dire = dire(1:end-1);
end

if (day < 10)
    daystring = ['0',num2str(day)];
else
    daystring = num2str(day);
end

animdirect = directoryname;

% Get time ranges from times file - No need
%----------------------------------
cd(dire);
dirlist = dir('*');
for i=2:length(dirlist)
    currname = dirlist(i).name;
    if strcmp(currname(1:2),daystring)
        daydir = currname;
        break;
    end
end
timesfile = sprintf('%s/%s/times.mat',dire,daydir);
load(timesfile);
nranges=1;
range1=ranges(epoch+1,:);

% --------------- Parameters ---------------


Fspos = 30; %Hz
respos = 1/30; % sec
Fseeg = 1500; %Hz
reseeg = 1/1500; % sec
Fsspikes = 10000; %Hz
resspikes = 1/10000; %sec

%% -----------------------------------------
% SET DATA
% -------------------------------------------
eegtets = riptetlist;
maineegtet = 1;
maineegidx=1;

% Also get a PFC eeg
peegtet = 15;
eegtets = [eegtets, peegtet];
peegidx=length(eegtets);

% HPa - day 2. CA1 cells 6+4; PFC cells 4+4
% CA1
% cellsi(1,:) = [day epoch 1 1];
% cellsi(2,:) = [day epoch 1 2];
% cellsi(3,:) = [day epoch 1 3];
% cellsi(4,:) = [day epoch 1 4];
% cellsi(5,:) = [day epoch 1 5];
% cellsi(6,:) = [day epoch 1 6];
% cellsi(7,:) = [day epoch 4 1];
% cellsi(8,:) = [day epoch 4 2];
% cellsi(9,:) = [day epoch 4 3];
% cellsi(10,:) = [day epoch 4 4];
% cellsi(11,:) = [day epoch 4 5];
% cellsi(12,:) = [day epoch 4 6];
% cellsi(13,:) = [day epoch 7 1];
% cellsi(14,:) = [day epoch 8 1];
% cellsi(15,:) = [day epoch 12 1];
% cellsi(16,:) = [day epoch 12 2];
% cellsi(17,:) = [day epoch 14 1];
% cellsi(18,:) = [day epoch 14 2];
% cellsi(19,:) = [day epoch 14 3];
% usecellsi=1:19;

% Putative 15 cells
% cellsi(1,:) = [day epoch 1 3];
% cellsi(2,:) = [day epoch 12 2]; % iHp
% cellsi(3,:) = [day epoch 8 1]; % iHp
% cellsi(4,:) = [day epoch 1 1];
% cellsi(5,:) = [day epoch 1 5];
% cellsi(6,:) = [day epoch 4 5];
% cellsi(7,:) = [day epoch 4 6];
% cellsi(8,:) = [day epoch 4 4];
% cellsi(9,:) = [day epoch 4 3];
% cellsi(10,:) = [day epoch 14 1]; %iHp
% cellsi(11,:) = [day epoch 1 6];
% cellsi(12,:) = [day epoch 1 4];
% cellsi(13,:) = [day epoch 4 1];
% cellsi(14,:) = [day epoch 7 1];
% cellsi(15,:) = [day epoch 1 2];
% usecellsi=1:15;


% % Ordered - For 7758.6-7722 (11 or 12 or 13 cells)
% %%%cellsi(1,:) = [day epoch 1 3]; % Too many
% cellsi(1,:) = [day epoch 12 2]; % iHp
% %%% cellsi(2,:) = [day epoch 8 1]; % iHp
% cellsi(2,:) = [day epoch 1 1];
% cellsi(3,:) = [day epoch 1 5];
% cellsi(4,:) = [day epoch 4 5];
% cellsi(5,:) = [day epoch 4 6];
% cellsi(6,:) = [day epoch 4 4];
% cellsi(7,:) = [day epoch 4 3];
% %%% cellsi(8,:) = [day epoch 14 1]; %iHp % Silent
% cellsi(8,:) = [day epoch 1 6];
% cellsi(9,:) = [day epoch 1 4]; % Can be in or out
% cellsi(11,:) = [day epoch 4 1]; 
% cellsi(10,:) = [day epoch 7 1];
% cellsi(12,:) = [day epoch 1 2];
% usecellsi=1:12;

% % Ordered - For 7758.6-7722 (13 0r 14 cells)
%cellsi(1,:) = [day epoch 1 3]; % Too many
cellsi(1,:) = [day epoch 12 2]; % iHp
cellsi(2,:) = [day epoch 8 1]; % iHp
cellsi(3,:) = [day epoch 1 1];
cellsi(4,:) = [day epoch 1 5];
cellsi(5,:) = [day epoch 4 5];
cellsi(6,:) = [day epoch 4 6];
cellsi(7,:) = [day epoch 4 4];
cellsi(8,:) = [day epoch 4 3];
%%% cellsi(8,:) = [day epoch 14 1]; %iHp % Silent
cellsi(9,:) = [day epoch 1 6];
cellsi(10,:) = [day epoch 1 4]; % Can be in or out
cellsi(12,:) = [day epoch 4 1]; 
cellsi(11,:) = [day epoch 7 1];
cellsi(13,:) = [day epoch 1 2];
usecellsi=1:13;



% % Ordered - For 7860-7970
% %cellsi(1,:) = [day epoch 1 3];
% cellsi(1,:) = [day epoch 12 2]; % iHp
% %cellsi(3,:) = [day epoch 8 1]; % iHp
% cellsi(2,:) = [day epoch 1 1];
% cellsi(3,:) = [day epoch 1 5];
% cellsi(4,:) = [day epoch 4 5];
% cellsi(5,:) = [day epoch 4 6];
% cellsi(6,:) = [day epoch 4 4];
% cellsi(7,:) = [day epoch 4 3];
% cellsi(8,:) = [day epoch 14 1]; %iHp
% cellsi(9,:) = [day epoch 1 6];
% %cellsi(12,:) = [day epoch 1 4];
% cellsi(10,:) = [day epoch 4 1];
% cellsi(11,:) = [day epoch 7 1];
% cellsi(12,:) = [day epoch 1 2];
% usecellsi=1:12;


% PFC
cellsp(1,:) = [day epoch 15 1];
cellsp(2,:) = [day epoch 15 2];
cellsp(3,:) = [day epoch 15 3];
cellsp(4,:) = [day epoch 15 4];
cellsp(5,:) = [day epoch 17 1];
cellsp(6,:) = [day epoch 17 2];
cellsp(7,:) = [day epoch 18 1];
cellsp(8,:) = [day epoch 18 2];
cellsp(9,:) = [day epoch 18 3];
cellsp(10,:) = [day epoch 18 4];

usecellsp=1:10;


% Spike data
%-----------
spikefile = sprintf('%s/%sspikes%02d.mat', animdirect, prefix, day);
load(spikefile);
for i=1:size(cellsi,1)
    eval(['spiketimei{',num2str(i),'}= spikes{cellsi(',num2str(i),',1)}{cellsi(',num2str(i),',2)}'...
        '{cellsi(',num2str(i),',3)}{cellsi(',num2str(i),',4)}.data(:,1);']);
    eval(['spikeposi{',num2str(i),'}= spikes{cellsi(',num2str(i),',1)}{cellsi(',num2str(i),',2)}'...
        '{cellsi(',num2str(i),',3)}{cellsi(',num2str(i),',4)}.data(:,2:3);']);
    eval(['spikeposidxi{',num2str(i),'}= spikes{cellsi(',num2str(i),',1)}{cellsi(',num2str(i),',2)}'...
        '{cellsi(',num2str(i),',3)}{cellsi(',num2str(i),',4)}.data(:,7);']);
end

for i=1:size(cellsp,1)
    eval(['spiketimep{',num2str(i),'}= spikes{cellsp(',num2str(i),',1)}{cellsp(',num2str(i),',2)}'...
        '{cellsp(',num2str(i),',3)}{cellsp(',num2str(i),',4)}.data(:,1);']);
    eval(['spikeposp{',num2str(i),'}= spikes{cellsp(',num2str(i),',1)}{cellsp(',num2str(i),',2)}'...
        '{cellsp(',num2str(i),',3)}{cellsp(',num2str(i),',4)}.data(:,2:3);']);
    eval(['spikeposidxp{',num2str(i),'}= spikes{cellsp(',num2str(i),',1)}{cellsp(',num2str(i),',2)}'...
        '{cellsp(',num2str(i),',3)}{cellsp(',num2str(i),',4)}.data(:,7);']);
end

% Extracted Ripples
%--------------
% Load ripple file
ripfile = sprintf('%s/%sripplall%02d.mat', directoryname, prefix, day);
load(ripfile);

% Get ripstarttime or midtime for current tetrodes in given day and epoch
rip_starttime = ripplesall{day}{epoch}{2}.starttime;   % in sec - Ntet=2 condition
ripsize = ripplesall{day}{epoch}{2}.maxthresh;
pt = rip_starttime;   % All ripple start time
discard = find(ripsize<=3);
pt(discard)=[];


%eind = lookup(pt(:,1), ti);

% EEg and Ripple cont data
% ------------------------
tets=[riptetlist,peegtet];
for te=1:length(tets)
    
    currtet=tets(te);
    cnteeg=0; cntrip=0;
    
    % Load EEG and ripple LFP file
    %-------------------------
    EEGfile = sprintf('%s/EEG/%seeggnd%02d-%01d-%02d.mat', directoryname, prefix, day,epoch,currtet);
    load(EEGfile);
    eeg=eeggnd;
    e = eeg{day}{epoch}{currtet};
    if te==1
        teeg = geteegtimes(e);
        eind = lookup(pt, teeg);
        e.samprate=round(e.samprate);
        eegstart = eeg{day}{epoch}{te}.starttime; % in secs - Epoch start
        eegend = teeg(end); % in secs - Epoch end
    end
    ripfile = sprintf('%s/EEG/%sripple%02d-%01d-%02d.mat', directoryname, prefix, day,epoch,currtet);
    load(ripfile);
    ripamp = ripple{day}{epoch}{currtet}.data(:,1);
    ripenv = ripple{day}{epoch}{currtet}.data(:,3);
    
    eegalltet{te} = e.data;
    
   
end % end tets


% Position Data
% --------------
linposfile = sprintf('%s/%slinpos%02d.mat', animdirect, prefix, day);
load(linposfile);
statematrix = linpos{day}{epoch}.statematrix;
linpostime = statematrix.time;

posfile = sprintf('%s/%spos%02d.mat', animdirect, prefix, day);
load(posfile);
absvel = abs(pos{day}{epoch}.data(:,5)); % Can also use field 9
posn = pos{day}{epoch}.data(:,2:3); % Can also use fields 6:7
postime = pos{day}{epoch}.data(:,1); % same as linpostime

% Get track, wells and traj info
trajwells = linpos{day}{epoch}.trajwells;
wellCoord = linpos{day}{epoch}.wellSegmentInfo.wellCoord;



% --------------------------------------------------------------------------------------------------
%% Plot example single-trial LFPs aligned to ripples. This is from sj_ripegs1. Plots LFP on all tets
% Here, use only main eggtet and plot EEG and Ripple. THen plot spike raster of all cells below
% -------------------------------------------------------------------------------------------------


% ------------------------------
% Figure and Font Sizes
forppr = 0;
% If yes, everything set to redimscreen_figforppr1
% If not, everything set to redimscreen_figforppt1

set(0,'defaultaxesfontweight','normal'); set(0,'defaultaxeslinewidth',2);

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
clr = {'b','m','g','y','c','k','r','b','g','y','b','m','g','y','c','k','r','b','g','y','b','m','g','y','c','k','r','b','g','y'};

clr1='k';
clr2='r';

datadir = '/data25/sjadhav/HPExpt/HPa_direct/ProcessedData/';
figdir = '/data25/sjadhav/HPExpt/HPa_direct/Figures/';
% ---------------------------------------

% winst = eegstart; winend = winst + win; % secs
winst = 7658.6; winend = 7722;
%winst = 7650; winend = 7950;
%winst = 7860; winend = 7970;
%winst = 7310; winend = 7440;
while winend <= eegend
    
    
    % For Manual
   
    
    % Only do if atleast one ripple is bigger than threshold sd
    %if max(ripsize_stim(i,:)) >= sd
    figure(n); hold on;
    %redimscreen_halfvert;
    redimscreen;
    
    % Get eeg times axis
    eind1 = lookup(winst, teeg);
    eind2 = lookup(winend, teeg);
    taxis = teeg(eind1:eind2);
    taxis = taxis - winst;
    
    winst_ms = winst*1000;
    winend_ms = winend*1000;
    
    %currriptime = pt(i)*1000;  % In ms Or use eeg time - currind = eind(i);
    %sttime = currriptime - pret; sttime_sec = sttime./1000;
    %endtime = currriptime + postt; endtime_sec = endtime./1000;
    
    % Look up position index at this time in stateatrix
    % posidx_rip = lookup(pt(i),statematrix.time); % All in secs
    posidx_st = lookup(winst,statematrix.time);
    posidx_end = lookup(winend,statematrix.time);
    currpostime = statematrix.time(posidx_st:posidx_end);
    currpostime = currpostime - currpostime(1); % Make time axis start at 0
    
    % Find ripples within this window
    ripidx = find(pt>=winst & pt<=winend);
    riptimes = pt(ripidx); % In secs
    riptimes_win = riptimes - winst;
    
    % First Spikes on bottom. Each tick has space of height2, and using 1.8 of it for line
    % -----------------------
    % PFC - SKIP PFC FOR NOW
    % ----
    cnt = 0; baseline = 0;
    
    % Now, CA1
    % --------
    
    %Divider line
    % If no PFC, No Need of the following
    % ----------------------------------
    %cnt = 0; baseline = baseline + size(cellsp,1)*2;
    %xpts = -pret:1:postt;
    %ypts = baseline*ones(size(xpts));
    %plot(xpts , ypts, 'k--','Linewidth',2);
    
    for c=usecellsi
        eval(['currspkt = spiketimei{',num2str(c),'};']);
        currspkt = currspkt;
        currspkt = currspkt(find(currspkt>=winst & currspkt<=winend ));
        
        % If spikes, subtract from subtract from start time and bin
        if ~isempty(currspkt)
            currspkt = currspkt - winst;
            %raster = starttime:0.001:endtime;
        end
        
        cnt=cnt+1;
        figure(n); hold on; %subplot(nplots,1,cnt+5); hold on;
        if ~isempty(currspkt)
            if size(currspkt,2)~=1, currspkt=currspkt'; end
            % Use plotraster or spikeTrain
            %plotraster(currspkt,(baseline+2*(cnt-1))*ones(size(currspkt)),1.8,[],'Color',[clr{cnt}],'LineWidth',2);
            plotraster(currspkt,(baseline+2*(cnt-1))*ones(size(currspkt)),1.8,[],'Color',[clr1],'LineWidth',2);
            %spikeTrain(currspkt,(baseline+(c-1))*ones(size(currspkt)),0.8,[],'Color',[clr{cnt}]);
        else
            plot(0,0,'k.');
        end
        %set(gca,'XLim',[0 endtime-starttime]);
    end
    
    Hp_baseline = baseline; % Save for linearized position
    % Update baseline
    baseline = baseline + size(cellsi,1)*2;
    up = baseline;
    %Divider line
    xpts = 0:1:win;
    ypts = baseline*ones(size(xpts));
    %plot(xpts , ypts, 'k--','Linewidth',2);
    
    % Linearized position (from well 1 and well 3) on top of Hp spikes
    plotscale = up-Hp_baseline;
    lindist1 = statematrix.linearDistanceToWells(posidx_st:posidx_end,1);
    lindist1 = lindist1-min(lindist1);
    distscale = max(lindist1)-min(lindist1);
    lindist1 = Hp_baseline + lindist1.*(plotscale/distscale);
    plot(currpostime,lindist1,'r-','LineWidth',4);
%     lindist3 = statematrix.linearDistanceToWells(posidx_st:posidx_end,3);
%     lindist3 = lindist3-min(lindist3);
%     distscale = max(lindist3)-min(lindist3);
%     lindist3 = Hp_baseline + lindist3.*(plotscale/distscale);
%     plot(currpostime,lindist3,'b-','LineWidth',2);
    
    
%     % EEG On main tet
%     % --------------
%     n = maineegidx;
%     eegtet = eegalltet{n};
%     curreeg = eegtet(eind1:eind2);
%     % Plot
%     % ----
%     eegscale = max(curreeg)-min(curreeg);
%     downeeg = baseline+0.5; upeeg = downeeg+4;
%     plotscale = 4;
%     curreeg = downeeg + (plotscale/2) + curreeg.*(plotscale/eegscale);
%     plot(taxis,curreeg,'k-','LineWidth',1);
%     
%     % Update baseline
%     baseline = baseline + 4;
%     
%     % EEG on PFC Tet instead of ripple band
%     % --------------------------------------
%     n = peegidx;
%     eegtet = eegalltet{n};
%     curreeg = eegtet(eind1:eind2);
%     % Plot
%     % ----
%     eegscale = max(curreeg)-min(curreeg);
%     downeeg = baseline+0.5; upeeg = downeeg+4;
%     plotscale = 4;
%     curreeg = downeeg + (plotscale/2) + curreeg.*(plotscale/eegscale);
%     plot(taxis,curreeg,'r-','LineWidth',1);
    
    
    %         % Rippleband
    %         % ------------
    %         ripscale = max(currrip)-min(currrip);
    %         downrip = upeeg;
    %         uprip = downrip + 8;
    %         plotscale = 8;
    %         currrip = downrip + (plotscale/2) + currrip.*(plotscale/eegscale);
    %         plot(taxis,currrip,'k-','LineWidth',1);
    %
    %         % Make a line at 0ms (start of ripple) and 100ms
    %         ypts = 0:1:uprip+1;
    %         xpts = 0*ones(size(ypts)); % Plot Line at 0ms
    %         plot(xpts , ypts, 'r:','Linewidth',2);
    %         xpts = lineat*ones(size(ypts)); % Plot lines at 100ms
    %         plot(xpts , ypts, 'r:','Linewidth',2);
    %
    %         set(gca, 'YLim',[0 uprip]);
    %         set(gca, 'XLim',[-pret-10 postt+10]);
    
    set(gca, 'YTick',[]);
    
    winsecs = [0:10:winend-winst]; %eg. 0:10:60
    secs = [winst:10:winend];
    secs = round(secs); %secs = roundn(secs,-1);
    msecs = [winst_ms:10000:winend_ms];
    set(gca,'XTick',winsecs,'XTickLabel',num2str(secs'));
    xlabel('Time (secs)','FontSize',18,'Fontweight','normal');
    title(['Win St Ttime: ' num2str(roundn(winst,-1))],'FontSize',18,'Fontweight','normal');
    %title(['SWR Time: ',num2str(roundn(pt(i),-1)), 'secs']);
    
    % Draw Lines
    ylim = get(gca,'YLim');
    ypts = ylim(1):ylim(2);
    for i=1:length(winsecs)
        xpts = winsecs(i)*ones(size(ypts));
        %plot(xpts , ypts, 'k--','Linewidth',0.5);
    end
    
    % Mark Ripple Times in Current Window
%     for s=1:length(riptimes_win)
%         DIOt = riptimes_win(s);
%         xaxis = DIOt:0.05:DIOt+0.3; % 300ms
%         jbfill(xaxis,ylim(2)*ones(size(xaxis)),ylim(1)*ones(size(xaxis)),'m','m',1,0.2);
%     end
        
    
    
    set(gca,'YLim',[0 size(cellsi,1)*2]);
    set(gca,'XLim',[0 winend-winst]);
    
    set(gcf,'Position',[25 700 1280 400]);
    
    saveg=0;
    if saveg==1
        figfile = [figdir,prefix,'_Day2Ep4_LongRasterEg1c_13cells'];
        print('-dpdf', figfile); print('-djpeg', figfile); saveas(gcf,figfile,'fig');
    end
    keyboard; % Pause after each plot
    
    close all
    
    % Update winst and winend
    winst = winst + win - overlap;
    winend = winst + win;
    
    
    %end % if ripsize > sd
end





% if saveg1==1
%     figfile = [figdir,'n',prefix,'_Day',num2str(day),'Ep',num2str(epoch),'Stim',num2str(i),'_RipEg1'];
%     print('-dpdf', figfile);
%     print('-djpeg', figfile);
%     saveas(gcf,figfile,'fig');
% end

% % Manually plot
% tet = 1; stim=75;
% figure; hold on;
% eegtet=e_stim{tet};
% x=eegtet(stim,:); % eg 75 tet4
% plot(taxis,x,'k-','LineWidth',2);
% set(gca,'XLim',[-50 50]);
% figfile = [figdir,prefix,'_Day',num2str(day),'Ep',num2str(epoch),'Tet',num2str(tet),'Stim',num2str(stim),'_lfpeg1'];
% print('-dpdf', figfile);
% print('-djpeg', figfile);
% saveas(gcf,figfile,'fig');
% 
% figure; hold on;
% riptet=ripamp_stim{tet};
% y=riptet(stim,:); % eg 75 tet4
% plot(taxis,y,'k-','LineWidth',2);
% set(gca,'XLim',[-50 50]);
% 
% figfile = [figdir,prefix,'_Day',num2str(day),'Ep',num2str(epoch),'Tet',num2str(tet),'Stim',num2str(stim),'_ripeg1'];
% print('-dpdf', figfile);
% print('-djpeg', figfile);
% saveas(gcf,figfile,'fig');

% *************************************************************************
cd(datadir);
cd(animdirect);
keyboard;

