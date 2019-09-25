
% Shantanu - Update to get ripples from getripplesdirectg with spiking condition. 
% So these are candidate replay events

%function [pret] = plotsj_rip_spikeraster(prefix, day, epoch, saveg1)
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

clear;
prefix='HPa';
%day=1; epoch=4;
day=2; epoch=4;
sd=4; %% SD threshold for ripples
% Pret and Post from start of ripple? Also plot start-middle and end of ripple
pret=1000;
postt=1000;




switch prefix
    case 'HPa'
        directoryname = '/data25/sjadhav/HPExpt/HPa_direct';
        animdirect = directoryname;
        dire = '/data25/sjadhav/HPExpt/HPa';
        riptetlist = [1,4,5,6,7,8,9,11,12,14];
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
peegtet = 16;
eegtets = [eegtets, peegtet];
peegidx=length(eegtets);

% HPa - day 2. CA1 cells 6+4; PFC cells 4+4
% CA1
cellsi(1,:) = [day epoch 1 1];  
cellsi(2,:) = [day epoch 1 2];  
cellsi(3,:) = [day epoch 1 3];  
cellsi(4,:) = [day epoch 1 4];  
cellsi(5,:) = [day epoch 1 5];
cellsi(6,:) = [day epoch 1 6]; 
cellsi(7,:) = [day epoch 4 1]; 
cellsi(8,:) = [day epoch 4 2];
cellsi(9,:) = [day epoch 4 3]; 
cellsi(10,:) = [day epoch 4 4];
cellsi(11,:) = [day epoch 4 5];
cellsi(12,:) = [day epoch 4 6];
cellsi(13,:) = [day epoch 7 1];
cellsi(14,:) = [day epoch 8 1];
cellsi(15,:) = [day epoch 12 1];
cellsi(16,:) = [day epoch 12 2];
cellsi(17,:) = [day epoch 14 1];
cellsi(18,:) = [day epoch 14 2];
cellsi(19,:) = [day epoch 14 3];

usecellsi=1:19;

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

% % HPa - day 1. 
% % CA1
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
% cellsi(13,:) = [day epoch 4 7];
% cellsi(14,:) = [day epoch 4 9];
% cellsi(15,:) = [day epoch 4 10];
% cellsi(16,:) = [day epoch 4 11];
% cellsi(17,:) = [day epoch 4 12];
% cellsi(18,:) = [day epoch 4 13];
% 
% usecellsi=1:18;
% 
% % PFC
% 
% cellsp(1,:) = [day epoch 16 2];
% cellsp(2,:) = [day epoch 15 1];  
% cellsp(3,:) = [day epoch 15 2];  
% cellsp(4,:) = [day epoch 16 1];  
% cellsp(5,:) = [day epoch 17 2];
% cellsp(6,:) = [day epoch 18 2];
% cellsp(7,:) = [day epoch 18 3]; 
% 
% usecellsp=1:7;



% Spike data
%-----------
spikefile = sprintf('%s/%sspikes%02d.mat', animdirect, prefix, day);
load(spikefile);
tetinfofile = sprintf('%s/%stetinfo.mat', directoryname, prefix);
load(tetinfofile);
cellinfofile = sprintf('%s/%scellinfo.mat', directoryname, prefix);
load(cellinfofile);

for i=1:size(cellsi,1)
    i;
    eval(['spiketimei{',num2str(i),'}= spikes{cellsi(',num2str(i),',1)}{cellsi(',num2str(i),',2)}'...
        '{cellsi(',num2str(i),',3)}{cellsi(',num2str(i),',4)}.data(:,1);']);
    eval(['spikeposi{',num2str(i),'}= spikes{cellsi(',num2str(i),',1)}{cellsi(',num2str(i),',2)}'...
        '{cellsi(',num2str(i),',3)}{cellsi(',num2str(i),',4)}.data(:,2:3);']);
    eval(['spikeposidxi{',num2str(i),'}= spikes{cellsi(',num2str(i),',1)}{cellsi(',num2str(i),',2)}'...
        '{cellsi(',num2str(i),',3)}{cellsi(',num2str(i),',4)}.data(:,7);']);
end

for i=1:size(cellsp,1)
    i;
    eval(['spiketimep{',num2str(i),'}= spikes{cellsp(',num2str(i),',1)}{cellsp(',num2str(i),',2)}'...
        '{cellsp(',num2str(i),',3)}{cellsp(',num2str(i),',4)}.data(:,1);']);
    eval(['spikeposp{',num2str(i),'}= spikes{cellsp(',num2str(i),',1)}{cellsp(',num2str(i),',2)}'...
        '{cellsp(',num2str(i),',3)}{cellsp(',num2str(i),',4)}.data(:,2:3);']);
    eval(['spikeposidxp{',num2str(i),'}= spikes{cellsp(',num2str(i),',1)}{cellsp(',num2str(i),',2)}'...
        '{cellsp(',num2str(i),',3)}{cellsp(',num2str(i),',4)}.data(:,7);']);
end

% Extracted Ripples
%--------------
% OLD METH
% ripfile = sprintf('%s/%sripplall%02d.mat', directoryname, prefix, day);
% load(ripfile);
% % Get ripstarttime or midtime for current tetrodes in given day and epoch
% %rip_starttime = ripples{day}{epoch}{maineegtet}.starttime;   % in sec - Ntet=2 condition
% %ripsize = ripples{day}{epoch}{maineegtet}.maxthresh;
% rip_starttime = ripplesall{day}{epoch}{2}.starttime;   % in sec - Ntet=2 condition
% ripsize = ripplesall{day}{epoch}{2}.maxthresh;
% pt = rip_starttime;   % All ripple start time
% discard = find(ripsize<=9);
% pt(discard)=[];

% Ripple times across tetrode
% ---------------------------------
cellcountthresh=5;
ripfile = sprintf('%s/%sripples%02d.mat', directoryname, prefix, day);
load(ripfile);
riptets = riptetlist;
riptimes = []; rip_starttime = []; triggers = [];
[riptimes] = getripples_direct([day, epoch], ripples, riptets,'minstd',3);
currriptet=riptets(1);
% SPIKE COUNT THRESHOLD - MOVE TO ANOTHER FILE?
    %filterString = '( strcmp($tag, ''CA1Pyr'') || strcmp($tag, ''iCA1Pyr'') || strcmp($tag, ''CA1Run'') || strcmp($tag, ''CA1Pyrpr'') || strcmp($tag, ''CA1Pyrp'') || strcmp($tag, ''iCA1Run'') || strcmp($tag, ''iCA1Pyrpr'') || strcmp($tag, ''iCA1Pyrp'') )';
    filterString = ' (   strcmp($tag2, ''CA1Pyr'') || strcmp($tag2, ''iCA1Pyr'')    )';
    cellindices = evaluatefilter(cellinfo{day}{epoch}, filterString);
    indices = [repmat([day epoch], size(cellindices,1),1 ), cellindices];
    spikecounts = [];   celldata = [];
    %go through each cell and calculate the binned spike counts
    for cellcount = 1:size(indices,1)
        index = indices(cellcount,:);
        if ~isempty(spikes{index(1)}{index(2)}{index(3)}{index(4)}.data)
            spiketimes = spikes{index(1)}{index(2)}{index(3)}{index(4)}.data(:,1);
        else
            spiketimes = [];
        end
        spikebins = periodAssign(spiketimes, riptimes(:,[1 2]));
        if ~isempty(spiketimes)
            validspikes = find(spikebins);
            spiketimes = spiketimes(validspikes);
            spikebins = spikebins(validspikes);
        end
        
        if ~isempty(spiketimes)
            tmpcelldata = [spiketimes spikebins];
            tmpcelldata(:,3) = cellcount;
        else
            tmpcelldata = [0 0 cellcount];
        end
        celldata = [celldata; tmpcelldata];
        spikecount = zeros(1,size(riptimes,1));
        for i = 1:length(spikebins)
            spikecount(spikebins(i)) = spikecount(spikebins(i))+1;
        end
        
        spikecounts = [spikecounts; spikecount];
    end
    celldata = sortrows(celldata,1); %sort all spikes by time
    cellcounts = sum((spikecounts > 0));
    %Find all events with enough cells
    eventindex = find(cellcounts >= cellcountthresh);
    riptimes_keep = riptimes(eventindex,:);
    rip_starttime = 1000*riptimes_keep(:,1);  % in ms
    rip_endtime = 1000*riptimes_keep(:,2);
    % Find ripples separated by atleast a second
    % --------------------------------------------
    [rip_starttime,sortidx] = sort(rip_starttime);
    triggers = rip_starttime; triggers_end = rip_endtime;
    iri = diff(triggers);
    keepidx = [1;find(iri>=1000)+1];
    triggers = triggers(keepidx); triggers_end = triggers_end(keepidx);
    pt = triggers;
    pt_end = triggers_end;


%eind = lookup(pt(:,1), ti);

% EEg and Ripple cont data
% ------------------------
eegpt = pt./1000; % in secs
tets=[riptetlist,peegtet];
for te=1:length(tets)
    
    currtet=tets(te);
    cnteeg=0; cntrip=0;
    
    % Load EEG and ripple LFP file
    %-------------------------
    EEGfile = sprintf('%s/EEG/%seeg%02d-%01d-%02d.mat', directoryname, prefix, day,epoch,currtet);
    load(EEGfile);
    e = eeg{day}{epoch}{currtet};
    if te==1
        teeg = geteegtimes(e);
        eind = lookup(eegpt, teeg);
        e.samprate=round(e.samprate);
        eegstart = eeg{day}{epoch}{te}.starttime; % Not used?
    end
    ripfile = sprintf('%s/EEG/%sripple%02d-%01d-%02d.mat', directoryname, prefix, day,epoch,currtet);
    load(ripfile);
    ripamp = ripple{day}{epoch}{currtet}.data(:,1);
    ripenv = ripple{day}{epoch}{currtet}.data(:,3);
    
    
             
    % Align EEG and Ripple Band to ripple time. 
    %------------------------------------------

    nelements = length(1000-round((pret/1000)*e.samprate):1000+round((postt/1000)*e.samprate));
    for i=1:length(eegpt) % Need to Skip initial and final indices?
        i;
        cnteeg=cnteeg+1;
        currriptime = eegpt(i); %currripsize = ripsize(i);
        currind = eind(i);
        if ( (currind-round((pret/1000)*e.samprate) <=0) || (currind+round((postt/1000)*e.samprate)>length(e.data)) )
            e_stim{te}(cnteeg,:)=0*(1:nelements);
            ripamp_stim{te}(cnteeg,:)=0*(1:nelements);
            ripenv_stim{te}(cnteeg,:)=0*(1:nelements);
        else
            e_stim{te}(cnteeg,:)=e.data(currind-round((pret/1000)*e.samprate):currind+round((postt/1000)*e.samprate));
            ripamp_stim{te}(cnteeg,:)=double(ripamp(currind-round((pret/1000)*e.samprate):currind+round((postt/1000)*e.samprate)));
            ripenv_stim{te}(cnteeg,:)=double(ripenv(currind-round((pret/1000)*e.samprate):currind+round((postt/1000)*e.samprate)));
            %ripsize_stim(cnteeg,te) = currripsize; 
        end
    end
    
end % end tets


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
 
datadir = '/data25/sjadhav/HPExpt/HPa_direct/ProcessedData/';
figdir = '/data25/sjadhav/HPExpt/HPa_direct/Figures/';
% ---------------------------------------

%nplots=length(tets); Use only main eegtet

taxis = [1:size(ripenv_stim{1},2)]*1000/e.samprate;
taxis = taxis-pret;
lineat=200;
    
for i=1:size(e_stim{1},1)
    
    % Only do if atleast one ripple is bigger than threshold sd
    %if max(ripsize_stim(i,:)) >= sd
        figure(i); hold on;
        %redimscreen_halfvert;
        redimscreen;
        
        %for n=1:nplots
        
        currriptime = pt(i);  % In ms Or use eeg time - currind = eind(i); 
        currripend = pt_end(i);
        
        
        % Find all ripples within this window
        winst = currriptime-pret; % in ms
        winend = currriptime+postt; % in ms
        ripidx = find(pt>=(winst) & pt<=(winend));
        riptimec = pt(ripidx); % In ms
        riptimes_win = riptimec - currriptime; % Ripple time in current time axis
    
        % First Spikes on bottom. Each tick has space of height2, and using 1.8 of it for line
        % -----------------------
        % PFC
        % ----
        cnt = 0; baseline = 0;
        for c=usecellsp
            eval(['currspkt = spiketimep{',num2str(c),'};']);
            % Convert to ms
            currspkt = currspkt*1000;
            currspkt = currspkt(find(currspkt>=currriptime-pret & currspkt<=currriptime+postt ));
            
            % If spikes, subtract from subtract from start time and bin
            if ~isempty(currspkt)
                currspkt = currspkt - currriptime;
                %raster = starttime:0.001:endtime;
            end
                    
            cnt=cnt+1;
            figure(i); hold on; %subplot(nplots,1,cnt+5); hold on;
            if ~isempty(currspkt)
                if size(currspkt,2)~=1, currspkt=currspkt'; end
                % Use plotraster or spikeTrain
                plotraster(currspkt,(baseline+2*(cnt-1))*ones(size(currspkt)),1.8,[],'Color','b','LineWidth',2);
                %plotraster(currspkt,(baseline+2*(cnt-1))*ones(size(currspkt)),1.8,[],'Color',[clr{cnt}],'LineWidth',2);
                %spikeTrain(currspkt,(baseline+(c-1))*ones(size(currspkt)),0.8,[],'Color',[clr{cnt}]);
            else
                plot(0,0,'k.');
            end
            %set(gca,'XLim',[0 endtime-starttime]);
        end
        
        % Now, CA1
        % --------
        cnt = 0; baseline = baseline + size(cellsp,1)*2;
        %Divider line
        xpts = -pret:1:postt;
        ypts = baseline*ones(size(xpts)); 
        plot(xpts , ypts, 'k--','Linewidth',2);
        
        for c=usecellsi
            eval(['currspkt = spiketimei{',num2str(c),'};']);
            currspkt = currspkt *1000;
            currspkt = currspkt(find(currspkt>=currriptime-pret & currspkt<=currriptime+postt ));
            
            % If spikes, subtract from subtract from start time and bin
            if ~isempty(currspkt)
                currspkt = currspkt - currriptime;
                %raster = starttime:0.001:endtime;
            end
            
            cnt=cnt+1;
            figure(i); hold on; %subplot(nplots,1,cnt+5); hold on;
            if ~isempty(currspkt)
                if size(currspkt,2)~=1, currspkt=currspkt'; end
                % Use plotraster or spikeTrain
                plotraster(currspkt,(baseline+2*(cnt-1))*ones(size(currspkt)),1.8,[],'Color','k','LineWidth',2);
                %plotraster(currspkt,(baseline+2*(cnt-1))*ones(size(currspkt)),1.8,[],'Color',[clr{cnt}],'LineWidth',2);
                %spikeTrain(currspkt,(baseline+(c-1))*ones(size(currspkt)),0.8,[],'Color',[clr{cnt}]);
            else
                plot(0,0,'k.');
            end
            %set(gca,'XLim',[0 endtime-starttime]);
        end
        
        % Update baseline
        baseline = baseline + size(cellsi,1)*2;
         %Divider line
        xpts = -pret:1:postt;
        ypts = baseline*ones(size(xpts)); 
        plot(xpts , ypts, 'k--','Linewidth',2);
        
        % EEG On main tet
        % --------------
        n = maineegidx;
        eegtet = e_stim{n};
        riptet = ripamp_stim{n};
        ripenvtet = ripenv_stim{n};
        curreeg = eegtet(i,:);
        currrip = riptet(i,:);
        %currsize = roundn(ripsize_stim(i,n),-2);
       
        % EEg
        % ----
        eegscale = max(curreeg)-min(curreeg);
        downeeg = baseline+2; upeeg = downeeg+8;
        plotscale = 8;
        curreeg = downeeg + (plotscale/2) + curreeg.*(plotscale/eegscale);   
        plot(taxis,curreeg,'k-','LineWidth',1);
        
         % Update baseline
        baseline = baseline + 8;
        
        % EEG on PFC Tet instead of ripple band
        % --------------------------------------
        n = peegidx;
        eegtet = e_stim{n};
        riptet = ripamp_stim{n};
        ripenvtet = ripenv_stim{n};
        curreeg = eegtet(i,:);
        currrip = riptet(i,:);
        %currsize = roundn(ripsize_stim(i,n),-2);
        eegscale = max(curreeg)-min(curreeg);
        downeeg = baseline+2; upeeg = downeeg+8;
        plotscale = 8;
        curreeg = downeeg + (plotscale/2) + curreeg.*(plotscale/eegscale);   
        plot(taxis,curreeg,'b-','LineWidth',1);
        
         % Update baseline
        baseline = baseline + 8;
        
        % Mark ripple start and end
        % --------------------------
        ypts = 0:1:baseline+1;
        xpts = 0*ones(size(ypts)); % Plot Line at 0ms
        plot(xpts , ypts, 'r:','Linewidth',2);
        lineat = currripend-currriptime;
        xpts = lineat*ones(size(ypts)); % Plot lines at 100ms
        plot(xpts , ypts, 'r:','Linewidth',2);
        
        
        
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
       

%         ylim = get(gca,'YLim');
%         ypts = ylim(1):ylim(2);
%         % Mark Ripple Times in Current Window
%         for s=1:length(riptimes_win)
%             DIOt = riptimes_win(s);
%             xaxis = DIOt:50:DIOt+300; % 300ms
%             jbfill(xaxis,ylim(2)*ones(size(xaxis)),ylim(1)*ones(size(xaxis)),'m','m',1,0.2);
%         end

        
        
        
        title(['SWR aligned'],'FontSize',18,'Fontweight','normal');
        xlabel('Time (ms)','FontSize',18,'Fontweight','normal');
        set(gca, 'YTick',[]);
               
        saveg=1;
        if saveg==1
            figfile = [figdir,prefix,'Day2Ep5RasterEEGegNo',num2str(i)]
            print('-dpdf', figfile); print('-djpeg', figfile); saveas(gcf,figfile,'fig');
        end
        keyboard; % Pause after each plot
        
        close all
        
    %end % if ripsize > sd
end





if saveg1==1
    figfile = [figdir,'n',prefix,'_Day',num2str(day),'Ep',num2str(epoch),'Stim',num2str(i),'_RipEg1'];
    print('-dpdf', figfile);
    print('-djpeg', figfile);
    saveas(gcf,figfile,'fig');
end
  
% Manually plot
tet = 1; stim=75;
figure; hold on; 
eegtet=e_stim{tet};
x=eegtet(stim,:); % eg 75 tet4
plot(taxis,x,'k-','LineWidth',2);
set(gca,'XLim',[-50 50]);
figfile = [figdir,prefix,'_Day',num2str(day),'Ep',num2str(epoch),'Tet',num2str(tet),'Stim',num2str(stim),'_lfpeg1'];
print('-dpdf', figfile);
print('-djpeg', figfile);
saveas(gcf,figfile,'fig');

figure; hold on; 
riptet=ripamp_stim{tet};
y=riptet(stim,:); % eg 75 tet4
plot(taxis,y,'k-','LineWidth',2);
set(gca,'XLim',[-50 50]);

figfile = [figdir,prefix,'_Day',num2str(day),'Ep',num2str(epoch),'Tet',num2str(tet),'Stim',num2str(stim),'_ripeg1'];
print('-dpdf', figfile);
print('-djpeg', figfile);
saveas(gcf,figfile,'fig');

% *************************************************************************
cd(datadir);
keyboard;

