
% MAke this a generic version so that it can plot any day-epoch. Load Ca1 and PFC cells automatically. 
% Can limit to ripple responsive PFC cells. Can plot position on top. 
% Finally, can also arrange CA1 cells by linear field alignment, but need to knwo which trajectory the animal
% is on during the ripple.

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
day=1; epoch=4;
sd=3; %% SD threshold for ripples - Dont need this forgetripples. Using 3 sds 
% Pret and Post from start of ripple? Also plot start-middle and end of ripple

pret=15000; postt=2000;
%pret=60000; postt=60000;


switch prefix
    case 'HPa'
        directoryname = '/data25/sjadhav/HPExpt/HPa_direct';
        animdirect = directoryname;
        dire = '/data25/sjadhav/HPExpt/HPa';
        riptetlist = [1,4,5,6,7,8,9,11,12,14];
        maineegtet = 1; maineegidx=1; % CA1 tet
        peegtet = 16; % PFCtet
    case 'HPb'
        directoryname = '/data25/sjadhav/HPExpt/HPb_direct';
        animdirect = directoryname;
        dire = '/data25/sjadhav/HPExpt/HPb';
        riptetlist = [1,3,4,5,6,16,17,18,20];
        maineegtet = 1; maineegidx=1; % CA1 tet
        peegtet = 9; % PFCtet
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
% Also get a PFC eeg
eegtets = [eegtets, peegtet];
peegidx= find(eegtets==peegtet);

% Get data files
% Spike data
%-----------
spikefile = sprintf('%s/%sspikes%02d.mat', animdirect, prefix, day);
load(spikefile);
tetinfofile = sprintf('%s/%stetinfo.mat', directoryname, prefix);
load(tetinfofile);
cellinfofile = sprintf('%s/%scellinfo.mat', directoryname, prefix);
load(cellinfofile);
ripfile = sprintf('%s/%sripples%02d.mat', directoryname, prefix, day);
load(ripfile);
posfile = sprintf('%s/%spos%02d.mat', animdirect, prefix, day);
load(posfile);
linposfile = sprintf('%s/%slinpos%02d.mat', animdirect, prefix, day);
load(linposfile);


% Pos - Linpos
% -------------
absvel = abs(pos{day}{epoch}.data(:,5)); % Can also use field 9
posn = pos{day}{epoch}.data(:,2:3); % Can also use fields 6:7
postime = pos{day}{epoch}.data(:,1); % same as linpostime
statematrix = linpos{day}{epoch}.statematrix;
linpostime = statematrix.time;
% Get track, wells and traj info
trajwells = linpos{day}{epoch}.trajwells;
wellCoord = linpos{day}{epoch}.wellSegmentInfo.wellCoord;


% Get cells
% ---------
% CA1 cells
%filterString = '( ($meanrate > 0.1) && ~strcmp($tag, ''iCA1Int'') && ~strcmp($tag, ''CA1Int'') && (strcmp($tag2, ''CA1Pyr'') || strcmp($tag2, ''iCA1Pyr'')) ) ';
filterString = '( (strcmp($tag2, ''CA1Pyr'') || strcmp($tag2, ''iCA1Pyr'')) && ($meanrate > 0.1) && ~strcmp($tag, ''iCA1Int'') && ~strcmp($tag, ''CA1Int'') ) ';
cellindices = evaluatefilter(cellinfo{day}{epoch}, filterString);
cellsi = [repmat([day epoch], size(cellindices,1),1 ), cellindices]; % day-epoch-tet-cell for CA1 cells
usecellsi = 1:size(cellsi,1);

% PFC cells
%filterString = 'strcmp($area, ''PFC'') && strcmp($ripmodtag, ''y'')';
filterString = 'strcmp($area, ''PFC'') && ($meanrate > 0.1)';
pcellindices = evaluatefilter(cellinfo{day}{epoch}, filterString);
cellsp = [repmat([day epoch], size(pcellindices,1),1 ), pcellindices]; % day-epoch-tet-cell for CA1 cells
usecellsp = 1:size(cellsp,1);



% Ripple times across tetrode
% ---------------------------------
cellcountthresh=8;

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
    rip_ncells = cellcounts(eventindex);
    riptimes_keep = riptimes(eventindex,:);
    rip_starttime = 1000*riptimes_keep(:,1);  % in ms
    rip_endtime = 1000*riptimes_keep(:,2);
    
    % Find ripples separated by atleast a second - SKIP THIS. Already looking across tetrodes and have cellcountthresh
    % ------------------------------------------------------------------------------------------------------------------
    [rip_starttime,sortidx] = sort(rip_starttime);
    triggers = rip_starttime; triggers_end = rip_endtime;
    % SKIPPING Inter-Rip-Interval
%     iri = diff(triggers);
%     keepidx = [1;find(iri>=1000)+1];
%     triggers = triggers(keepidx); triggers_end = triggers_end(keepidx);
    pt = triggers;
    pt_end = triggers_end;

    nrip = length(pt),
    maxcells = max(rip_ncells)
    

%eind = lookup(pt(:,1), ti);


lastidx = size(cellsi,1); % For HPa d2 ep4, last cell not in linfields, so skip for now.
cellsi(lastidx,:)=[];
usecellsi = 1:size(cellsi,1); 


% SORT CA1 CELLS BY LINFIELDS As in plotsj_linearpositions... - HAVE TO CHOOSE TRAJ
% -------------------------------------------------

traj1=1; traj2=2; % Out Left and In Left

linfieldfile = sprintf('%s/%slinfields%02d.mat', animdirect, prefix, day);
load(linfieldfile);

for i=1:size(cellsi,1)    
    try
        currlf = linfields{day}{epoch}{cellsi(i,3)}{cellsi(i,4)}; %day, epoch, tet, cell  
    catch
        cellsi(i,3), cellsi(i,4),
        keyboard;
    end
    lf1 = currlf{traj1}(:,5); % Norm Occ Rate for traj 1.
    lf2 = currlf{traj2}(:,5);
    lfcomb = mean([lf1,lf2],2); % Combine In and Out. Dominant Trajectory will define position in order
    % Save trajs
    lin_traj1{i}=lf1;  lin_traj2{i}=lf2;  lin_trajcomb{i}=lfcomb; 
    
    % Get peak positions
    [peak_traj1(i), pos_traj1(i)] = max(lf1); 
    [peak_traj2(i), pos_traj2(i)] = max(lf2); 
    [peak_trajcomb(i), pos_trajcomb(i)] = max(lfcomb); 
end 
    
% Get Sort order
% --------------
[~,sort1] = sort(pos_traj1); % By traj1 - eg. 3=Out Right or 1=Out Left
[~,sort2] = sort(pos_traj2); % By traj2 - eg. 4=In Right or 2=In Left
[~,sortcomb] = sort(pos_trajcomb); % By trajcomb - eg. Combine Out and In Right/ Left

% Return ordered list
cellsi_sorttraj1 = cellsi(sort1',:);
cellsi_sorttraj2 = cellsi(sort2',:);
cellsi_sorttrajcomb = cellsi(sortcomb',:);

% Done ordering
% --------------

% Replace cellsi by sorted order
cellsi = cellsi_sorttrajcomb


% spiketimes and spikeposidxs

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





% EEg and Ripple cont data
% ------------------------
eegpt = pt./1000; % in secs
tets=[riptetlist,peegtet];
for te=1:length(tets)
    
    currtet=tets(te);
    cnteeg=0; cntrip=0;
    
    % Load EEG and ripple LFP file and theta LFP file
    %-------------------------
    EEGfile = sprintf('%s/EEG/%seeggnd%02d-%01d-%02d.mat', directoryname, prefix, day,epoch,currtet);
    load(EEGfile); eeg = eeggnd;
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
    
    thetafile = sprintf('%s/EEG/%stheta%02d-%01d-%02d.mat', directoryname, prefix, day,epoch,currtet);
    load(thetafile);
    thetaamp = theta{day}{epoch}{currtet}.data(:,1);
    thetaenv = theta{day}{epoch}{currtet}.data(:,3);
    % theta is downsampled = maybe time axis is different?
    etheta = theta{day}{epoch}{currtet};
    if te==1
        ttheta = geteegtimes(etheta);
        thetaind = lookup(eegpt, ttheta);
        etheta_samprate=round(etheta.samprate);
    end
    
             
    % Align EEG and Ripple Band to ripple time. 
    %------------------------------------------

    nelements = length(1000-round((pret/1000)*e.samprate):1000+round((postt/1000)*e.samprate));
    nelements_theta = length(1000-round((pret/1000)*etheta_samprate):1000+round((postt/1000)*etheta_samprate));
    for i=1:length(eegpt) % Need to Skip initial and final indices?
        i;
        cnteeg=cnteeg+1;
        currriptime = eegpt(i); %currripsize = ripsize(i);
        currind = eind(i);
        currthetaind = thetaind(i);
        if ( (currind-round((pret/1000)*e.samprate) <=0) || (currind+round((postt/1000)*e.samprate)>length(e.data)) )
            e_stim{te}(cnteeg,:)=0*(1:nelements);
            ripamp_stim{te}(cnteeg,:)=0*(1:nelements);
            ripenv_stim{te}(cnteeg,:)=0*(1:nelements);
            thetaamp_stim{te}(cnteeg,:)=0*(1:nelements_theta);
            thetaenv_stim{te}(cnteeg,:)=0*(1:nelements_theta);
        else
            e_stim{te}(cnteeg,:)=e.data(currind-round((pret/1000)*e.samprate):currind+round((postt/1000)*e.samprate));
            ripamp_stim{te}(cnteeg,:)=double(ripamp(currind-round((pret/1000)*e.samprate):currind+round((postt/1000)*e.samprate)));
            ripenv_stim{te}(cnteeg,:)=double(ripenv(currind-round((pret/1000)*e.samprate):currind+round((postt/1000)*e.samprate)));
            
            thetaamp_stim{te}(cnteeg,:)=double(thetaamp(currthetaind-round((pret/1000)*etheta_samprate):currthetaind+round((postt/1000)*etheta_samprate)));
            thetaenv_stim{te}(cnteeg,:)=double(thetaenv(currthetaind-round((pret/1000)*etheta_samprate):currthetaind+round((postt/1000)*etheta_samprate)));
            %ripsize_stim(cnteeg,te) = currripsize; 
        end
    end
    
end % end tets


% --------------------------------------------------------------------------------------------------
%% Plot example single-trial LFPs aligned to ripples. This is from sj_ripegs1. Plots LFP 
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
 
datadir = '/data25/sjadhav/HPExpt/ProcessedData/';
figdir = '/data25/sjadhav/HPExpt/Figures/Replay/';
% ---------------------------------------

%nplots=length(tets); Use only main eegtet

taxis = [1:size(ripenv_stim{1},2)]*1000/e.samprate;
taxis = taxis-pret;
thetaaxis = [1:size(thetaenv_stim{1},2)]*1000/etheta_samprate;
thetaaxis = thetaaxis-pret;
lineat=200;
    
for i=1:size(e_stim{1},1)
    
    
    if ( i ==7 || i == 11 || i ==12 || i ==18 || i ==23)
    
    % Only do if atleast one ripple is bigger than threshold sd
    %if max(ripsize_stim(i,:)) >= sd
      
        
        %for n=1:nplots
        
        % Ripple time
        %------------
        currriptime = pt(i);  % In ms Or use eeg time - currind = eind(i); 
        currripend = pt_end(i);    
        currncells = rip_ncells(i)
        % Find all ripples within this window - Not really using this now. Using getripples_direct
        winst = currriptime-pret; % in ms
        winend = currriptime+postt; % in ms
        ripidx = find(pt>=(winst) & pt<=(winend));
        riptimec = pt(ripidx); % In ms
        riptimes_win = riptimec - currriptime; % Ripple time in current time axis
        
        % Pos and speed
        % -------------
        swinst = winst./1000; %sec
        swinend = winend./1000; %sec
        % Look up position index at this time in stateatrix. Same as index in pos file.
        posidx_st = lookup(swinst,statematrix.time); posidx_end = lookup(swinend,statematrix.time);
        currpostime = statematrix.time(posidx_st:posidx_end);
        currpostime = currpostime - currpostime(1); % Make time axis start at 0
        % Get speed from start to end of window
        posidx_st = lookup(swinst,postime); posidx_end = lookup(swinend,postime);
        currspeed = absvel(posidx_st:posidx_end);
        currposn = posn(posidx_st:posidx_end,:); % xy posn in window
        % Get posn during SWR, and n secs before and n secs after
        ripposidx_st = lookup(currriptime./1000,postime); ripposidx_end = lookup(currripend./1000,postime);
        ripposn = posn(ripposidx_st:ripposidx_end,:);
        surrposidx_st = lookup((currriptime./1000)-(pret/1000),postime); surrposidx_end = lookup((currripend./1000)+((postt/1000)),postime);
        preposn = posn(surrposidx_st:ripposidx_st,:);
        postposn = posn(ripposidx_end:surrposidx_end,:);
        
        
        
 

        % Main Figure
        % ----------
        
        figure(i); hold on;
        %redimscreen_halfvert;
        redimscreen;
    
        % First Speed on bottom
        % ----------------------
%         speedscale = max(currspeed)-min(currspeed);
%         plotscale = 5;
%         plotcurrspeed = currspeed.*(plotscale/speedscale); thrs = 5.*(plotscale/speedscale);
%         xpts = -pret:1:postt;
%         plot()        
%         %Divider line       
%         ypts = baseline*ones(size(xpts)); 
%         plot(xpts , ypts, 'k--','Linewidth',2);
        
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
        
        Hp_baseline = baseline; % Save for linearized position
        % Update baseline
        baseline = baseline + size(cellsi,1)*2;
        up = baseline;
        % Linearized position (from well 1 and well 3) on top of Hp spikes
        % Have to choose wells based on which trajectory the animal is on
%         plotscale = up-Hp_baseline;
%         lindist1 = statematrix.linearDistanceToWells(posidx_st:posidx_end,1);
%         lindist1 = lindist1-min(lindist1);
%         distscale = max(lindist1)-min(lindist1);
%         lindist1 = Hp_baseline + lindist1.*(plotscale/distscale);
%         plot(currpostime,lindist1,'r-','LineWidth',4);
        
        %Divider line
        xpts = -pret:1:postt;
        ypts = baseline*ones(size(xpts)); 
        plot(xpts , ypts, 'k--','Linewidth',2);
        
        % EEG On main tet
        % --------------
        n = maineegidx;
        eegtet = e_stim{n};
        curreeg = eegtet(i,:);
        eegscale = max(curreeg)-min(curreeg);
        downeeg = baseline+2; upeeg = downeeg+8;
        plotscale = 8;
        curreeg = downeeg + (plotscale/2) + curreeg.*(plotscale/eegscale);   
        plot(taxis,curreeg,'k-','LineWidth',2);
        
         % Update baseline
        baseline = baseline + 8;
        
         % EEG in ripple band
        % --------------------------------------
        n = maineegidx;
        riptet = ripamp_stim{n};
        ripenvtet = ripenv_stim{n};
        currrip = riptet(i,:);
        %currsize = roundn(ripsize_stim(i,n),-2);
        eegscale = max(currrip)-min(currrip);
        downeeg = baseline+2; upeeg = downeeg+8;
        plotscale = 8;
        currrip = downeeg + (plotscale/2) + currrip.*(plotscale/eegscale);   
        plot(taxis,currrip,'r-','LineWidth',2);
                    
        % Update baseline
        baseline = baseline + 8;       
        
          % EEG in theta band 
        % --------------------------------------
        n = maineegidx;
        thetatet = thetaamp_stim{n};
        thetaenvtet = thetaenv_stim{n};
        currtheta = thetatet(i,:);
        %currsize = roundn(ripsize_stim(i,n),-2);
        eegscale = max(currtheta)-min(currtheta);
        downeeg = baseline+2; upeeg = downeeg+8;
        plotscale = 8;
        currtheta = downeeg + (plotscale/2) + currtheta.*(plotscale/eegscale);   
        plot(thetaaxis,currtheta,'c-','LineWidth',1);
      
        % Update baseline
        baseline = baseline + 8;
        
        
%         % EEG on PFC Tet instead of ripple band
%         % --------------------------------------
%         n = peegidx;
%         eegtet = e_stim{n};
%         riptet = ripamp_stim{n};
%         ripenvtet = ripenv_stim{n};
%         curreeg = eegtet(i,:);
%         currrip = riptet(i,:);
%         %currsize = roundn(ripsize_stim(i,n),-2);
%         eegscale = max(curreeg)-min(curreeg);
%         downeeg = baseline+2; upeeg = downeeg+8;
%         plotscale = 8;
%         curreeg = downeeg + (plotscale/2) + curreeg.*(plotscale/eegscale);
%         plot(taxis,curreeg,'b-','LineWidth',1);
%         
%         % Update baseline
%         baseline = baseline + 8;
        
        
        
        
        
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

        
        
        
        title(sprintf('SWR aligned: ncells = %d', currncells),'FontSize',18,'Fontweight','normal');
        xlabel('Time (ms)','FontSize',18,'Fontweight','normal');
        set(gca, 'YTick',[]);
               
        saveg=0;
        if saveg==1
            figfile = [figdir,prefix,'Day2Ep5RasterEEGegNo',num2str(i)]
            print('-dpdf', figfile); print('-djpeg', figfile); saveas(gcf,figfile,'fig');
        end
        
        
          % Plot Posn/Speed in separate window
        % ----------------------------
        figure(100*i); hold on; Clgy = [0.5 0.5 0.5];
        plot(posn(:,1),posn(:,2),'Color',Clgy,'LineWidth',0.5);
        % plot(posn(:,1),posn(:,2),'MarkerFaceColor',Clgy,'MarkerEdgeColor',Clgy,'MarkerSize',1);
        plot(preposn(:,1),preposn(:,2),'k-','Linewidth',2); plot(postposn(:,1),postposn(:,2),'b-','LineWidth',2); % Pre in black
        plot(ripposn(:,1),ripposn(:,2),'r.','MarkerSize',24); % Post in blue
%        
        figure(200*i); hold on;
        currpostime = currpostime - (pret/1000);
        plot(currpostime, currspeed,'k-');
        
        keyboard; % Pause after each plot
        
        
    end % end if i
        
        
        
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

