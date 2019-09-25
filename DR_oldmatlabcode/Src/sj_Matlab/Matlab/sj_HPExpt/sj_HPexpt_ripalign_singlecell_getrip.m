function [basefr0] = sj_HPexpt_ripalign_singlecell_getrip(prefix, day, epoch, riptets, Hptets, PFCtet, PFCcell, binsize, saveg1, dospeed)
% Use getripples to look across multiple tetrodes

%  sj_HPexpt_ripalign_singlecell_getrip('HPa', 2, 4, [],[], 15, 1, 10, 0);

%  sj_HPexpt_ripalign_singlecell_getrip('HPa', 2, 4, [1 4 6 7 12 14],[1 4 6 7 12 14], 17, 2, 10, 0);
%  sj_HPexpt_ripalign_singlecell_getrip('HPa', 3, 2, [1 4 7 8 14],[1 4 7 8 14], 18,3, 10, 0);
%  sj_HPexpt_ripalign_singlecell_getrip('HPa', 4, 4, [1 4 7 8 14],[1 4 7 8 14], 15,3, 10, 0);

%  sj_HPexpt_ripalign_singlecell_getrip('HPb', 2, 4, [1 3 4 6 18],[1 3 4 6 18], 9, 1, 10, 0);
%  sj_HPexpt_ripalign_singlecell_getrip('HPb', 3, 4, [1 3 4 6 18 20],[1 3 4 6 18 20], 9, 1, 10, 0);
%  sj_HPexpt_ripalign_singlecell_getrip('HPb', 1, 4, [1,3,4,6,18],[1,3,4,6,18],9,1,10,0);

% Shantanu Nov 2012: From sj_HPexpt_burstaligntorip
% Instead of all cells on a tet, just do single PFCcell 
% sj_HPexpt_burstaligntorip_singlecell('HPb',1,4,4,[1,3,4,6,18],[9],[1],10,0);
% sj_HPexpt_burstaligntorip_singlecell('HPa',2,4,1,[1 4],[15],[1],10,0);

% Shantanu - Aug 2012. From RippleDisruption / sj_multiunitandspike_align
% Aligns bursts around ripples in increasing order of size for Hp and PFC. Similar to Wierzynski, 2009 - Fig. 8
% eg. sj_HPexpt_burstaligntorip('HPa',2,4,1,[1 4],[15 18],1,1,[],10,0);
% sj_HPexpt_burstaligntorip('HPb',1,5,4,[1 4 18],[9 10 14],1,1,[],10,0);

%[ Shantanu - May 2012. Also used to align to optical stimulation for animal JW7.
% From sj_multiunit_align, which is derived from sj_rippower_multiunit_calib.
% Change to getting times directly from epoch file. From given tetrode, plot mulitunit align, and spike align
% for either given cell or all cells on tetrode.
% Also, adding ability to use MU cluster rather than the multi file.] - Old

% eg. sj_multiunitandspike_align('JW7', 5, 1, 8, 1, [], 200, 0);


% Only aligning multiunit activity to DIO times.
% Enter times manually. Implement loading times file later

% From sj_rippower_stimampl3, which split into sj_rippower_multiunit_run and sj_rippower_multiunit_calib
% Only for calibration now, unlike sj_rippower_stimampl3  - needs manual input of time ranges



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
if nargin<4,
    riptets=[]; %
end
if nargin<5,
    Hptets=[]; %
end
if nargin<6,
    PFCtet=[]; %
end
if nargin<7,
    PFCcell=0; %
end
if nargin<8,
    binsize=10;
end
if nargin<9
    saveg1=0;
end
if nargin<10
    dospeed=[];
end

doraster_mu=0;
plotline1 = 100; plotline2=[]; % Where to plot lines after stimulation
Screen = get(0,'Screensize');

% --------------- Parameters ---------------
pret=550; postt=550; %% Times to plot - For 2 sec
trim = 50;
plot_ex=0; % For plotting example single-trial MU rate plots
sd=2; %% SD for ripples - Not Really using this now. Taking all ripples in "ripples" file. Usually sd is 2.

cellcountthresh = 3;

smwin=10; %Smoothing Window

twoplot_idx=1;
if isempty(binsize)
    binsize = 5;  %% ms, for MU Fir Rate
end
binsize_plot=binsize; %% If you put binsize_plot=1000, then units are Nspikes/binsize, not inst. firing rate in Hz
timeaxis = -pret:binsize:postt;
bins_center = find(abs(timeaxis)<=100);
bins_resp = find(timeaxis>=0 & timeaxis<=100);

% Speed parameters
lowsp_thrs = 2; % cm/sec
highsp_thrs = 7; % cm/sec

% Sleep
if ismember(epoch,[1 3 5 7]),
    lowsp_thrs = 0.5; %cm/sec
    highsp_thrs = 2; % cm/sec
end


% ------------------------------
% Figure and Font Sizes

forppr = 1;
% If yes, everything set to redimscreen_figforppr1
% If not, everything set to redimscreen_figforppt1

figdir = '/data25/sjadhav/';
%figdir = '/data25/sjadhav/RippleInterruption/Figures/01AugSep11_RippleDisFigs/DisruptnCalibrationAndEgs/CalibrationEgs/';
datadir = '/data25/sjadhav/';
summdir = figdir;

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
clr = {'b','r','g','c','m','y','k','r'};

% ---------------------------------------


%% -----------------------------------------
% SET DATA
% -------------------------------------------

switch prefix
    case 'HPa'
        rawdir = '/data25/sjadhav/HPExpt/HPa/';
        directoryname = '/data25/sjadhav/HPExpt/HPa_direct/';
%         rawdir = 'D:\1ShantanuWork\HPExpt\HPa';
%         directoryname = 'D:\1ShantanuWork\HPExpt\HPa_direct';
    case 'HPb'
        rawdir = '/data25/sjadhav/HPExpt/HPb/';
        directoryname = '/data25/sjadhav/HPExpt/HPb_direct/';
end


clr = {'b','g','c','m','y','k','r'};

%% --------------------------------------------------
%  Align MU and Cell-combine to Ripples
%  %%Align MU and Single Cell Firing Rate to Stimulations and Ripples
% -------------------------------------------------

% Load times file
% -----------------
currdir = pwd;
cd(rawdir);
dayfolders = dir;
daystr = sprintf('%02d', day);
for i = 3:length(dayfolders)
    if dayfolders(i).isdir
        if strcmp(daystr,dayfolders(i).name(1:2))
            disp(upper(dayfolders(i).name))
            cd(dayfolders(i).name);
            load times;
        end
    end
end
cd(currdir);
% Now Getting Range directly from times file
userange = ranges(epoch+1,:); % 1st row is allepochs. So you need +1
usename = names{epoch+1}(end-15:end);
amp=70; amp1=amp;
nranges=1;

% % Load dio file
% %------------------
% DIOfile = sprintf('%s/%sDIO%02d.mat', directoryname, prefix, day);
% load(DIOfile);
% stim = DIO{day}{epoch}{16};
% stim_starttime = stim.pulsetimes(:,1)./10; %ms
% stim_endtime = stim.pulsetimes(:,2)./10; %ms
% stim_length = stim.pulselength;
% stim_isi = stim.timesincelast(2:end)./10; %ms

% Load extracted ripple file, tetinfo and cellinfo -
% -------------------------------------------------
% SHOULD SWITCH THIS TO ALLTET. If alltet does not exist, then use given tet/tets.
ripfile = sprintf('%s/%sripples%02d.mat', directoryname, prefix, day);
load(ripfile);
spikefile = sprintf('%s/%sspikes%02d.mat', directoryname, prefix, day);
load(spikefile);
tetinfofile = sprintf('%s/%stetinfo.mat', directoryname, prefix);
load(tetinfofile);
cellinfofile = sprintf('%s/%scellinfo.mat', directoryname, prefix);
load(cellinfofile);

% Riptet
if isempty(riptets) % get from tetinfo if not given
    ntets = length(tetinfo{day}{epoch});
    for i=1:ntets
        if isfield(tetinfo{day}{epoch}{i},'descrip')
            if strcmp(tetinfo{day}{epoch}{i}.descrip,'riptet')
                riptets = [riptets; i];
            end
        end
    end     
end

% Ripple times across tetrode
% ---------------------------------
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
 

% Legacy way of doing single ripple
% ---------------------------------
% rip_starttime=[]; rip_sizes=[];
% for i=1:length(riptet)
%     currriptet=riptet(i);
%     rip_starttime = [rip_starttime; 1000* ripples{day}{epoch}{currriptet}.starttime];   % in msec
%     rip_sizes = [rip_sizes; ripples{day}{epoch}{currriptet}.maxthresh];   % in units of std dev
% end
[rip_starttime,sortidx] = sort(rip_starttime);
%rip_sizes = rip_sizes(sortidx);
% Find ripples within xxms of each other
%remidx = find(diff(rip_starttime)<=50)+1;
%rip_starttime(remidx)=[]; rip_sizes(remidx)=[];
triggers = rip_starttime;

% Find ripples separated by atleast a second
% --------------------------------------------
iri = diff(triggers);
keepidx = [1;find(iri>=1000)+1];
triggers = triggers(keepidx);


% Implement speed criterion
if ~isempty(dospeed)
    posfile = sprintf('%s/%spos%02d.mat', directoryname, prefix, day);
    load(posfile);
    absvel = abs(pos{day}{epoch}.data(:,5)); % Can also use field 9
    postime = pos{day}{epoch}.data(:,1); % in secs
    
    pidx = lookup(triggers,postime*1000);
    speed_atrip = absvel(pidx);
    lowsp_idx = find(speed_atrip <= lowsp_thrs);
    highsp_idx = find(speed_atrip >= highsp_thrs);
    
    if strcmp(dospeed,'low')==1
        triggers = triggers(lowsp_idx);
    end
    
    if strcmp(dospeed,'high')==1
        triggers = triggers(highsp_idx);
    end
end
rip_starttime = triggers;


% Hippocampal tetrodes for cells, if empty
% ---------------------------------
filterString = ' (   strcmp($area, ''CA1'') || strcmp($area, ''iCA1'')    )';
tetindices = evaluatefilter(tetinfo{day}{epoch}, filterString);
for tet=1:length(tetindices)
    if ~isempty(tetinfo{day}{epoch}{tet})
        if tetinfo{day}{epoch}{tet}.numcells>0
            Hptets = [Hptets; tet];
        end
    end 
end


% Load Spikes and MU file - MU file used in case MU cluster not defined in spikes file
% ------------------------------------------------------------------------------------

spikefile = sprintf('%s/%sspikes%02d.mat', directoryname, prefix, day);
if exist(spikefile)==2
    load(spikefile);
end
multifile = sprintf('%s/%smulti%02d.mat', directoryname, prefix, day);
load(multifile);

spikeu_Hp = []; spikeu_PFC=[]; %multiu_Hp=[]; multiu_PFC=[];
nHpcells=0; nPFCcells=0;

% Hp DATA
% ------
for ct=1:length(Hptets)
    tet=Hptets(ct);
    % Spikes
    % ------
    
    % Get Spike timestamps (in secs - mult by 1000 for ms)
    for i=1:length(spikes{day}{epoch}{tet})
        if ~isempty(spikes{day}{epoch}{tet}{i})
            %tet, i
            %if tet==4 & i==7, keyboard; end
            if ~isempty(spikes{day}{epoch}{tet}{i}.data)
                cmd=sprintf('spike%d = spikes{day}{epoch}{tet}{%d}.data(:,1)*1000;',i,i); eval(cmd);
                cmd=sprintf('spikeu_Hp = [spikeu_Hp;spike%d];',i); eval(cmd);
                nHpcells=nHpcells+1;
            end
        end
    end
end

% PFC DATA
% ------
spikeu_PFC = spikes{day}{epoch}{PFCtet}{PFCcell}.data(:,1)*1000;




% ALINGING TO RIPPLES
% ----------------------

% Smooth over time axis binsize*3 ms (30ms for 10ms; 15ms for 5 ms)
% nstd = round(binsize*3/binsize); g1 = gaussian(nstd, 5*nstd+1);
% gaussian(sigma,numpoints);
% nstd=1: gaussian of length 4. nstd = 2: gaussian of length 7, nstd=3: gaussian of length 10.
nstd = round(binsize*2/binsize); g1 = gaussian(nstd, 3*nstd+1);

% Align spikes to ripples
cntrip=0;
for i=2:length(rip_starttime)-1
    i;
    cntrip=cntrip+1;
    currrip = rip_starttime(i);
    %ripsize_cell(cntrip) = rip_sizes(i); ripsize_multi(cntrip) = rip_sizes(i);
    
    % Hp
    currspks =  spikeu_Hp(find( (spikeu_Hp>=(currrip-pret)) & (spikeu_Hp<=(currrip+postt)) ));
    currspks = currspks-(currrip);
    histspks = histc(currspks,[-pret:binsize:postt]);
    rip_spks_cell_Hp{cntrip}=currspks;
    rip_spkshist_cell_Hp(cntrip,:) = histspks;
    histspks = smoothvect(histspks, g1);
    rip_spkshist_cell_Hp(cntrip,:) = histspks;
    % PFC
    currspks =  spikeu_PFC(find( (spikeu_PFC>=(currrip-pret)) & (spikeu_PFC<=(currrip+postt)) ));
    currspks = currspks-(currrip);
    histspks = histc(currspks,[-pret:binsize:postt]);
    rip_spks_cell_PFC{cntrip}=currspks;
    rip_spkshist_cell_PFC(cntrip,:) = histspks;
    histspks = smoothvect(histspks, g1);
    rip_spkshist_cell_PFC(cntrip,:) = histspks;
end

% % TRIM ENDS - Due to Time Smoothing
rip_spkshist_cell_Hp = rip_spkshist_cell_Hp(:,trim/binsize:(pret+postt-trim)/binsize);
rip_spkshist_cell_PFC = rip_spkshist_cell_PFC(:,trim/binsize:(pret+postt-trim)/binsize);
% Update pret and postt due to end trims
pret=pret-trim; postt=postt-trim;

% Redo bins_resp
timeaxis = -pret:binsize:postt;
bins_center = find(abs(timeaxis)<=100);
bins_resp = find(timeaxis>=0 & timeaxis<=100);

% Try Raw rates instead of Z-scores
% Sort by Fir rate in -100 to 100 ms for Hp, and apply same order to PFC/sort PFC by itself
% AND CONVERT TO FIRING RATES
% Keep a separate matrix for samae order for PFC to see trends over time - 
% this is rip_spkshist_cell_PFC and rip_spkshist_cell_PFC
% % ----------------------------------------------------------------------
cellresp = sum(rip_spkshist_cell_Hp(:,bins_resp),2);
[cellresp_sort,cellsortidx] = sort(cellresp);
rip_spkshist_cellsort_Hp = rip_spkshist_cell_Hp(cellsortidx,:).*(1000/binsize);
%cellresp = sum(rip_spkshist_cellZ_PFC(:,bins_resp),2);
%[cellresp_sort,cellsortidx] = sort(cellresp);
rip_spkshist_cellsort_PFC = rip_spkshist_cell_PFC(cellsortidx,:).*(1000/binsize);

% Original matrices - convert to firing rate
rip_spkshist_cell_Hp = rip_spkshist_cell_Hp.*(1000/binsize);
rip_spkshist_cell_PFC = rip_spkshist_cell_PFC.*(1000/binsize);


% Re-order raster as well
% -----------------------
rip_spks_cellsort_PFC=[];
for i=1:length(rip_spks_cell_PFC)
    rip_spks_cellsort_PFC{i}=rip_spks_cell_PFC{cellsortidx(i)};
end


% -------------------------------------------------
% Smooth Matrix in Y-direction (along ripples). Only for Matrix
% Smooth along time already done
% -------------------------------------------------

% Sorted by Hp popln response
rip_spkshist_cellsort_Hp_matrix = rip_spkshist_cellsort_Hp;
rip_spkshist_cellsort_PFC_matrix = rip_spkshist_cellsort_PFC;

% Sorted by time
rip_spkshist_cell_Hp_matrix = rip_spkshist_cell_Hp;
rip_spkshist_cell_PFC_matrix = rip_spkshist_cell_PFC;

for i=1:size(rip_spkshist_cellsort_Hp_matrix,1)
    winst = i-smwin/2; winend = i+smwin/2;
    if winst<1, winst=1; end
    if winend>size(rip_spkshist_cellsort_Hp,1), winend = size(rip_spkshist_cellsort_Hp,1); end
    rip_spkshist_cellsort_Hp_matrix(i,:) = mean(rip_spkshist_cellsort_Hp_matrix(winst:winend,:));
    rip_spkshist_cellsort_PFC_matrix(i,:) = mean(rip_spkshist_cellsort_PFC_matrix(winst:winend,:));
    
    % Also smooth original matrix
    rip_spkshist_cell_Hp_matrix(i,:) = mean(rip_spkshist_cell_Hp_matrix(winst:winend,:));
    rip_spkshist_cell_PFC_matrix(i,:) = mean(rip_spkshist_cell_PFC_matrix(winst:winend,:));
  
end



% % Convert to Z-scores
% % -------------------
% if docellcomb==1
%     % 1) For cellspks
%     cellmean_Hp = mean(rip_spkshist_cell_Hp,2); cellstd_Hp = std(rip_spkshist_cell_Hp,[],2);
%     cellmean_PFC = mean(rip_spkshist_cell_PFC,2); cellstd_PFC = std(rip_spkshist_cell_PFC,[],2);
%     x=find(cellstd_Hp==0 | cellstd_PFC==0); % Remove rows which are all zero / 
%     % OR can set all these rows to zero in Z-score: Make means for these rows 0, and make std any n (say 1). Z-score wil then be 0 
%     rip_spkshist_cell_Hp(x,:)=[]; cellmean_Hp(x)=[]; cellstd_Hp(x)=[]; %rip_spks_cell{x}=[];
%     rip_spkshist_cell_PFC(x,:)=[]; cellmean_PFC(x)=[]; cellstd_PFC(x)=[]; ripsize_cell(x) = [];
%     % 1a) Hp
%     cellmean_mat = repmat(cellmean_Hp,1,size(rip_spkshist_cell_Hp,2)); cellstd_mat = repmat(cellstd_Hp,1,size(rip_spkshist_cell_Hp,2));
%     rip_spkshist_cellZ_Hp = (rip_spkshist_cell_Hp-cellmean_mat)./cellstd_mat;
%     % 1b) PFC
%     cellmean_mat = repmat(cellmean_PFC,1,size(rip_spkshist_cell_PFC,2)); cellstd_mat = repmat(cellstd_PFC,1,size(rip_spkshist_cell_PFC,2));
%     rip_spkshist_cellZ_PFC = (rip_spkshist_cell_PFC-cellmean_mat)./cellstd_mat;
%     
%     % A) SORT BY RIPPLE SIZE. Do either A or B
%     % -----------------------------------------
%     % [cellresp_sort,cellsortidx] = sort(ripsize_cell);
%     % rip_spkshist_cellZsort_Hp = rip_spkshist_cellZ_Hp(cellsortidx,:);
%     % rip_spkshist_cellZsort_PFC = rip_spkshist_cellZ_PFC(cellsortidx,:);
%     
%     % % B )Sort by Fir rate in -100 to 100 ms for Hp, and apply same order to PFC/sort PFC by itself
%     % % ----------------------------------------------------------------------
%     cellresp = sum(rip_spkshist_cellZ_Hp(:,bins_resp),2);
%     [cellresp_sort,cellsortidx] = sort(cellresp);
%     rip_spkshist_cellZsort_Hp = rip_spkshist_cellZ_Hp(cellsortidx,:);
%     %cellresp = sum(rip_spkshist_cellZ_PFC(:,bins_resp),2);
%     %[cellresp_sort,cellsortidx] = sort(cellresp);
%     rip_spkshist_cellZsort_PFC = rip_spkshist_cellZ_PFC(cellsortidx,:);
%     
%     % -------------------------------------------------
%     % Smooth Matrix
%     % -------------------------------------------------
%     
%     for i=1:size(rip_spkshist_cellZsort_Hp,1)
%         winst = i-smwin/2; winend = i+smwin/2;
%         if winst<1, winst=1; end
%         if winend>size(rip_spkshist_cellZsort_Hp,1), winend = size(rip_spkshist_cellZsort_Hp,1); end
%         rip_spkshist_cellZsort_Hp(i,:) = mean(rip_spkshist_cellZsort_Hp(winst:winend,:));
%         rip_spkshist_cellZsort_PFC(i,:) = mean(rip_spkshist_cellZsort_PFC(winst:winend,:));
%     end
%     % TRIM ENDS
%     rip_spkshist_cellZsort_Hp = rip_spkshist_cellZsort_Hp(:,50/binsize:(pret+postt-50)/binsize);
%     rip_spkshist_cellZsort_PFC = rip_spkshist_cellZsort_PFC(:,50/binsize:(pret+postt-50)/binsize);
%     
% end




%% ------------------------------------------------
% PLOT
% -------------------------------------------------


    
%     %% Imagesc-Hp Cell
%     % -------------
%     
%     figure; hold on;
%     %redimscreen_figforppt1;
%     %set(gcf,'Position',[0 Screen(4)*0.55 Screen(3)*0.2 Screen(4)*0.4])
%     set(gcf,'Position',[16 44 900 1045]);
%     subplot(3,1,[1 2]); hold on;
%     xaxis = -pret:binsize:postt; yaxis = 1:size(rip_spkshist_cellsort_Hp_matrix,1);
%     imagesc(xaxis,yaxis,flipud(rip_spkshist_cellsort_Hp_matrix));
%     title(['Hp CellComb Firing aligned to ripples-' num2str(binsize) 'ms bins'],...
%         'FontSize',20,'Fontweight','normal');
%     %ylabel('Stim Amp / Laser Power','FontSize',20,'Fontweight','normal');
%     %axis([0 ((pret+postt)/binsize)+2 0.5 nranges+0.5]);
%     set(gca,'XLim',[-pret postt]);
%     set(gca,'YLim',[0 size(rip_spkshist_cellsort_Hp_matrix,1)]);
%     %xlabel('Time(ms)','FontSize',20,'Fontweight','normal');
%     
%     ypts = 0:1:size(rip_spkshist_cellsort_Hp_matrix,1);
%     xpts = 0*ones(size(ypts));
%     
%     % Plot Line at 0 ms - Onset of stimulation
%     plot(xpts , ypts, 'k--','Linewidth',3);
%     % Plot lines at 100ms and 2000ms
%     xpts = (plotline1)*ones(size(ypts));
%     plot(xpts , ypts, 'k--','Linewidth',2);
%     if ~isempty(plotline2)
%         xpts = (plotline2)*ones(size(ypts));
%         plot(xpts , ypts, 'k--','Linewidth',2);
%     end
%     % Ticks every 100ms
% %     set(gca,'xtick',[-pret/binsize:(100/binsize):postt/binsize],'xticklabel',{num2str([-pret:100:postt]')},...
% %         'FontSize',14,'Fontweight','normal');
%     
%     subplot(3,1,3); hold on;
%     xaxis = -pret:binsize:postt;
%     plot(xaxis,mean(rip_spkshist_cellsort_Hp),'Linewidth',3);
%     set(gca,'XLim',[-pret postt]);
%     xlabel('Time(ms)','FontSize',20,'Fontweight','normal');
%     set(gca,'YLim',[min(mean(rip_spkshist_cellsort_Hp))-2 max(mean(rip_spkshist_cellsort_Hp))+2]);
%     ypts = min(mean(rip_spkshist_cellsort_Hp))-2:0.1:max(mean(rip_spkshist_cellsort_Hp))+2;
%     xpts = 0*ones(size(ypts));
%     % Plot Line at 0 ms - Onset of stimulation
%     plot(xpts , ypts, 'k--','Linewidth',2);
%     % Plot lines at 100ms and 2000ms
%     xpts = plotline1*ones(size(ypts));
%     plot(xpts , ypts, 'k--','Linewidth',1);
%     
%     
%     if saveg1==1,
%         figfile = [figdir,prefix,'_Day',num2str(day),'_Tet',num2str(tetu),'_HpCellMatrix'];
%         print('-dpdf', figfile);
%         print('-djpeg', figfile);
%         saveas(gcf,figfile,'fig');
%     end
    
    %% Imagesc-PFC Cell
    % -------------
    
    figure; hold on;
    %redimscreen_figforppt1;
    %set(gcf,'Position',[Screen(3)*0.205 Screen(4)*0.55 Screen(3)*0.2 Screen(4)*0.4])
    set(gcf,'Position',[984 42 900 1045]);
       
    % Imagesc
    %--------
    title(['PFC CellComb Firing Z aligned to ripples-' num2str(binsize) 'ms bins'],...
        'FontSize',20,'Fontweight','normal');
    subplot(3,1,[1 2]); hold on;
    xaxis = -pret:binsize:postt; yaxis = 1:size(rip_spkshist_cellsort_PFC_matrix,1);
    imagesc(xaxis,yaxis,flipud(rip_spkshist_cellsort_PFC_matrix));
    
    %ylabel('Stim Amp / Laser Power','FontSize',20,'Fontweight','normal');
    %axis([0 ((pret+postt)/binsize)+2 0.5 nranges+0.5]);
    set(gca,'XLim',[-pret postt]);
    set(gca,'YLim',[0 size(rip_spkshist_cellsort_PFC_matrix,1)]);
    %xlabel('Time(ms)','FontSize',20,'Fontweight','normal');
    
    ypts = 0:1:size(rip_spkshist_cellsort_PFC_matrix,1);
    xpts = 0*ones(size(ypts));
    
    % Plot Line at 0 ms - Onset of stimulation
    plot(xpts , ypts, 'k--','Linewidth',3);
    % Plot lines at 100ms and 2000ms
    xpts = (plotline1)*ones(size(ypts));
    plot(xpts , ypts, 'k--','Linewidth',2);
    if ~isempty(plotline2)
        xpts = (plotline2)*ones(size(ypts));
        plot(xpts , ypts, 'k--','Linewidth',2);
    end
    % Ticks every 100ms
%     set(gca,'xtick',[-pret:(100/binsize):postt],'xticklabel',{num2str([-pret:100:postt]')},...
%         'FontSize',14,'Fontweight','normal');
    
    subplot(3,1,3); hold on;
    xaxis = -pret:binsize:postt;
    plot(xaxis,mean(rip_spkshist_cellsort_PFC),'Linewidth',3);
    set(gca,'XLim',[-pret postt]);
    xlabel('Time(ms)','FontSize',20,'Fontweight','normal');
    set(gca,'YLim',[min(mean(rip_spkshist_cellsort_PFC))-0.5 max(mean(rip_spkshist_cellsort_PFC))+0.5]);
    ypts = min(mean(rip_spkshist_cellsort_PFC))-0.5:0.1:max(mean(rip_spkshist_cellsort_PFC))+0.5;
    xpts = 0*ones(size(ypts));
    % Plot Line at 0 ms - Onset of stimulation
    plot(xpts , ypts, 'k--','Linewidth',2);
    % Plot lines at 100ms and 2000ms
    xpts = plotline1*ones(size(ypts));
    plot(xpts , ypts, 'k--','Linewidth',1);
    
    if saveg1==1,
        figfile = [figdir,prefix,'_Day',num2str(day),'_Tet',num2str(tetu),'_PFCCellMatrix'];
        print('-dpdf', figfile);
        print('-djpeg', figfile);
        saveas(gcf,figfile,'fig');
    end


   
    % Raster
    % ------
     figure; hold on;
     set(gcf, 'Position',[250 465 640 590]);
%     title(['PFC CellComb Firing Z aligned to ripples-' num2str(binsize) 'ms bins'],...
%         'FontSize',20,'Fontweight','normal');
    subplot(2,1,1); hold on;
    spkcount=[];
    for c=1:length(rip_spks_cell_PFC)
        tmps = rip_spks_cell_PFC{c};
        if ~isempty(tmps)
            subset_tmps = find(tmps>=-200 & tmps<=200);
            spkcount = [spkcount; length(subset_tmps)];
        end
        plotraster(tmps,c*ones(size(tmps)),0.8,[],'LineWidth',3,'Color','b');
    end
    set(gca,'XLim',[-pret postt]);
    set(gca,'YLim',[0 size(rip_spkshist_cellsort_PFC,1)]);
    ypts = 0:1:size(rip_spkshist_cellsort_PFC,1);
    xpts = 0*ones(size(ypts));
    % Plot Line at 0 ms and 100 ms - Onset of stimulation
    plot(xpts , ypts, 'k--','Linewidth',2);
    xpts = 100*ones(size(ypts));
    plot(xpts , ypts, 'k--','Linewidth',2);
    title(sprintf('Original order. NSpikes -200 to 200: %d ', sum(spkcount)));
    
    subplot(2,1,2); hold on;
    for c=1:length(rip_spks_cellsort_PFC)
        tmps = rip_spks_cellsort_PFC{c};
        plotraster(tmps,(length(rip_spks_cellsort_PFC)-c+1)*ones(size(tmps)),0.8,[],'LineWidth',3,'Color','b');
    end
    set(gca,'XLim',[-pret postt]);
    set(gca,'YLim',[0 size(rip_spkshist_cellsort_PFC,1)]);
    ypts = 0:1:size(rip_spkshist_cellsort_PFC,1);
    xpts = 0*ones(size(ypts));
    % Plot Line at 0 ms and 100 ms - Onset of stimulation
    plot(xpts , ypts, 'k--','Linewidth',2);
    xpts = 100*ones(size(ypts));
    plot(xpts , ypts, 'k--','Linewidth',2);
    title('Sorted');
    %axis off
    % Plot lines at 100ms and 2000ms
    %xpts = (plotline1)*ones(size(ypts));
    %plot(xpts , ypts, 'k--','Linewidth',2);


    
% % Original matrices: Not Sorted by Popln response - order is time during epoch
%  %% Imagesc-Hp Cell
%     % -------------
%     
%     figure; hold on;
%     %redimscreen_figforppt1;
%     %set(gcf,'Position',[0 Screen(4)*0.55 Screen(3)*0.2 Screen(4)*0.4])
%     set(gcf,'Position',[16 44 900 1045]);
%     subplot(3,1,[1 2]); hold on;
%     xaxis = -pret:binsize:postt; yaxis = 1:size(rip_spkshist_cell_Hp_matrix,1);
%     imagesc(xaxis,yaxis,flipud(rip_spkshist_cell_Hp_matrix));
%     title(['Hp CellComb Firing aligned to ripples-' num2str(binsize) 'ms bins'],...
%         'FontSize',20,'Fontweight','normal');
%     %ylabel('Stim Amp / Laser Power','FontSize',20,'Fontweight','normal');
%     %axis([0 ((pret+postt)/binsize)+2 0.5 nranges+0.5]);
%     set(gca,'XLim',[-pret postt]);
%     set(gca,'YLim',[0 size(rip_spkshist_cell_Hp_matrix,1)]);
%     %xlabel('Time(ms)','FontSize',20,'Fontweight','normal');
%     
%     ypts = 0:1:size(rip_spkshist_cell_Hp_matrix,1);
%     xpts = 0*ones(size(ypts));
%     
%     % Plot Line at 0 ms - Onset of stimulation
%     plot(xpts , ypts, 'k--','Linewidth',3);
%     % Plot lines at 100ms and 2000ms
%     xpts = (plotline1)*ones(size(ypts));
%     plot(xpts , ypts, 'k--','Linewidth',2);
%     if ~isempty(plotline2)
%         xpts = (plotline2)*ones(size(ypts));
%         plot(xpts , ypts, 'k--','Linewidth',2);
%     end
%     % Ticks every 100ms
% %     set(gca,'xtick',[-pret/binsize:(100/binsize):postt/binsize],'xticklabel',{num2str([-pret:100:postt]')},...
% %         'FontSize',14,'Fontweight','normal');
%     
%     subplot(3,1,3); hold on;
%     xaxis = -pret:binsize:postt;
%     plot(xaxis,mean(rip_spkshist_cell_Hp),'Linewidth',3);
%     set(gca,'XLim',[-pret postt]);
%     xlabel('Time(ms)','FontSize',20,'Fontweight','normal');
%     set(gca,'YLim',[min(mean(rip_spkshist_cell_Hp))-2 max(mean(rip_spkshist_cell_Hp))+2]);
%     ypts = min(mean(rip_spkshist_cell_Hp))-2:0.1:max(mean(rip_spkshist_cell_Hp))+2;
%     xpts = 0*ones(size(ypts));
%     % Plot Line at 0 ms - Onset of stimulation
%     plot(xpts , ypts, 'k--','Linewidth',2);
%     % Plot lines at 100ms and 2000ms
%     xpts = plotline1*ones(size(ypts));
%     plot(xpts , ypts, 'k--','Linewidth',1);
%     
%     
%     if saveg1==1,
%         figfile = [figdir,prefix,'_Day',num2str(day),'_Tet',num2str(tetu),'_HpCellMatrix'];
%         print('-dpdf', figfile);
%         print('-djpeg', figfile);
%         saveas(gcf,figfile,'fig');
%     end
%     
%     %% Imagesc-PFC Cell
%     % -------------
%     
%     figure; hold on;
%     %redimscreen_figforppt1;
%     %set(gcf,'Position',[Screen(3)*0.205 Screen(4)*0.55 Screen(3)*0.2 Screen(4)*0.4])
%     set(gcf,'Position',[984 42 900 1045]);
%        
%     % Imagesc
%     %--------
%     title(['PFC CellComb Firing Z aligned to ripples-' num2str(binsize) 'ms bins'],...
%         'FontSize',20,'Fontweight','normal');
%     subplot(3,1,[1 2]); hold on;
%     xaxis = -pret:binsize:postt; yaxis = 1:size(rip_spkshist_cell_PFC_matrix,1);
%     imagesc(xaxis,yaxis,flipud(rip_spkshist_cell_PFC_matrix));
%     
%     %ylabel('Stim Amp / Laser Power','FontSize',20,'Fontweight','normal');
%     %axis([0 ((pret+postt)/binsize)+2 0.5 nranges+0.5]);
%     set(gca,'XLim',[-pret postt]);
%     set(gca,'YLim',[0 size(rip_spkshist_cell_PFC_matrix,1)]);
%     %xlabel('Time(ms)','FontSize',20,'Fontweight','normal');
%     
%     ypts = 0:1:size(rip_spkshist_cell_PFC_matrix,1);
%     xpts = 0*ones(size(ypts));
%     
%     % Plot Line at 0 ms - Onset of stimulation
%     plot(xpts , ypts, 'k--','Linewidth',3);
%     % Plot lines at 100ms and 2000ms
%     xpts = (plotline1)*ones(size(ypts));
%     plot(xpts , ypts, 'k--','Linewidth',2);
%     if ~isempty(plotline2)
%         xpts = (plotline2)*ones(size(ypts));
%         plot(xpts , ypts, 'k--','Linewidth',2);
%     end
%     % Ticks every 100ms
% %     set(gca,'xtick',[-pret:(100/binsize):postt],'xticklabel',{num2str([-pret:100:postt]')},...
% %         'FontSize',14,'Fontweight','normal');
%     
%     subplot(3,1,3); hold on;
%     xaxis = -pret:binsize:postt;
%     plot(xaxis,mean(rip_spkshist_cell_PFC),'Linewidth',3);
%     set(gca,'XLim',[-pret postt]);
%     xlabel('Time(ms)','FontSize',20,'Fontweight','normal');
%     set(gca,'YLim',[min(mean(rip_spkshist_cell_PFC))-0.5 max(mean(rip_spkshist_cell_PFC))+0.5]);
%     ypts = min(mean(rip_spkshist_cell_PFC))-0.5:0.1:max(mean(rip_spkshist_cell_PFC))+0.5;
%     xpts = 0*ones(size(ypts));
%     % Plot Line at 0 ms - Onset of stimulation
%     plot(xpts , ypts, 'k--','Linewidth',2);
%     % Plot lines at 100ms and 2000ms
%     xpts = plotline1*ones(size(ypts));
%     plot(xpts , ypts, 'k--','Linewidth',1);
%     
%     if saveg1==1,
%         figfile = [figdir,prefix,'_Day',num2str(day),'_Tet',num2str(tetu),'_PFCCellMatrix'];
%         print('-dpdf', figfile);
%         print('-djpeg', figfile);
%         saveas(gcf,figfile,'fig');
%     end
    
    
    
    
i=1;
%set(gcf,'Position',[Screen(3)*0.205 Screen(4)*0.03 Screen(3)*0.2 Screen(4)*0.4])

%figfile = ['Figures/JW7_d7t1c1_inh_hist100_1sec'];saveas(gcf,figfile,'fig');print('-djpeg', figfile);print('-dpdf', figfile);

































%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



