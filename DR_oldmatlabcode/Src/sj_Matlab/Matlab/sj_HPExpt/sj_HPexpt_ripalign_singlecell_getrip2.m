function [basefr0] = sj_HPexpt_ripalign_singlecell_getrip2(prefix, day, epoch, riptets, Hptets, PFCtet, PFCcell, binsize, saveg1, dospeed)

% Version2 - add statistics
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



% --------------- Parameters ---------------
pret=550; postt=550; %% Times to plot - For 2 sec
trim = 50;

cellcountthresh = 3;
smwin=10; %Smoothing Window

if isempty(binsize)
    binsize = 5;  %% ms, for MU Fir Rate
end

rwin = [0 150];
bwin = [-500 -300];
push = 500; % either bwin(2) or postt-trim=500
%type ='exc'; 
type='inh';
%type='com'; %complex

savedir = '/data25/sjadhav/HPExpt/ProcessedData/';
savefile = [savedir 'HP_ripplemod_PFC'];
savedata = 0;
saveg1 = 0;

binsize_plot=binsize; %% If you put binsize_plot=1000, then units are Nspikes/binsize, not inst. firing rate in Hz
timeaxis = -pret:binsize:postt;
bins_center = find(abs(timeaxis)<=100);
bins_resp = find(timeaxis>=rwin(1) & timeaxis<=rwin(2));
bins_bck = find(timeaxis>=bwin(1) & timeaxis<=bwin(2));

bins_respHp = find(timeaxis>=0 & timeaxis<=100); % For quantifying Hippocampal response


% Speed parameters
lowsp_thrs = 2; % cm/sec
highsp_thrs = 7; % cm/sec

% Sleep
if ismember(epoch,[1 3 5 7]),
    lowsp_thrs = 0.5; %cm/sec
    highsp_thrs = 2; % cm/sec
end


plot_ex=0; % For plotting example single-trial MU rate plots
sd=2; %% SD for ripples - Not Really using this now. Taking all ripples in "ripples" file. Usually sd is 2.


doraster_mu=0;
plotline1 = 100; plotline2=[]; % Where to plot lines after stimulation
Screen = get(0,'Screensize');

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
% ----------------------
cntrip=0; bckresp=[]; rresp=[];
len_bckwin = bwin(2)-bwin(1); % eg 200ms

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
    
    % Get bin-values is resp window and back win - Raw or smoothened?
    try
        bckresp(:,cntrip) = histspks(bins_bck); % Stack rather than vectorize
        rresp(:,cntrip) = histspks(bins_resp);
    catch
        keyboard;
    end
    
    
    
end


% Align to similar number of triggers in bck window randomly
% ------------------------------------------------------------------
len_bckwin = bwin(2)-bwin(1); % eg 200ms
len_respwin = rwin(2)-rwin(1); % eg. 400ms 
% Random points have to at least len_respwin = 400ms away?
% CAN ALSO DO - NEED TO BE RWIN(2) AWAY - TIME AFTER 0
% Look in a 100ms window which is at least len_respwin away

%push = max(len_respwin,abs(bwin(2))); %eg if len_respwin=100, and bwin(2)=-300, then push by 300.
% OR you can just do postt-trim (usually 500 ms)
%push = 500;


cntbck=0;
for i=2:length(rip_starttime)-1
    cntbck=cntbck+1;
    %rdmidx = randperm(len_bckwin); % permute len_bckwin 1ms bins 
    %rdmidx = rdmidx(1);
    %currtrig = bwin(1)+rdmidx; % This is time in ms. No need for binsize here. Some random time in bck window
    rdmidx = randperm(100); % permute 100 1ms bins 
    rdmidx = rdmidx(1);
    currtrig = push+rdmidx; % This is time in ms. No need for binsize here. Some random time in bck window
    
    currrip = rip_starttime(i);
    currrip = currrip - currtrig; % Get rdm time behind actual ripple time in given window

    currspks =  spikeu_PFC(find( (spikeu_PFC>=(currrip-pret)) & (spikeu_PFC<=(currrip+postt)) ));
    currspks = currspks-(currrip);
    histspks = histc(currspks,[-pret:binsize:postt]);
    rdm_spks_cell_PFC{cntbck}=currspks;
    rdm_spkshist_cell_PFC(cntbck,:) = histspks;
    histspks = smoothvect(histspks, g1);
    rdm_spkshist_cell_PFC(cntbck,:) = histspks;

    % Get bin-values is resp window and back win - Raw or smoothened?
    rdm_bckresp(:,cntbck) = histspks(bins_bck);
    rdm_rresp(:,cntbck) = histspks(bins_resp);
end


% % TRIM ENDS - Due to Time Smoothing
rip_spkshist_cell_Hp = rip_spkshist_cell_Hp(:,trim/binsize:(pret+postt-trim)/binsize);
rip_spkshist_cell_PFC = rip_spkshist_cell_PFC(:,trim/binsize:(pret+postt-trim)/binsize);
% Also do the rdm hist
rdm_spkshist_cell_PFC = rdm_spkshist_cell_PFC(:,trim/binsize:(pret+postt-trim)/binsize);

% Update pret and postt due to end trims
pret=pret-trim; postt=postt-trim;

% Redo bins_resp
timeaxis = -pret:binsize:postt;
bins_center = find(abs(timeaxis)<=100);
bins_resp = find(timeaxis>=rwin(1) & timeaxis<=rwin(2));
bins_bck = find(timeaxis>=bwin(1) & timeaxis<=bwin(2));


% Try Raw rates instead of Z-scores. For Z-scores, see commented part in Version 1
% Sort by Fir rate in 0 to 100 ms for Hp, and apply same order to PFC/sort PFC by itself
% AND CONVERT TO FIRING RATES
% Keep a separate matrix for samae order for PFC to see trends over time -
% this is rip_spkshist_cell_PFC and rip_spkshist_cell_PFC
% % ----------------------------------------------------------------------
cellresp = sum(rip_spkshist_cell_Hp(:,bins_respHp),2);
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





%%-------------------------------------------------
% Stats
% -------------------------------------------------

sig_ttest = 0; sig_shuf = 0;

% Bck
avgbckresp_trial = mean(bckresp,1); % Mean in bck for each ripple
avgbckhist = mean(bckresp,2); % Avg bck histogram
mean_bckresp = mean(mean(bckresp)); % Single value
distr_bckresp = bckresp(:); %All values taken by bins in background

% Response
avgresp_trial = mean(rresp,1); % Mean for each ripple
avgresphist = mean(rresp,2); % Avg resp histogram
mean_rresp = mean(mean(rresp)); % Single value


% 1) Significance test
% --------------------

% First find locatin of peak or trough in response window
% Simply find max or mean in average histogram, or do a fit
% -------------------------------------------------

% Response type - Exc or inhibited. How to deal with complicated response - Shuffle method
if strcmp(type,'exc'),
    [peak, peakloc] = max(avgresphist);
    peakloc = peakloc(1); % incase there are multiple locations
    peakresp_trial = max(rresp,[],1); % Peak for each ripple. No need to use peak loc
end
if strcmp(type,'inh'),
    [peak, peakloc] = min(avgresphist);
    peakloc = peakloc(1); % incase there are multiple locations
    peakresp_trial = rresp(peakloc,:); % Use peak locn to get min value for each ripple
end

% Alternate method: Fit histogram
% window = (rwin(2)-rwin(1))/2;
% [~,fmu,confl,confu] = barfit_psth1(avgresphist(2:end)', window, binsize);
% % Need nbars. See in sj_Utility, PhDCodes/Chr_Analysis/Utility, or smkin 


% Modln and stats
modln = abs(100*(mean(peakresp_trial)-mean_bckresp)/mean_bckresp); % %tage change above/below baseline
[hpeak,ppeak] = ttest2(avgbckresp_trial,peakresp_trial); % Test if peak is significantly diff from backgnd in each trial  

if (ppeak < 0.05), 
    sig_ttest=1; 
end


% 2) Stats - Shuffle
% ------------------
% Rdm resp
avgrdmresp_trial = mean(rdm_rresp,1); % Mean for each ripple
avgrdmresphist = mean(rdm_rresp,2); % Avg resp histogram
mean_rdm_rresp = mean(mean(rdm_rresp)); % Single value
%figure; hold on; plot(avgresphist,'b'); plot(avgrdmresphist,'g')

% Get Distance between avg response and random response
D = dist(avgresphist,avgrdmresphist); % Global distance between vectors. Can also do ptwise

% Shuffle trials between real and rdm triggers, and generate D-distribution
combresp = [rresp, rdm_rresp];
ntr = size(combresp,2); % Ntr rresp + Ntr rdm_rresp. Both must be the same length
nshuffles = 1000;
for shufidx = 1:1000
    order = randperm(ntr);
    shuffle = combresp(:,order);
    
    shufresp = shuffle(:,1:ntr/2); shufavgresp = mean(shufresp,2);
    shufrdm = shuffle(:,(ntr/2)+1:ntr); shufavgrdmresp = mean(shufrdm,2);
    %figure; hold on; plot(shufavgresp,'b'); plot(shufavgrdmresp,'g')
    Dshuf(shufidx) = dist(shufavgresp, shufavgrdmresp);  
end
%histD = histc(Dshuf,min(Dshuf):0.01:max(Dshuf));
%figure; hold on; plot([min(Dshuf):0.01:max(Dshuf)],histD)
if D > prctile(Dshuf,95)
    sig_shuf = 1;
end

% Get the p-value of shuffle. The modulation index for shuffle will be the prctile value
pshuf = length(find(Dshuf>D))/nshuffles;
modln_shuf = 100 - (pshuf*100); %eg p=0.05 => prctile is 95%


% -----------
% SAVING DATA 
% ------------

if (savedata)
    load savefile;
    
    i = length(allripplemod)+1;  % idx to save
    
    % Index
    if strcmp(prefix,'HPa'), 
        anim=1;
    end
    if strcmp(prefix,'HPb'), 
        anim=2;
    end
    
    animdaytetcell = [anim day PFCtet PFCcell];
    allripplemod_idx(i,:) = animdaytetcell; 
    allripplemod(i).index = animdaytetcell;
    
    % Histograms
    allripplemod(i).ripplehist = rip_spkshist_cell_PFC;
    allripplemod(i).ripplehistsort = rip_spkshist_cellsort_PFC;
    allripplemod(i).ripplehist_Hp = rip_spkshist_cell_Hp;
    
    % Rasters - raw spiketimes within -pret to postt
    allripplemod(i).rippleraster = rip_spks_cell_PFC;
    allripplemod(i).ripplerastersort = rip_spks_cellsort_PFC;
    allripplemod(i).rippleraster_Hp = rip_spks_cell_Hp;
    
    % Stats
    allripplemod(i).modln = modln;
    allripplemod(i).sig_ttest = sig_ttest;
    allripplemod(i).rresp = rresp;
    allripplemod(i).bckresp = bckresp;
    
    % Shuffle
    allripplemod(i).modln_shuf = modln_shuf;
    allripplemod(i).sig_shuf = sig_shuf;
    allripplemod(i).pshuf = pshuf;
    allripplemod(i).rdm_rresp = rdm_resp;
    allripplemod(i).D = D;
    allripplemod(i).Dshuf = Dshuf;
    allripplemod(i).push = push;
    allripplemod(i).rdm_spks_cell_PFC = rdm_spks_cell_PFC;
    allripplemod(i).rdm_spkshist_cell_PFC = rdm_spkshist_cell_PFC;
    
    %Parameters
    allripplemod(i).rwin = rwin;
    allripplemod(i).bwin = bwin;
    allripplemod(i).pret = pret;
    allripplemod(i).postt = postt;
    allripplemod(i).trim = trim;
    allripplemod(i).cellcountthresh = cellcountthresh;
     
end





% 3) Stats - Zscore

% 4) Stats - prctile method

% sig=0;
% if mean_rresp >= mean_bckresp   % compare means?
%     type = 'exc';
% else
%     type = 'inh';
% end
% 
% if strcmp(type,'exc'),
%     avgresp = mean(rresp,2); % Histogram in response window
%     peak = max(avgresp);
%     if peak > prctile(distr_bckresp,95)
%         sig = 1;
%     end
% end
% 
% if strcmp(type,'inh'),
%     avgresp = mean(rresp,2); % Histogram in response window
%     trough = min(avgresp);
%     if trough < prctile(distr_bckresp,95)
%         sig = 1;
%     end
% end
    
   





%% ------------------------------------------------
% PLOT
% -------------------------------------------------

forppr = 0; 
% If yes, everything set to redimscreen_figforppr1
% If not, everything set to redimscreen_figforppt1

figdir = '/data25/sjadhav/HPExpt/Figures/RippleMod/';
summdir = figdir;
set(0,'defaultaxesfontweight','normal'); set(0,'defaultaxeslinewidth',2);

if forppr==1
    set(0,'defaultaxesfontsize',16); tfont = 18; xfont = 16; yfont = 16;
else
    set(0,'defaultaxesfontsize',24); tfont = 28; xfont = 20; yfont = 20;
end



% 1) Histogram-PFC Cell
% -------------

figure; hold on; redimscreen_figforppt1;
set(gcf, 'Position',[205 136 723 446]);

xaxis = -pret:binsize:postt;
plot(xaxis,mean(rip_spkshist_cellsort_PFC),'Linewidth',3);
plot(xaxis,mean(rip_spkshist_cellsort_PFC)+sem(rip_spkshist_cellsort_PFC),'b--','Linewidth',1);
plot(xaxis,mean(rip_spkshist_cellsort_PFC)-sem(rip_spkshist_cellsort_PFC),'b--','Linewidth',1);

set(gca,'XLim',[-pret postt]);
xlabel('Time(ms)','FontSize',xfont,'Fontweight','normal');
set(gca,'XTick',[-pret:200:postt],'XTickLabel',num2str([-pret:200:postt]'));

ylow = min(mean(rip_spkshist_cellsort_PFC)-sem(rip_spkshist_cellsort_PFC));
yhigh = max(mean(rip_spkshist_cellsort_PFC)+sem(rip_spkshist_cellsort_PFC));
set(gca,'YLim',[ylow-0.1 yhigh+0.1]);

ypts = ylow-0.1:0.1:yhigh+0.1;
xpts = 0*ones(size(ypts));
% Plot Line at 0 ms - Onset of stimulation
plot(xpts , ypts, 'k--','Linewidth',2);
% Plot lines at rwin and bckwi
xpts = rwin(1)*ones(size(ypts)); plot(xpts , ypts, 'k--','Linewidth',1);
xpts = rwin(2)*ones(size(ypts)); plot(xpts , ypts, 'k--','Linewidth',1);
xpts = bwin(1)*ones(size(ypts)); plot(xpts , ypts, 'g--','Linewidth',1);
xpts = bwin(2)*ones(size(ypts)); plot(xpts , ypts, 'g--','Linewidth',1);

if sig_ttest ==1, str = '*'; else, str = ''; end
if sig_shuf ==1, str_shuf = '*'; else, str_shuf = ''; end
title(sprintf('%s Day%d Tet%d Cell%d: Modln %g%s PrctileShuf %g%s', prefix, day, PFCtet, PFCcell, modln, str, modln_shuf, str_shuf),...
    'FontSize',tfont,'Fontweight','normal');

if saveg1==1,
    figfile = [figdir,'RippleAlignHist_',prefix,'_Day',num2str(day),'_Tet',num2str(PFCtet),num2str(PFCcell)]; 
    print('-dpdf', figfile); print('-djpeg', figfile); saveas(gcf,figfile,'fig');
end



% Raster
% ------
figure; hold on; redimscreen_figforppt1;
set(gcf, 'Position',[205 658 723 446]);
spkcount = [];
for c=1:length(rip_spks_cellsort_PFC)
    tmps = rip_spks_cellsort_PFC{c};
    plotraster(tmps,(length(rip_spks_cellsort_PFC)-c+1)*ones(size(tmps)),0.8,[],'LineWidth',3,'Color','b');
    % Get count of spikes in response window
    if ~isempty(tmps)
        subset_tmps = find(tmps>=rwin(1) & tmps<=rwin(2));
        spkcount = [spkcount; length(subset_tmps)];
    end
end
set(gca,'XLim',[-pret postt]);
set(gca,'XTick',[-pret:200:postt],'XTickLabel',num2str([-pret:200:postt]'));

set(gca,'YLim',[0 size(rip_spkshist_cellsort_PFC,1)]);
% Plot Line at 0 ms and rwin 
ypts = 0:1:size(rip_spkshist_cellsort_PFC,1);
xpts = 0*ones(size(ypts));
plot(xpts , ypts, 'k--','Linewidth',2);
% Plot lines at rwin 
xpts = rwin(1)*ones(size(ypts)); plot(xpts , ypts, 'k--','Linewidth',1);
xpts = rwin(2)*ones(size(ypts)); plot(xpts , ypts, 'k--','Linewidth',1);
xpts = bwin(1)*ones(size(ypts)); plot(xpts , ypts, 'g--','Linewidth',1);
xpts = bwin(2)*ones(size(ypts)); plot(xpts , ypts, 'g--','Linewidth',1);

title(sprintf('Ripple align: %s Day %d Tet %d Cell %d Nspk in rwin %d', prefix, day, PFCtet, PFCcell, sum(spkcount)),...
    'FontSize',tfont,'Fontweight','normal');

if saveg1==1,
    figfile = [figdir,'RippleAlignRaster_',prefix,'_Day',num2str(day),'_Tet',num2str(PFCtet),num2str(PFCcell)]; 
    print('-dpdf', figfile); print('-djpeg', figfile); saveas(gcf,figfile,'fig');
end

% Histogram with Random used for shuffling
% ----------------------------------------
figure; hold on; redimscreen_figforppt1;
set(gcf, 'Position',[960 658 723 446]);

xaxis = -pret:binsize:postt;
plot(xaxis,mean(rip_spkshist_cellsort_PFC),'Linewidth',3);
plot(xaxis,mean(rip_spkshist_cellsort_PFC)+sem(rip_spkshist_cellsort_PFC),'b--','Linewidth',1);
plot(xaxis,mean(rip_spkshist_cellsort_PFC)-sem(rip_spkshist_cellsort_PFC),'b--','Linewidth',1);

set(gca,'XLim',[-pret postt]);
xlabel('Time(ms)','FontSize',xfont,'Fontweight','normal');
set(gca,'XTick',[-pret:200:postt],'XTickLabel',num2str([-pret:200:postt]'));

ylow = min(mean(rip_spkshist_cellsort_PFC)-sem(rip_spkshist_cellsort_PFC));
yhigh = max(mean(rip_spkshist_cellsort_PFC)+sem(rip_spkshist_cellsort_PFC));
set(gca,'YLim',[ylow-0.1 yhigh+0.1]);

ypts = ylow-0.1:0.1:yhigh+0.1;
xpts = 0*ones(size(ypts));
% Plot Line at 0 ms - Onset of stimulation
plot(xpts , ypts, 'k--','Linewidth',2);
% Plot lines at rwin and bckwi
xpts = rwin(1)*ones(size(ypts)); plot(xpts , ypts, 'k--','Linewidth',1);
xpts = rwin(2)*ones(size(ypts)); plot(xpts , ypts, 'k--','Linewidth',1);
xpts = bwin(1)*ones(size(ypts)); plot(xpts , ypts, 'g--','Linewidth',1);
xpts = bwin(2)*ones(size(ypts)); plot(xpts , ypts, 'g--','Linewidth',1);

% Plot Random
rdm_spkshist_cell_PFC = rdm_spkshist_cell_PFC.*(1000/binsize);
plot(xaxis,mean(rdm_spkshist_cell_PFC),'g','Linewidth',1.5);
plot(xaxis,mean(rdm_spkshist_cell_PFC)+sem(rip_spkshist_cell_PFC),'g--','Linewidth',1);
plot(xaxis,mean(rdm_spkshist_cell_PFC)-sem(rip_spkshist_cell_PFC),'g--','Linewidth',1);

if sig_ttest ==1, str = '*'; else, str = ''; end
if sig_shuf ==1, str_shuf = '*'; else, str_shuf = ''; end
title(sprintf('%s Day%d Tet%d Cell%d: Modln %g%s PrctileShuf %g%s', prefix, day, PFCtet, PFCcell, modln, str, modln_shuf, str_shuf),...
    'FontSize',tfont,'Fontweight','normal');

if saveg1==1,
    figfile = [figdir,prefix,'_Day',num2str(day),'_Tet',num2str(tetu),'_PFCCellMatrix']; print('-dpdf', figfile); print('-djpeg', figfile); saveas(gcf,figfile,'fig');
end



% keyboard;



% %subplot(2,1,1); hold on;
% spkcount=[];
% for c=1:length(rip_spks_cell_PFC)
%     tmps = rip_spks_cell_PFC{c};
%     if ~isempty(tmps)
%         subset_tmps = find(tmps>=-200 & tmps<=200);
%         spkcount = [spkcount; length(subset_tmps)];
%     end
%     plotraster(tmps,c*ones(size(tmps)),0.8,[],'LineWidth',3,'Color','b');
% end
% set(gca,'XLim',[-pret postt]);
% set(gca,'YLim',[0 size(rip_spkshist_cellsort_PFC,1)]);
% % Plot Line at 0 ms and 100 ms - Onset of stimulation
% ypts = 0:1:size(rip_spkshist_cellsort_PFC,1);
% xpts = 0*ones(size(ypts));
% plot(xpts , ypts, 'k--','Linewidth',2);
% xpts = 100*ones(size(ypts));
% plot(xpts , ypts, 'k--','Linewidth',2);
% title(sprintf('Original order. NSpikes -200 to 200: %d ', sum(spkcount)));
















% ---------------
% From Version 1
% ---------------


% figure; hold on;
% %redimscreen_figforppt1;
% %set(gcf,'Position',[Screen(3)*0.205 Screen(4)*0.55 Screen(3)*0.2 Screen(4)*0.4])
% set(gcf,'Position',[984 42 900 1045]);
%
% % Imagesc
% %--------
% title(['PFC CellComb Firing Z aligned to ripples-' num2str(binsize) 'ms bins'],...
%     'FontSize',20,'Fontweight','normal');
% subplot(3,1,[1 2]); hold on;
% xaxis = -pret:binsize:postt; yaxis = 1:size(rip_spkshist_cellsort_PFC_matrix,1);
% imagesc(xaxis,yaxis,flipud(rip_spkshist_cellsort_PFC_matrix));
%
% %ylabel('Stim Amp / Laser Power','FontSize',20,'Fontweight','normal');
% %axis([0 ((pret+postt)/binsize)+2 0.5 nranges+0.5]);
% set(gca,'XLim',[-pret postt]);
% set(gca,'YLim',[0 size(rip_spkshist_cellsort_PFC_matrix,1)]);
% %xlabel('Time(ms)','FontSize',20,'Fontweight','normal');
%
% ypts = 0:1:size(rip_spkshist_cellsort_PFC_matrix,1);
% xpts = 0*ones(size(ypts));
%
% % Plot Line at 0 ms - Onset of stimulation
% plot(xpts , ypts, 'k--','Linewidth',3);
% % Plot lines at 100ms and 2000ms
% xpts = (plotline1)*ones(size(ypts));
% plot(xpts , ypts, 'k--','Linewidth',2);
% if ~isempty(plotline2)
%     xpts = (plotline2)*ones(size(ypts));
%     plot(xpts , ypts, 'k--','Linewidth',2);
% end
% % Ticks every 100ms
% %     set(gca,'xtick',[-pret:(100/binsize):postt],'xticklabel',{num2str([-pret:100:postt]')},...
% %         'FontSize',14,'Fontweight','normal');
%
% subplot(3,1,3); hold on;
% xaxis = -pret:binsize:postt;
% plot(xaxis,mean(rip_spkshist_cellsort_PFC),'Linewidth',3);
% plot(xaxis,mean(rip_spkshist_cellsort_PFC)+sem((rip_spkshist_cellsort_PFC)),'b--','Linewidth',1);
% plot(xaxis,mean(rip_spkshist_cellsort_PFC)-sem((rip_spkshist_cellsort_PFC)),'b--','Linewidth',1);
%
% set(gca,'XLim',[-pret postt]);
% xlabel('Time(ms)','FontSize',20,'Fontweight','normal');
% set(gca,'YLim',[min(mean(rip_spkshist_cellsort_PFC))-0.5 max(mean(rip_spkshist_cellsort_PFC))+0.5]);
% ypts = min(mean(rip_spkshist_cellsort_PFC))-0.5:0.1:max(mean(rip_spkshist_cellsort_PFC))+0.5;
% xpts = 0*ones(size(ypts));
% % Plot Line at 0 ms - Onset of stimulation
% plot(xpts , ypts, 'k--','Linewidth',2);
% % Plot lines at 100ms and 2000ms
% xpts = plotline1*ones(size(ypts));
% plot(xpts , ypts, 'k--','Linewidth',1);
%
% if saveg1==1,
%     figfile = [figdir,prefix,'_Day',num2str(day),'_Tet',num2str(tetu),'_PFCCellMatrix'];
%     print('-dpdf', figfile);
%     print('-djpeg', figfile);
%     saveas(gcf,figfile,'fig');
% end




% % Raster
% % ------
% figure; hold on;
% set(gcf, 'Position',[250 465 640 590]);
% %     title(['PFC CellComb Firing Z aligned to ripples-' num2str(binsize) 'ms bins'],...
% %         'FontSize',20,'Fontweight','normal');
% subplot(2,1,1); hold on;
% spkcount=[];
% for c=1:length(rip_spks_cell_PFC)
%     tmps = rip_spks_cell_PFC{c};
%     if ~isempty(tmps)
%         subset_tmps = find(tmps>=-200 & tmps<=200);
%         spkcount = [spkcount; length(subset_tmps)];
%     end
%     plotraster(tmps,c*ones(size(tmps)),0.8,[],'LineWidth',3,'Color','b');
% end
% set(gca,'XLim',[-pret postt]);
% set(gca,'YLim',[0 size(rip_spkshist_cellsort_PFC,1)]);
% ypts = 0:1:size(rip_spkshist_cellsort_PFC,1);
% xpts = 0*ones(size(ypts));
% % Plot Line at 0 ms and 100 ms - Onset of stimulation
% plot(xpts , ypts, 'k--','Linewidth',2);
% xpts = 100*ones(size(ypts));
% plot(xpts , ypts, 'k--','Linewidth',2);
% title(sprintf('Original order. NSpikes -200 to 200: %d ', sum(spkcount)));
%
% subplot(2,1,2); hold on;
% for c=1:length(rip_spks_cellsort_PFC)
%     tmps = rip_spks_cellsort_PFC{c};
%     plotraster(tmps,(length(rip_spks_cellsort_PFC)-c+1)*ones(size(tmps)),0.8,[],'LineWidth',3,'Color','b');
% end
% set(gca,'XLim',[-pret postt]);
% set(gca,'YLim',[0 size(rip_spkshist_cellsort_PFC,1)]);
% ypts = 0:1:size(rip_spkshist_cellsort_PFC,1);
% xpts = 0*ones(size(ypts));
% % Plot Line at 0 ms and 100 ms - Onset of stimulation
% plot(xpts , ypts, 'k--','Linewidth',2);
% xpts = 100*ones(size(ypts));
% plot(xpts , ypts, 'k--','Linewidth',2);
% title('Sorted');
% %axis off
% % Plot lines at 100ms and 2000ms
% %xpts = (plotline1)*ones(size(ypts));
% %plot(xpts , ypts, 'k--','Linewidth',2);






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



