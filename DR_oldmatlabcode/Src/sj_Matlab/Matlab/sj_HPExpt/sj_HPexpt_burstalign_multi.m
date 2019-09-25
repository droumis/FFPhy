function [basefr0] = sj_HPexpt_burstalign(prefix, day, epoch, Hptets, PFCtets, docellcomb, domulti, mucluster, binsize, saveg1)

% Shantanu - Aug 2012. From busrtaligntorip. Here - get MU bursts and align to them, rather than aligning to ripples
% Shantanu - Aug 2012. From RippleDisruption / sj_multiunitandspike_align
% Aligns bursts around ripples in increasing order of size for Hp and PFC. Similar to Wierzynski, 2009 - Fig. 8
% eg. sj_HPexpt_burstalign('HPa',2,5,[12 14],[15 18],1,0,[],5,0);
% sj_HPexpt_burstalign('HPb',1,4,[4 16 18],[9 14],1,0,[],5,0);


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
    Hptets=[]; %
end
if nargin<5,
    PFCtets=[]; %
end
if nargin<6,
    docellcomb=0; %
end
if nargin<7,
    domulti=1; %
end
if nargin<8,
    mucluster=[]; % If no MU cluster is specified, then will look for Mu tag in spikes. If does not exist, will use MU file
end
if nargin<9,
    binsize=5;
end
if nargin<10
    saveg1=0;
end

doraster_mu=0;
plotline1 = 100; plotline2=[]; % Where to plot lines after stimulation
Screen = get(0,'Screensize');

% --------------- Parameters ---------------
pret=550; postt=550; %% Times to plot - For 2 sec
plot_ex=0; % For plotting example single-trial MU rate plots
dorip = 1; % SET TO 1 IF YOU ALSO WANT TO ALIGN TO RIPPLES
sd=2; %% SD for ripples - Not Really using this now. Taking all ripples in "ripples" file. Usually sd is 2.

smwin=100; %Smoothing Window

twoplot_idx=1;
if isempty(binsize)
    binsize = 5;  %% ms, for MU Fir Rate
end
binsize_plot=binsize; %% If you put binsize_plot=1000, then units are Nspikes/binsize, not inst. firing rate in Hz
timeaxis = -pret:binsize:postt;
bins_center = find(abs(timeaxis)<=100);
bins_resp = find(timeaxis>=0 & timeaxis<=100);

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

% Load extracted ripple file -
% --------------------------
% SHOULD SWITCH THIS TO ALLTET. If alltet does not exist, then use given tet.
%ripfile = sprintf('%s/%sripples%02d.mat', directoryname, prefix, day);
%load(ripfile);
%rip_starttime = 1000* ripples{day}{epoch}{Hptets(1)}.starttime;   % in msec
%rip_sizes = 1000* ripples{day}{epoch}{Hptets(1)}.maxthresh;   % in units of std dev

% Load Spikes and MU file - MU file used in case MU cluster not defined in spikes file
% ------------------------------------------------------------------------------------

spikefile = sprintf('%s/%sspikes%02d.mat', directoryname, prefix, day);
if exist(spikefile)==2
    load(spikefile);
end
multifile = sprintf('%s/%smulti%02d.mat', directoryname, prefix, day);
load(multifile);

spikeu_Hp = []; spikeu_PFC=[]; multiu_Hp=[]; multiu_PFC=[];
nHpcells=0; nPFCcells=0;

% Hp DATA
% ------
for ct=1:length(Hptets)
    tet=Hptets(ct);
    % Spikes
    % ------
    if docellcomb==1
        % Get Spike timestamps (in secs - mult by 1000 for ms)
        for i=1:length(spikes{day}{epoch}{tet})
            if ~isempty(spikes{day}{epoch}{tet}{i})
                if ~isempty(spikes{day}{epoch}{tet}{i}.data)
                    cmd=sprintf('spike%d = spikes{day}{epoch}{tet}{%d}.data(:,1)*1000;',i,i); eval(cmd);
                    cmd=sprintf('spikeu_Hp = [spikeu_Hp;spike%d];',i); eval(cmd);
                    nHpcells=nHpcells+1;
                end
            end
        end
    end % end docellcomb for Hptets
    % MULTI
    % -----
    if domulti==1
        if ~isempty(mucluster)
            multiu_Hp = [multiu_Hp;spikes{day}{epoch}{tet}{mucluster}.data(:,1)*1000]; %sec to ms
        else
            if exist('spikes')
                if length(spikes{day}{epoch})>=tet
                    for i=1:length(spikes{day}{epoch}{tet})
                        if ~isempty(spikes{day}{epoch}{tet}{i})
                            cmd=sprintf('currtag = spikes{day}{epoch}{tet}{%d}.tag;',i); eval(cmd);
                            if strcmp(currtag,'MU')
                                multiu_Hp = [multiu_Hp;spikes{day}{epoch}{tet}{i}.data(:,1)*1000]; %sec to ms
                                break
                            end
                        end
                    end
                end
            end
        end
    end % end domulti for Hptets
end % end Hp tets
% If MU is still empty from spikes field after looping tets, use multi file
% --------------------------------------------------------------------------
if domulti==1 && isempty(multiu_Hp)
    for ct=1:length(Hptets)
        tet=Hptets(ct);
        cmd=sprintf('multiu_Hp = [multiu_Hp;multi{day}{epoch}{%d}/10];',tet); eval(cmd); % from 0.1 ms units to ms
    end
end

% PFC Data
% --------
for ct=1:length(PFCtets)
    tet=PFCtets(ct);
    % Spikes
    % ------
    if docellcomb==1
        % Get Spike timestamps (in secs - mult by 1000 for ms)
        if length(spikes{day}{epoch})>=tet     
            for i=1:length(spikes{day}{epoch}{tet})
                if ~isempty(spikes{day}{epoch}{tet}{i})
                    if ~isempty(spikes{day}{epoch}{tet}{i}.data)
                        cmd=sprintf('spike%d = spikes{day}{epoch}{tet}{%d}.data(:,1)*1000;',i,i); eval(cmd);
                        cmd=sprintf('spikeu_PFC = [spikeu_PFC;spike%d];',i); eval(cmd);
                        nPFCcells=nPFCcells+1;
                    end
                end
            end
        end
    end % end docellcomb for PFCtets
    % MULTI
    % -----
    if domulti==1
        if ~isempty(mucluster)
            multiu_PFC = [multiu_PFC;spikes{day}{epoch}{tet}{mucluster}.data(:,1)*1000]; %sec to ms
        else
            if exist('spikes')
                if length(spikes{day}{epoch})>=tet
                    for i=1:length(spikes{day}{epoch}{tet})
                        if ~isempty(spikes{day}{epoch}{tet}{i})
                            cmd=sprintf('currtag = spikes{day}{epoch}{tet}{%d}.tag;',i); eval(cmd);
                            if strcmp(currtag,'MU')
                                multiu_PFC = [multiu_PFC;spikes{day}{epoch}{tet}{i}.data(:,1)*1000]; %sec to ms
                                break
                            end
                        end
                    end
                end
            end
        end
    end % end domulti for PFCtets
end % end PFC tets
% If MU is still empty from spikes field after looping tets, use multi file
% --------------------------------------------------------------------------
if domulti==1 && isempty(multiu_PFC)
    for ct=1:length(PFCtets)
        tet=PFCtets(ct);
        cmd=sprintf('multiu_PFC = [multiu_PFC;multi{day}{epoch}{%d}/10];',tet); eval(cmd); % from 0.1 ms units to ms
    end
end


% EXTRACT EVENTS in MU VECTOR
% ----------------------------
if docellcomb==1
    % Hp
    histvec1_Hp = histc(spikeu_Hp,[min(spikeu_Hp):binsize:max(spikeu_Hp)]);
    histvec1_PFC = histc(spikeu_PFC,[min(spikeu_Hp):binsize:max(spikeu_Hp)]); % Same edges as Hp
    base1_Hp = mean(histvec1_Hp); std1_Hp = std(histvec1_Hp); thresh1_Hp = base1_Hp + 4*std1_Hp; % n s.d.s above mean
    tmpHp1 = extractevents(histvec1_Hp, thresh1_Hp, base1_Hp, 0, 1, 0)';
    %    [event] = extractevents(data, threshold, baseline, min_separation,min_duration(in samples), allowstartend)
    spikeu_Hp_bur.startind = tmpHp1(:,1);
    spikeu_Hp_bur.peakval = tmpHp1(:,3);
    spikeu_Hp_bur.peakind = tmpHp1(:,4);
    spikeu_Hp_bur.starttime = min(spikeu_Hp) + spikeu_Hp_bur.startind*binsize; %What is the time of the burst
    spikeu_Hp_bur.peaktime = min(spikeu_Hp) + spikeu_Hp_bur.peakind*binsize;
    
    %      % PFC
    %      histvec1_PFC = histc(spikeu_PFC,[min(spikeu_PFC):binsize:max(spikeu_PFC));
    %      base1_PFC = mean(histvec1_PFC); std1_PFC = std(histvec1_PFC); thresh1_PFC = base1_PFC + 2*std1_PFC; % 2s.d.s above mean
    %      tmpPFC1 = extractevents(histvec1_PFC, thresh1_PFC, base1_PFC, 0, 1, 0)';
    % %    [event] = extractevents(data, threshold, baseline, min_separation,min_duration(in samples), allowstartend)
    %      spikeu_PFC_bur.startind = tmpPFC(:,1);
    %      spikeu_PFC_bur.peakval = tmpPFC(:,3);
    %      spikeu_PFC_bur.peakind = tmpPFC(:,4);
    %      spikeu_PFC_bur.starttime = spikeu_PFC(spikeu_PFC_bur.startind);
    %      spikeu_PFC_bur.peaktime = spikeu_PFC(spikeu_PFC_bur.peakind);
    
end

if domulti==1
    % Hp
    histvec2_Hp = histc(multiu_Hp,[min(multiu_Hp):binsize:max(multiu_Hp)]);
    histvec2_PFC = histc(multiu_PFC,[min(multiu_Hp):binsize:max(multiu_Hp)]); % Same edges as Hp
    base1_Hp = mean(histvec2_Hp); std1_Hp = std(histvec2_Hp); thresh1_Hp = base1_Hp + 4*std1_Hp; % n s.d.s above mean
    tmpHp1 = extractevents(histvec2_Hp, thresh1_Hp, base1_Hp, 0, 1, 0)';
    %    [event] = extractevents(data, threshold, baseline, min_separation,min_duration(in samples), allowstartend)
    multiu_Hp_bur.startind = tmpHp1(:,1);
    multiu_Hp_bur.peakval = tmpHp1(:,3);
    multiu_Hp_bur.peakind = tmpHp1(:,4);
    multiu_Hp_bur.starttime = min(multiu_Hp) + multiu_Hp_bur.startind*binsize; %What is the time of the burst
    multiu_Hp_bur.peaktime = min(multiu_Hp) + multiu_Hp_bur.peakind*binsize;
end




% ALINGING TO BURST PEAK EVENTS
% -----------------------------

% Align spikes-cellscomb and MU to bursts

if docellcomb==1
    %     cntbur=0;
    %     for i=2:length(spikeu_Hp_bur.peaktime)-1
    %         i;
    %         cntbur=cntbur+1;
    %         currbur = spikeu_Hp_bur.peaktime(i); bursize_cell(cntbur) = spikeu_Hp_bur.peakval(i);
    %         % Hp
    %         currspks =  spikeu_Hp(find( (spikeu_Hp>=(currbur-pret)) & (spikeu_Hp<=(currbur+postt)) ));
    %         currspks = currspks-(currbur-pret);
    %         histspks = histc(currspks,[0:binsize:pret+postt]);
    %         bur_spks_cell_Hp{cntbur}=currspks;
    %         bur_spkshist_cell_Hp(cntbur,:) = histspks;
    %         % Smooth over binsize*3 ms (30ms for 10ms; 15ms for 5 ms)
    %         nstd = round(binsize*3/binsize); g1 = gaussian(nstd, 5*nstd+1);
    %         histspks = smoothvect(histspks, g1);
    %         bur_spkshist_cell_Hp(cntbur,:) = histspks;
    %         % PFC
    %         currspks =  spikeu_PFC(find( (spikeu_PFC>=(currbur-pret)) & (spikeu_PFC<=(currbur+postt)) ));
    %         currspks = currspks-(currbur-pret);
    %         histspks = histc(currspks,[0:binsize:pret+postt]);
    %         bur_spks_cell_PFC{cntbur}=currspks;
    %         bur_spkshist_cell_PFC(cntbur,:) = histspks;
    %         % Smooth over binsize*3 ms (30ms for 10ms; 15ms for 5 ms)
    %         nstd = round(binsize*3/binsize); g1 = gaussian(nstd, 5*nstd+1);
    %         histspks = smoothvect(histspks, g1);
    %         bur_spkshist_cell_PFC(cntbur,:) = histspks;
    %     end
    
    % Alternatively - use the vector already made
    cntbur=0;
    for i=5:length(spikeu_Hp_bur.peakind)-5
        i;
        cntbur=cntbur+1;
        currbur = spikeu_Hp_bur.peakind(i); bursize_cell(cntbur) = spikeu_Hp_bur.peakval(i);
        % Hp
        histspks = histvec1_Hp(spikeu_Hp_bur.peakind(i)-pret/binsize:spikeu_Hp_bur.peakind(i)+postt/binsize);
        bur_spkshist_cell_Hp(cntbur,:) = histspks;
        % Smooth over binsize*3 ms (30ms for 10ms; 15ms for 5 ms)
        nstd = round(binsize*3/binsize); g1 = gaussian(nstd, 5*nstd+1);
        histspks = smoothvect(histspks, g1);
        bur_spkshist_cell_Hp(cntbur,:) = histspks;
        % PFC
        histspks = histvec1_PFC(spikeu_Hp_bur.peakind(i)-pret/binsize:spikeu_Hp_bur.peakind(i)+postt/binsize);
        bur_spkshist_cell_PFC(cntbur,:) = histspks;
        % Smooth over binsize*3 ms (30ms for 10ms; 15ms for 5 ms)
        nstd = round(binsize*3/binsize); g1 = gaussian(nstd, 5*nstd+1);
        histspks = smoothvect(histspks, g1);
        bur_spkshist_cell_PFC(cntbur,:) = histspks;
    end
    
    % Convert to Z-scores
    % -------------------
    cellmean_Hp = mean(bur_spkshist_cell_Hp,2); cellstd_Hp = std(bur_spkshist_cell_Hp,[],2);
    cellmean_PFC = mean(bur_spkshist_cell_PFC,2); cellstd_PFC = std(bur_spkshist_cell_PFC,[],2);
    x=find(cellstd_Hp==0 | cellstd_PFC==0); % Remove rows which are all zero
    bur_spkshist_cell_Hp(x,:)=[]; cellmean_Hp(x)=[]; cellstd_Hp(x)=[]; %bur_spks_cell{x}=[];
    bur_spkshist_cell_PFC(x,:)=[]; cellmean_PFC(x)=[]; cellstd_PFC(x)=[]; bursize_cell(x) = [];
    % 1a) Hp
    cellmean_mat = repmat(cellmean_Hp,1,size(bur_spkshist_cell_Hp,2)); cellstd_mat = repmat(cellstd_Hp,1,size(bur_spkshist_cell_Hp,2));
    bur_spkshist_cellZ_Hp = (bur_spkshist_cell_Hp-cellmean_mat)./cellstd_mat;
    % 1b) PFC
    cellmean_mat = repmat(cellmean_PFC,1,size(bur_spkshist_cell_PFC,2)); cellstd_mat = repmat(cellstd_PFC,1,size(bur_spkshist_cell_PFC,2));
    bur_spkshist_cellZ_PFC = (bur_spkshist_cell_PFC-cellmean_mat)./cellstd_mat;
    
    %save temp
    
    % % Sort by MU rate in -100 to 100 ms for Hp, and apply same order to PFC/sort PFC by itself
    % % ----------------------------------------------------------------------
    cellresp = sum(bur_spkshist_cellZ_Hp(:,bins_resp),2);
    [cellresp_sort,cellsortidx] = sort(cellresp);
    bur_spkshist_cellZsort_Hp = bur_spkshist_cellZ_Hp(cellsortidx,:);
    %cellresp = sum(bur_spkshist_cellZ_PFC(:,bins_resp),2);
    %[cellresp_sort,cellsortidx] = sort(cellresp);
    bur_spkshist_cellZsort_PFC = bur_spkshist_cellZ_PFC(cellsortidx,:);
    
    % -------------------------------------------------
    % Smooth Matrix
    % -------------------------------------------------
    for i=1:length(bur_spkshist_cellZsort_Hp)
        winst = i-smwin/2; winend = i+smwin/2;
        if winst<1, winst=1; end
        if winend>length(bur_spkshist_cellZsort_Hp), winend = length(bur_spkshist_cellZsort_Hp); end
        bur_spkshist_cellZsort_Hp(i,:) = mean(bur_spkshist_cellZsort_Hp(winst:winend,:));
        bur_spkshist_cellZsort_PFC(i,:) = mean(bur_spkshist_cellZsort_PFC(winst:winend,:));
    end
    % TRIM ENDS
    bur_spkshist_cellZsort_Hp = bur_spkshist_cellZsort_Hp(:,50/binsize:(pret+postt-50)/binsize);
    bur_spkshist_cellZsort_PFC = bur_spkshist_cellZsort_PFC(:,50/binsize:(pret+postt-50)/binsize);
    
end


if domulti==1
    %     cntbur=0;
    %     for i=2:length(spikeu_Hp_bur.peaktime)-1
    %         i;
    %         cntbur=cntbur+1;
    %         currbur = multiu_Hp_bur.peaktime(i); bursize_cell(cntbur) = multiu_Hp_bur.peakval(i);
    %         % Hp
    %         bursize_multi(cntbur) = bur_sizes(i);
    %         currspks =  multiu_Hp(find( (multiu_Hp>=(currbur-pret)) & (multiu_Hp<=(currbur+postt)) ));
    %         currspks = currspks-(currbur-pret);
    %         histspks = histc(currspks,[0:binsize:pret+postt]);
    %         bur_spks_multi_Hp{cntbur}=currspks;
    %         bur_spkshist_multi_Hp(cntbur,:) = histspks;
    %         % Smooth over binsize*3 ms (30ms for 10ms; 15ms for 5 ms) or binsize*2 ms
    %         nstd = round(binsize*3/binsize); g1 = gaussian(nstd, 5*nstd+1);
    %         histspks = smoothvect(histspks, g1);
    %         bur_spkshist_multi_Hp(cntbur,:) = histspks;
    %         % PFC
    %         currspks =  multiu_PFC(find( (multiu_PFC>=(currbur-pret)) & (multiu_PFC<=(currbur+postt)) ));
    %         currspks = currspks-(currbur-pret);
    %         histspks = histc(currspks,[0:binsize:pret+postt]);
    %         bur_spks_multi_PFC{cntbur}=currspks;
    %         bur_spkshist_multi_PFC(cntbur,:) = histspks;
    %         % Smooth over binsize*3 ms (30ms for 10ms; 15ms for 5 ms) or binsize*2 ms
    %         nstd = round(binsize*3/binsize); g1 = gaussian(nstd, 5*nstd+1);
    %         histspks = smoothvect(histspks, g1);
    %         bur_spkshist_multi_PFC(cntbur,:) = histspks;
    %     end
    
    % Alternatively - use the vector already made
    cntbur=0;
    for i=5:length(multiu_Hp_bur.peakind)-5
        i;
        cntbur=cntbur+1;
        currbur = multiu_Hp_bur.peakind(i); bursize_multi(cntbur) = multiu_Hp_bur.peakval(i);
        % Hp
        histspks = histvec2_Hp(multiu_Hp_bur.peakind(i)-pret/binsize:multiu_Hp_bur.peakind(i)+postt/binsize);
        bur_spkshist_multi_Hp(cntbur,:) = histspks;
        % Smooth over binsize*3 ms (30ms for 10ms; 15ms for 5 ms)
        nstd = round(binsize*3/binsize); g1 = gaussian(nstd, 5*nstd+1);
        histspks = smoothvect(histspks, g1);
        bur_spkshist_multi_Hp(cntbur,:) = histspks;
        % PFC
        histspks = histvec2_PFC(multiu_Hp_bur.peakind(i)-pret/binsize:multiu_Hp_bur.peakind(i)+postt/binsize);
        bur_spkshist_multi_PFC(cntbur,:) = histspks;
        % Smooth over binsize*3 ms (30ms for 10ms; 15ms for 5 ms)
        nstd = round(binsize*3/binsize); g1 = gaussian(nstd, 5*nstd+1);
        histspks = smoothvect(histspks, g1);
        bur_spkshist_multi_PFC(cntbur,:) = histspks;
    end
    
    
    % Convert to Z-scores
    % -------------------
    multimean_Hp = mean(bur_spkshist_multi_Hp,2); multistd_Hp = std(bur_spkshist_multi_Hp,[],2);
    multimean_PFC = mean(bur_spkshist_multi_PFC,2); multistd_PFC = std(bur_spkshist_multi_PFC,[],2);
    x=find(multistd_Hp==0 | multistd_PFC==0); % Remove rows which are all zero
    bur_spkshist_multi_Hp(x,:)=[]; multimean_Hp(x)=[]; multistd_Hp(x)=[]; %bur_spks_multi{x}=[];
    bur_spkshist_multi_PFC(x,:)=[]; multimean_PFC(x)=[]; multistd_PFC(x)=[]; bursize_multi(x) = [];
    % 2a) Hp
    multimean_mat = repmat(multimean_Hp,1,size(bur_spkshist_multi_Hp,2)); multistd_mat = repmat(multistd_Hp,1,size(bur_spkshist_multi_Hp,2));
    bur_spkshist_multiZ_Hp = (bur_spkshist_multi_Hp-multimean_mat)./multistd_mat;
    % 2b) PFC
    multimean_mat = repmat(multimean_PFC,1,size(bur_spkshist_multi_PFC,2)); multistd_mat = repmat(multistd_PFC,1,size(bur_spkshist_multi_PFC,2));
    bur_spkshist_multiZ_PFC = (bur_spkshist_multi_PFC-multimean_mat)./multistd_mat;
    
    % % Sort by MU rate in -100 to 100 ms for Hp, and apply same order to PFC/sort PFC by itself
    % % ----------------------------------------------------------------------
    multiresp = sum(bur_spkshist_multiZ_Hp(:,bins_resp),2);
    [multiresp_sort,multisortidx] = sort(multiresp);
    bur_spkshist_multiZsort_Hp = bur_spkshist_multiZ_Hp(multisortidx,:);
    %multiresp = sum(bur_spkshist_multiZ_PFC(:,bins_resp),2);
    %[multiresp_sort,multisortidx] = sort(multiresp);
    bur_spkshist_multiZsort_PFC = bur_spkshist_multiZ_PFC(multisortidx,:);
    
    % -------------------------------------------------
    % Smooth Matrix
    % -------------------------------------------------
    for i=1:length(bur_spkshist_multiZsort_Hp)
        winst = i-smwin/2; winend = i+smwin/2;
        if winst<1, winst=1; end
        if winend>length(bur_spkshist_multiZsort_Hp), winend = length(bur_spkshist_multiZsort_Hp); end
        bur_spkshist_multiZsort_Hp(i,:) = mean(bur_spkshist_multiZsort_Hp(winst:winend,:));
        bur_spkshist_multiZsort_PFC(i,:) = mean(bur_spkshist_multiZsort_PFC(winst:winend,:));
    end
    
    % TRIM ENDS
    bur_spkshist_multiZsort_Hp = bur_spkshist_multiZsort_Hp(:,50/binsize:(pret+postt-50)/binsize);
    bur_spkshist_multiZsort_PFC = bur_spkshist_multiZsort_PFC(:,50/binsize:(pret+postt-50)/binsize);
    %bur_spkshist_multiZsort_Hp(:,1:50/binsize) = bur_spkshist_multiZsort_Hp(:,50/binsize+1:100/binsize);
    %bur_spkshist_multiZsort_Hp(:,(pret+postt-50)/binsize:end) = bur_spkshist_multiZsort_Hp(:,(pret+postt-100)/binsize-2:(pret+postt-50)/binsize-1);
    %bur_spkshist_multiZsort_PFC(:,1:50/binsize) = bur_spkshist_multiZsort_PFC(:,50/binsize+1:100/binsize);
    %bur_spkshist_multiZsort_PFC(:,(pret+postt-50)/binsize:end) = bur_spkshist_multiZsort_PFC(:,(pret+postt-100)/binsize-2:(pret+postt-50)/binsize-1);
end


% Update pret and postt due to end trims
pret=500; postt=500;







%% ------------------------------------------------
% PLOT
% -------------------------------------------------

if docellcomb == 1
    
    %% Imagesc-Hp Cell
    % -------------
    
    figure; hold on;
    %redimscreen_figforppt1;
    set(gcf,'Position',[0 Screen(4)*0.55 Screen(3)*0.2 Screen(4)*0.4])
    subplot(3,1,[1 2]); hold on;
    imagesc(flipud(bur_spkshist_cellZsort_Hp));
    title(['Hp CellComb Firing Z aligned to bursts-' num2str(binsize) 'ms bins'],...
        'FontSize',20,'Fontweight','normal');
    %ylabel('Stim Amp / Laser Power','FontSize',20,'Fontweight','normal');
    %axis([0 ((pret+postt)/binsize)+2 0.5 nranges+0.5]);
    set(gca,'XLim',[1 ((pret+postt)/binsize)+1]);
    set(gca,'YLim',[0 size(bur_spkshist_cellZsort_Hp,1)]);
    %xlabel('Time(ms)','FontSize',20,'Fontweight','normal');
    
    ypts = 0:1:size(bur_spkshist_cellZsort_Hp,1);
    xpts = ((pret/binsize)+1)*ones(size(ypts));
    
    % Plot Line at 0 ms - Onset of stimulation
    plot(xpts , ypts, 'k--','Linewidth',3);
    % Plot lines at 100ms and 2000ms
    xpts = ((pret+plotline1+binsize)/binsize)*ones(size(ypts));
    plot(xpts , ypts, 'k--','Linewidth',2);
    if ~isempty(plotline2)
        xpts = ((pret+plotline2-binsize)/binsize)*ones(size(ypts));
        plot(xpts , ypts, 'k--','Linewidth',2);
    end
    % set(gca,'xtick',[-1,(pret/binsize)-2,(pret+postt)/binsize-1],'xticklabel',{num2str([-pret,0,postt]')},...
    %      'FontSize',20,'Fontweight','normal');
    % Ticks every 100ms
    set(gca,'xtick',[1:(100/binsize):((pret+postt)/binsize)+1],'xticklabel',{num2str([-pret:100:postt]')},...
        'FontSize',14,'Fontweight','normal');
    
    subplot(3,1,3); hold on;
    plot(mean(bur_spkshist_cellZsort_Hp),'Linewidth',3);
    set(gca,'XLim',[1 ((pret+postt)/binsize)+1]);
    xlabel('Time(ms)','FontSize',20,'Fontweight','normal');
    set(gca,'xtick',[1:(100/binsize):((pret+postt)/binsize)+1],'xticklabel',{num2str([-pret:100:postt]')},...
        'FontSize',14,'Fontweight','normal');
    set(gca,'YLim',[min(mean(bur_spkshist_cellZsort_Hp)) max(mean(bur_spkshist_cellZsort_Hp))+0.5]);
    ypts = min(mean(bur_spkshist_cellZsort_Hp)):0.1:max(mean(bur_spkshist_cellZsort_Hp))+0.5;
    xpts = ((pret/binsize)+1)*ones(size(ypts));
    % Plot Line at 0 ms - Onset of stimulation
    plot(xpts , ypts, 'k--','Linewidth',2);
    % Plot lines at 100ms and 2000ms
    xpts = ((pret+plotline1+binsize)/binsize)*ones(size(ypts));
    plot(xpts , ypts, 'k--','Linewidth',1);
    
    if saveg1==1,
        figfile = [figdir,prefix,'_Day',num2str(day),'_Tet',num2str(tetu),'_HpCellMatrix'];
        print('-dpdf', figfile);
        print('-djpeg', figfile);
        saveas(gcf,figfile,'fig');
    end
    
    
    %% Imagesc-PFC Cell
    % -------------
    
    figure; hold on;
    %redimscreen_figforppt1;
    set(gcf,'Position',[Screen(3)*0.205 Screen(4)*0.55 Screen(3)*0.2 Screen(4)*0.4])
    subplot(3,1,[1 2]); hold on;
    imagesc(flipud(bur_spkshist_cellZsort_PFC));
    title(['PFC CellComb Firing Z aligned to bursts-' num2str(binsize) 'ms bins'],...
        'FontSize',20,'Fontweight','normal');
    %ylabel('Stim Amp / Laser Power','FontSize',20,'Fontweight','normal');
    %axis([0 ((pret+postt)/binsize)+2 0.5 nranges+0.5]);
    set(gca,'XLim',[1 ((pret+postt)/binsize)+1]);
    set(gca,'YLim',[0 size(bur_spkshist_cellZsort_PFC,1)]);
    %xlabel('Time(ms)','FontSize',20,'Fontweight','normal');
    
    ypts = 0:1:size(bur_spkshist_cellZsort_PFC,1);
    xpts = ((pret/binsize)+1)*ones(size(ypts));
    
    % Plot Line at 0 ms - Onset of stimulation
    plot(xpts , ypts, 'k--','Linewidth',3);
    % Plot lines at 100ms and 2000ms
    xpts = ((pret+plotline1+binsize)/binsize)*ones(size(ypts));
    plot(xpts , ypts, 'k--','Linewidth',2);
    if ~isempty(plotline2)
        xpts = ((pret+plotline2-binsize)/binsize)*ones(size(ypts));
        plot(xpts , ypts, 'k--','Linewidth',2);
    end
    % set(gca,'xtick',[-1,(pret/binsize)-2,(pret+postt)/binsize-1],'xticklabel',{num2str([-pret,0,postt]')},...
    %      'FontSize',20,'Fontweight','normal');
    % Ticks every 100ms
    set(gca,'xtick',[1:(100/binsize):((pret+postt)/binsize)+1],'xticklabel',{num2str([-pret:100:postt]')},...
        'FontSize',14,'Fontweight','normal');
    
    subplot(3,1,3); hold on;
    plot(mean(bur_spkshist_cellZsort_PFC),'Linewidth',3);
    set(gca,'XLim',[1 ((pret+postt)/binsize)+1]);
    xlabel('Time(ms)','FontSize',20,'Fontweight','normal');
    set(gca,'xtick',[1:(100/binsize):((pret+postt)/binsize)+1],'xticklabel',{num2str([-pret:100:postt]')},...
        'FontSize',14,'Fontweight','normal');
    set(gca,'YLim',[min(mean(bur_spkshist_cellZsort_PFC)) max(mean(bur_spkshist_cellZsort_PFC))+0.5]);
    ypts = min(mean(bur_spkshist_cellZsort_PFC)):0.1:max(mean(bur_spkshist_cellZsort_PFC))+0.5;
    xpts = ((pret/binsize)+1)*ones(size(ypts));
    % Plot Line at 0 ms - Onset of stimulation
    plot(xpts , ypts, 'k--','Linewidth',2);
    % Plot lines at 100ms and 2000ms
    xpts = ((pret+plotline1+binsize)/binsize)*ones(size(ypts));
    plot(xpts , ypts, 'k--','Linewidth',1);
    
    if saveg1==1,
        figfile = [figdir,prefix,'_Day',num2str(day),'_Tet',num2str(tetu),'_PFCCellMatrix'];
        print('-dpdf', figfile);
        print('-djpeg', figfile);
        saveas(gcf,figfile,'fig');
    end
    
end

if domulti==1
    
    %% Imagesc-Hp Multi
    % -------------
    
    figure; hold on;
    %redimscreen_figforppt1;
    set(gcf,'Position',[0 Screen(4)*0.03 Screen(3)*0.2 Screen(4)*0.4])
    subplot(3,1,[1 2]); hold on;
    imagesc(flipud(bur_spkshist_multiZsort_Hp));
    title(['Hp Multi Firing Z aligned to bursts-' num2str(binsize) 'ms bins'],...
        'FontSize',20,'Fontweight','normal');
    %ylabel('Stim Amp / Laser Power','FontSize',20,'Fontweight','normal');
    %axis([0 ((pret+postt)/binsize)+2 0.5 nranges+0.5]);
    set(gca,'XLim',[1 ((pret+postt)/binsize)+1]);
    set(gca,'YLim',[0 size(bur_spkshist_multiZsort_Hp,1)])
    %xlabel('Time(ms)','FontSize',20,'Fontweight','normal');
    
    ypts = 0:1:size(bur_spkshist_multiZsort_Hp,1);
    xpts = ((pret/binsize)+1)*ones(size(ypts));
    
    % Plot Line at 0 ms - Onset of stimulation
    plot(xpts , ypts, 'k--','Linewidth',3);
    % Plot lines at 100ms and 2000ms
    xpts = ((pret+plotline1+binsize)/binsize)*ones(size(ypts));
    plot(xpts , ypts, 'k--','Linewidth',2);
    if ~isempty(plotline2)
        xpts = ((pret+plotline2-binsize)/binsize)*ones(size(ypts));
        plot(xpts , ypts, 'k--','Linewidth',2);
    end
    % set(gca,'xtick',[-1,(pret/binsize)-2,(pret+postt)/binsize-1],'xticklabel',{num2str([-pret,0,postt]')},...
    %      'FontSize',20,'Fontweight','normal');
    % Ticks every 100ms
    set(gca,'xtick',[1:(100/binsize):((pret+postt)/binsize)+1],'xticklabel',{num2str([-pret:100:postt]')},...
        'FontSize',14,'Fontweight','normal');
    
    subplot(3,1,3); hold on;
    plot(mean(bur_spkshist_multiZsort_Hp),'Linewidth',3);
    set(gca,'XLim',[1 ((pret+postt)/binsize)+1]);
    xlabel('Time(ms)','FontSize',20,'Fontweight','normal');
    set(gca,'xtick',[1:(100/binsize):((pret+postt)/binsize)+1],'xticklabel',{num2str([-pret:100:postt]')},...
        'FontSize',14,'Fontweight','normal');
    set(gca,'YLim',[min(mean(bur_spkshist_multiZsort_Hp)) max(mean(bur_spkshist_multiZsort_Hp))+0.5]);
    ypts = min(mean(bur_spkshist_multiZsort_Hp)):0.1:max(mean(bur_spkshist_multiZsort_Hp))+0.5;
    xpts = ((pret/binsize)+1)*ones(size(ypts));
    % Plot Line at 0 ms - Onset of stimulation
    plot(xpts , ypts, 'k--','Linewidth',2);
    % Plot lines at 100ms and 2000ms
    xpts = ((pret+plotline1+binsize)/binsize)*ones(size(ypts));
    plot(xpts , ypts, 'k--','Linewidth',1);
    
    
    
    if saveg1==1,
        figfile = [figdir,prefix,'_Day',num2str(day),'_Tet',num2str(tetu),'_HpMultiMatrix'];
        print('-dpdf', figfile);
        print('-djpeg', figfile);
        saveas(gcf,figfile,'fig');
    end
    
    
    
    
    %% Imagesc-PFC Multi
    % -------------
    
    figure; hold on;
    %redimscreen_figforppt1;
    set(gcf,'Position',[Screen(3)*0.205 Screen(4)*0.03 Screen(3)*0.2 Screen(4)*0.4])
    subplot(3,1,[1 2]); hold on;
    imagesc(flipud(bur_spkshist_multiZsort_PFC));
    title(['PFC Multi Firing Z aligned to bursts-' num2str(binsize) 'ms bins'],...
        'FontSize',20,'Fontweight','normal');
    %ylabel('Stim Amp / Laser Power','FontSize',20,'Fontweight','normal');
    %axis([0 ((pret+postt)/binsize)+2 0.5 nranges+0.5]);
    set(gca,'XLim',[1 ((pret+postt)/binsize)+1]);
    set(gca,'YLim',[0 size(bur_spkshist_multiZsort_PFC,1)])
    %xlabel('Time(ms)','FontSize',20,'Fontweight','normal');
    
    ypts = 0:1:size(bur_spkshist_multiZsort_PFC,1);
    xpts = ((pret/binsize)+1)*ones(size(ypts));
    
    % Plot Line at 0 ms - Onset of stimulation
    plot(xpts , ypts, 'k--','Linewidth',3);
    % Plot lines at 100ms and 2000ms
    xpts = ((pret+plotline1+binsize)/binsize)*ones(size(ypts));
    plot(xpts , ypts, 'k--','Linewidth',2);
    if ~isempty(plotline2)
        xpts = ((pret+plotline2-binsize)/binsize)*ones(size(ypts));
        plot(xpts , ypts, 'k--','Linewidth',2);
    end
    % set(gca,'xtick',[-1,(pret/binsize)-2,(pret+postt)/binsize-1],'xticklabel',{num2str([-pret,0,postt]')},...
    %      'FontSize',20,'Fontweight','normal');
    % Ticks every 100ms
    set(gca,'xtick',[1:(100/binsize):((pret+postt)/binsize)+1],'xticklabel',{num2str([-pret:100:postt]')},...
        'FontSize',14,'Fontweight','normal');
    
    subplot(3,1,3); hold on;
    plot(mean(bur_spkshist_multiZsort_PFC),'Linewidth',3);
    set(gca,'XLim',[1 ((pret+postt)/binsize)+1]);
    xlabel('Time(ms)','FontSize',20,'Fontweight','normal');
    set(gca,'xtick',[1:(100/binsize):((pret+postt)/binsize)+1],'xticklabel',{num2str([-pret:100:postt]')},...
        'FontSize',14,'Fontweight','normal');
    set(gca,'YLim',[min(mean(bur_spkshist_multiZsort_PFC)) max(mean(bur_spkshist_multiZsort_PFC))+0.5]);
    ypts = min(mean(bur_spkshist_multiZsort_PFC)):0.1:max(mean(bur_spkshist_multiZsort_PFC))+0.5;
    xpts = ((pret/binsize)+1)*ones(size(ypts));
    % Plot Line at 0 ms - Onset of stimulation
    plot(xpts , ypts, 'k--','Linewidth',2);
    % Plot lines at 100ms and 2000ms
    xpts = ((pret+plotline1+binsize)/binsize)*ones(size(ypts));
    plot(xpts , ypts, 'k--','Linewidth',1);
    
    if saveg1==1,
        figfile = [figdir,prefix,'_Day',num2str(day),'_Tet',num2str(tetu),'_PFCMultiMatrix'];
        print('-dpdf', figfile);
        print('-djpeg', figfile);
        saveas(gcf,figfile,'fig');
    end
    
end
i=1;
%set(gcf,'Position',[Screen(3)*0.205 Screen(4)*0.03 Screen(3)*0.2 Screen(4)*0.4])








%figfile = ['Figures/JW7_d7t1c1_inh_hist100_1sec'];saveas(gcf,figfile,'fig');print('-djpeg', figfile);print('-dpdf', figfile);

































%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



