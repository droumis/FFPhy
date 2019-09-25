function [basefr0] = sj_HPexpt_burstaligntorip(prefix, day, epoch, riptet, Hptets, PFCtets, docellcomb, domulti, mucluster, binsize, saveg1, dospeed)
% Shantanu - Aug 2012. From RippleDisruption / sj_multiunitandspike_align
% Aligns bursts around ripples in increasing order of size for Hp and PFC. Similar to Wierzynski, 2009 - Fig. 8
% eg. sj_HPexpt_burstaligntorip('HPa',2,5,1,[1 4],[15 18],1,1,[],10,0);
% sj_HPexpt_burstaligntorip('HPb',1,7,4,[1 3 4 6 18],[8 9 10 12],1,1,[],10,0);

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
    riptet=1; %
end
if nargin<5,
    Hptets=[]; %
end
if nargin<6,
    PFCtets=[]; %
end
if nargin<7,
    docellcomb=0; %
end
if nargin<8,
    domulti=1; %
end
if nargin<9,
    mucluster=[]; % If no MU cluster is specified, then will look for Mu tag in spikes. If does not exist, will use MU file
end
if nargin<10,
    binsize=5;
end
if nargin<11
    saveg1=0;
end
if nargin<12
    dospeed=[];
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
% currdir = pwd;
% cd(rawdir);
% dayfolders = dir;
% daystr = sprintf('%02d', day);
% for i = 3:length(dayfolders)
%     if dayfolders(i).isdir
%         if strcmp(daystr,dayfolders(i).name(1:2))
%             disp(upper(dayfolders(i).name))
%             cd(dayfolders(i).name);
%             load times;
%         end
%     end
% end
% cd(currdir);
% % Now Getting Range directly from times file
% userange = ranges(epoch+1,:); % 1st row is allepochs. So you need +1
% usename = names{epoch+1}(end-15:end);
% amp=70; amp1=amp;
% nranges=1;

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
% SHOULD SWITCH THIS TO ALLTET. If alltet does not exist, then use given tet/tets.
ripfile = sprintf('%s/%sripples%02d.mat', directoryname, prefix, day);
load(ripfile);
rip_starttime=[]; rip_sizes=[];
for i=1:length(riptet)
    currriptet=riptet(i);
    rip_starttime = [rip_starttime; 1000* ripples{day}{epoch}{currriptet}.starttime];   % in msec
    rip_sizes = [rip_sizes; ripples{day}{epoch}{currriptet}.maxthresh];   % in units of std dev
end
[rip_starttime,sortidx] = sort(rip_starttime);
rip_sizes = rip_sizes(sortidx);
% Find ripples within xxms of each other
%remidx = find(diff(rip_starttime)<=50)+1;
%rip_starttime(remidx)=[]; rip_sizes(remidx)=[];



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
                %tet, i
                %if tet==4 & i==7, keyboard; end
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
if isempty(multiu_Hp)
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
        for i=1:length(spikes{day}{epoch}{tet})
            if ~isempty(spikes{day}{epoch}{tet}{i})
                if ~isempty(spikes{day}{epoch}{tet}{i}.data)
                    cmd=sprintf('spike%d = spikes{day}{epoch}{tet}{%d}.data(:,1)*1000;',i,i); eval(cmd);
                    cmd=sprintf('spikeu_PFC = [spikeu_PFC;spike%d];',i); eval(cmd);
                    nPFCcells=nPFCcells+1;
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
if isempty(multiu_PFC)
    for ct=1:length(PFCtets)
        tet=PFCtets(ct);
        cmd=sprintf('multiu_PFC = [multiu_PFC;multi{day}{epoch}{%d}/10];',tet); eval(cmd); % from 0.1 ms units to ms
    end
end


% ALINGING TO RIPPLES
% ----------------------

% Align spikes-cellscomb and MU to ripples
cntrip=0;
for i=2:length(rip_starttime)-1
    i;
    cntrip=cntrip+1;
    currrip = rip_starttime(i);
    ripsize_cell(cntrip) = rip_sizes(i); ripsize_multi(cntrip) = rip_sizes(i);
    if docellcomb==1
        % Hp
        currspks =  spikeu_Hp(find( (spikeu_Hp>=(currrip-pret)) & (spikeu_Hp<=(currrip+postt)) ));
        currspks = currspks-(currrip);
        histspks = histc(currspks,[-pret:binsize:postt]);
        rip_spks_cell_Hp{cntrip}=currspks;
        rip_spkshist_cell_Hp(cntrip,:) = histspks;
        % Smooth over binsize*3 ms (30ms for 10ms; 15ms for 5 ms)
        nstd = round(binsize*3/binsize); g1 = gaussian(nstd, 5*nstd+1);
        histspks = smoothvect(histspks, g1);
        rip_spkshist_cell_Hp(cntrip,:) = histspks;
        % PFC
        currspks =  spikeu_PFC(find( (spikeu_PFC>=(currrip-pret)) & (spikeu_PFC<=(currrip+postt)) ));
        currspks = currspks-(currrip);
        histspks = histc(currspks,[-pret:binsize:postt]);
        rip_spks_cell_PFC{cntrip}=currspks;
        rip_spkshist_cell_PFC(cntrip,:) = histspks;
        % Smooth over binsize*3 ms (30ms for 10ms; 15ms for 5 ms)
        nstd = round(binsize*3/binsize); g1 = gaussian(nstd, 5*nstd+1);
        histspks = smoothvect(histspks, g1);
        rip_spkshist_cell_PFC(cntrip,:) = histspks;
    end
    
    if domulti==1
        % Hp
        currspks =  multiu_Hp(find( (multiu_Hp>=(currrip-pret)) & (multiu_Hp<=(currrip+postt)) ));
        currspks = currspks-(currrip);
        histspks = histc(currspks,[-pret:binsize:postt]);
        rip_spks_multi_Hp{cntrip}=currspks;
        rip_spkshist_multi_Hp(cntrip,:) = histspks;
        % Smooth over binsize*3 ms (30ms for 10ms; 15ms for 5 ms) or binsize*2 ms
        nstd = round(binsize*3/binsize); g1 = gaussian(nstd, 5*nstd+1);
        histspks = smoothvect(histspks, g1);
        rip_spkshist_multi_Hp(cntrip,:) = histspks;
        % PFC
        currspks =  multiu_PFC(find( (multiu_PFC>=(currrip-pret)) & (multiu_PFC<=(currrip+postt)) ));
        currspks = currspks-(currrip);
        histspks = histc(currspks,[-pret:binsize:postt]);
        rip_spks_multi_PFC{cntrip}=currspks;
        rip_spkshist_multi_PFC(cntrip,:) = histspks;
        % Smooth over binsize*3 ms (30ms for 10ms; 15ms for 5 ms) or binsize*2 ms
        nstd = round(binsize*3/binsize); g1 = gaussian(nstd, 5*nstd+1);
        histspks = smoothvect(histspks, g1);
        rip_spkshist_multi_PFC(cntrip,:) = histspks;
    end
end


% Convert to Z-scores
% -------------------
if docellcomb==1
    % 1) For cellspks
    cellmean_Hp = mean(rip_spkshist_cell_Hp,2); cellstd_Hp = std(rip_spkshist_cell_Hp,[],2);
    cellmean_PFC = mean(rip_spkshist_cell_PFC,2); cellstd_PFC = std(rip_spkshist_cell_PFC,[],2);
    x=find(cellstd_Hp==0 | cellstd_PFC==0); % Remove rows which are all zero / 
    % OR can set all these rows to zero in Z-score: Make means for these rows 0, and make std any n (say 1). Z-score wil then be 0 
    rip_spkshist_cell_Hp(x,:)=[]; cellmean_Hp(x)=[]; cellstd_Hp(x)=[]; %rip_spks_cell{x}=[];
    rip_spkshist_cell_PFC(x,:)=[]; cellmean_PFC(x)=[]; cellstd_PFC(x)=[]; ripsize_cell(x) = [];
    % 1a) Hp
    cellmean_mat = repmat(cellmean_Hp,1,size(rip_spkshist_cell_Hp,2)); cellstd_mat = repmat(cellstd_Hp,1,size(rip_spkshist_cell_Hp,2));
    rip_spkshist_cellZ_Hp = (rip_spkshist_cell_Hp-cellmean_mat)./cellstd_mat;
    % 1b) PFC
    cellmean_mat = repmat(cellmean_PFC,1,size(rip_spkshist_cell_PFC,2)); cellstd_mat = repmat(cellstd_PFC,1,size(rip_spkshist_cell_PFC,2));
    rip_spkshist_cellZ_PFC = (rip_spkshist_cell_PFC-cellmean_mat)./cellstd_mat;
    
    % A) SORT BY RIPPLE SIZE. Do either A or B
    % -----------------------------------------
    % [cellresp_sort,cellsortidx] = sort(ripsize_cell);
    % rip_spkshist_cellZsort_Hp = rip_spkshist_cellZ_Hp(cellsortidx,:);
    % rip_spkshist_cellZsort_PFC = rip_spkshist_cellZ_PFC(cellsortidx,:);
    
    % % B )Sort by MU rate in -100 to 100 ms for Hp, and apply same order to PFC/sort PFC by itself
    % % ----------------------------------------------------------------------
    cellresp = sum(rip_spkshist_cellZ_Hp(:,bins_resp),2);
    [cellresp_sort,cellsortidx] = sort(cellresp);
    rip_spkshist_cellZsort_Hp = rip_spkshist_cellZ_Hp(cellsortidx,:);
    %cellresp = sum(rip_spkshist_cellZ_PFC(:,bins_resp),2);
    %[cellresp_sort,cellsortidx] = sort(cellresp);
    rip_spkshist_cellZsort_PFC = rip_spkshist_cellZ_PFC(cellsortidx,:);
    
    % -------------------------------------------------
    % Smooth Matrix
    % -------------------------------------------------
    
    for i=1:size(rip_spkshist_cellZsort_Hp,1)
        winst = i-smwin/2; winend = i+smwin/2;
        if winst<1, winst=1; end
        if winend>size(rip_spkshist_cellZsort_Hp,1), winend = size(rip_spkshist_cellZsort_Hp,1); end
        rip_spkshist_cellZsort_Hp(i,:) = mean(rip_spkshist_cellZsort_Hp(winst:winend,:));
        rip_spkshist_cellZsort_PFC(i,:) = mean(rip_spkshist_cellZsort_PFC(winst:winend,:));
    end
    % TRIM ENDS
    rip_spkshist_cellZsort_Hp = rip_spkshist_cellZsort_Hp(:,50/binsize:(pret+postt-50)/binsize);
    rip_spkshist_cellZsort_PFC = rip_spkshist_cellZsort_PFC(:,50/binsize:(pret+postt-50)/binsize);
    
end

% 2) For multi spks
if domulti==1
    multimean_Hp = mean(rip_spkshist_multi_Hp,2); multistd_Hp = std(rip_spkshist_multi_Hp,[],2);
    multimean_PFC = mean(rip_spkshist_multi_PFC,2); multistd_PFC = std(rip_spkshist_multi_PFC,[],2);
    x=find(multistd_Hp==0 | multistd_PFC==0); % Remove rows which are all zero
    rip_spkshist_multi_Hp(x,:)=[]; multimean_Hp(x)=[]; multistd_Hp(x)=[]; %rip_spks_multi{x}=[];
    rip_spkshist_multi_PFC(x,:)=[]; multimean_PFC(x)=[]; multistd_PFC(x)=[]; ripsize_multi(x) = [];
    % 2a) Hp
    multimean_mat = repmat(multimean_Hp,1,size(rip_spkshist_multi_Hp,2)); multistd_mat = repmat(multistd_Hp,1,size(rip_spkshist_multi_Hp,2));
    rip_spkshist_multiZ_Hp = (rip_spkshist_multi_Hp-multimean_mat)./multistd_mat;
    % 2b) PFC
    multimean_mat = repmat(multimean_PFC,1,size(rip_spkshist_multi_PFC,2)); multistd_mat = repmat(multistd_PFC,1,size(rip_spkshist_multi_PFC,2));
    rip_spkshist_multiZ_PFC = (rip_spkshist_multi_PFC-multimean_mat)./multistd_mat;
    
    % A) SORT BY RIPPLE SIZE. Do either A or B
    % -----------------------------------------
    % [multiresp_sort,multisortidx] = sort(ripsize_multi);
    % rip_spkshist_multiZsort_Hp = rip_spkshist_multiZ_Hp(multisortidx,:);
    % rip_spkshist_multiZsort_PFC = rip_spkshist_multiZ_PFC(multisortidx,:);
    
    % % B )Sort by MU rate in -100 to 100 ms for Hp, and apply same order to PFC/sort PFC by itself
    % % ----------------------------------------------------------------------
    multiresp = sum(rip_spkshist_multiZ_Hp(:,bins_resp),2);
    [multiresp_sort,multisortidx] = sort(multiresp);
    rip_spkshist_multiZsort_Hp = rip_spkshist_multiZ_Hp(multisortidx,:);
    %multiresp = sum(rip_spkshist_multiZ_PFC(:,bins_resp),2);
    %[multiresp_sort,multisortidx] = sort(multiresp);
    rip_spkshist_multiZsort_PFC = rip_spkshist_multiZ_PFC(multisortidx,:);
    
    % -------------------------------------------------
    % Smooth Matrix
    % -------------------------------------------------
    
    for i=1:size(rip_spkshist_multiZsort_Hp,1)
        winst = i-smwin/2; winend = i+smwin/2;
        if winst<1, winst=1; end
        if winend>size(rip_spkshist_multiZsort_Hp,1), winend = size(rip_spkshist_multiZsort_Hp,1); end
        rip_spkshist_multiZsort_Hp(i,:) = mean(rip_spkshist_multiZsort_Hp(winst:winend,:));
        rip_spkshist_multiZsort_PFC(i,:) = mean(rip_spkshist_multiZsort_PFC(winst:winend,:));
    end

      % TRIM ENDS
    rip_spkshist_multiZsort_Hp = rip_spkshist_multiZsort_Hp(:,50/binsize:(pret+postt-50)/binsize);
    rip_spkshist_multiZsort_PFC = rip_spkshist_multiZsort_PFC(:,50/binsize:(pret+postt-50)/binsize);
end

% Update pret and postt due to end trims
pret=500; postt=500;


% For d2 ep4 tet 17 - remove first 300

rip_spkshist_cellZsort_Hp(1:500,:)=[]; 
rip_spkshist_cellZsort_PFC(1:500,:)=[]; 
rip_spkshist_multiZsort_Hp(1:500,:)=[]; 
rip_spkshist_multiZsort_PFC(1:500,:)=[]; 




%% ------------------------------------------------
% PLOT
% -------------------------------------------------


if docellcomb==1
    
    %% Imagesc-Hp Cell
    % -------------
    
    figure; hold on;
    %redimscreen_figforppt1;
    set(gcf,'Position',[0 Screen(4)*0.55 Screen(3)*0.2 Screen(4)*0.4])
    subplot(3,1,[1 2]); hold on;
    xaxis = -pret:binsize:postt; yaxis = 1:size(rip_spkshist_cellZsort_Hp,1);
    imagesc(xaxis,yaxis,flipud(rip_spkshist_cellZsort_Hp));
    title(['Hp CellComb Firing Z aligned to ripples-' num2str(binsize) 'ms bins'],...
        'FontSize',20,'Fontweight','normal');
    set(gca,'XLim',[-pret postt]);
    set(gca,'YLim',[0 size(rip_spkshist_cellZsort_Hp,1)]);
    %xlabel('Time(ms)','FontSize',20,'Fontweight','normal');
    
    ypts = 0:1:size(rip_spkshist_cellZsort_Hp,1);
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
    % set(gca,'xtick',[-1,(pret/binsize)-2,(pret+postt)/binsize-1],'xticklabel',{num2str([-pret,0,postt]')},...
    %      'FontSize',20,'Fontweight','normal');
    % Ticks every 100ms
  %     set(gca,'xtick',[-pret/binsize:(100/binsize):postt/binsize],'xticklabel',{num2str([-pret:100:postt]')},...
  %         'FontSize',14,'Fontweight','normal');
    
    subplot(3,1,3); hold on;
    xaxis = -pret:binsize:postt;
    plot(xaxis,mean(rip_spkshist_cellZsort_Hp),'Linewidth',3);
    set(gca,'XLim',[-pret postt]);
    xlabel('Time(ms)','FontSize',20,'Fontweight','normal');
%     set(gca,'xtick',[1:(100/binsize):((pret+postt)/binsize)+1],'xticklabel',{num2str([-pret:100:postt]')},...
%         'FontSize',14,'Fontweight','normal');
    set(gca,'YLim',[min(mean(rip_spkshist_cellZsort_Hp))-0.05 max(mean(rip_spkshist_cellZsort_Hp))+0.05]);
    ypts = min(mean(rip_spkshist_cellZsort_Hp))-0.05:0.01:max(mean(rip_spkshist_cellZsort_Hp))+0.05;
    xpts = (0)*ones(size(ypts));
    % Plot Line at 0 ms - Onset of stimulation
    plot(xpts , ypts, 'k--','Linewidth',2);
    % Plot lines at 100ms and 2000ms
    xpts = (plotline1)*ones(size(ypts));
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
    xaxis = -pret:binsize:postt; yaxis = 1:size(rip_spkshist_cellZsort_PFC,1);
    imagesc(xaxis,yaxis,flipud(rip_spkshist_cellZsort_PFC)); 
    title(['PFC CellComb Firing Z aligned to ripples-' num2str(binsize) 'ms bins'],...
        'FontSize',20,'Fontweight','normal');  
    set(gca,'XLim',[-pret postt]);
    set(gca,'YLim',[0 size(rip_spkshist_cellZsort_PFC,1)]);
    %xlabel('Time(ms)','FontSize',20,'Fontweight','normal');
    
    ypts = 0:1:size(rip_spkshist_cellZsort_PFC,1);
    xpts = (0)*ones(size(ypts));
    
    % Plot Line at 0 ms - Onset of stimulation
    plot(xpts , ypts, 'k--','Linewidth',3);
    % Plot lines at 100ms and 2000ms
    xpts = (plotline1)*ones(size(ypts));
    plot(xpts , ypts, 'k--','Linewidth',2);
    if ~isempty(plotline2)
        xpts = (plotline2)*ones(size(ypts));
        plot(xpts , ypts, 'k--','Linewidth',2);
    end
    % set(gca,'xtick',[-1,(pret/binsize)-2,(pret+postt)/binsize-1],'xticklabel',{num2str([-pret,0,postt]')},...
    %      'FontSize',20,'Fontweight','normal');
    % Ticks every 100ms
%     set(gca,'xtick',[1:(100/binsize):((pret+postt)/binsize)+1],'xticklabel',{num2str([-pret:100:postt]')},...
%         'FontSize',14,'Fontweight','normal');
    
    subplot(3,1,3); hold on;
    xaxis = -pret:binsize:postt;
    plot(xaxis, mean(rip_spkshist_cellZsort_PFC),'Linewidth',3);
    set(gca,'XLim',[-pret postt]);
    xlabel('Time(ms)','FontSize',20,'Fontweight','normal');
    set(gca,'YLim',[min(mean(rip_spkshist_cellZsort_PFC))-0.05 max(mean(rip_spkshist_cellZsort_PFC))+0.05]);
    ypts = min(mean(rip_spkshist_cellZsort_PFC))-0.05:0.1:max(mean(rip_spkshist_cellZsort_PFC))+0.05;
    xpts = (0)*ones(size(ypts));
    % Plot Line at 0 ms - Onset of stimulation
    plot(xpts , ypts, 'k--','Linewidth',2);
    % Plot lines at 100ms and 2000ms
    xpts = (plotline1)*ones(size(ypts));
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
    xaxis = -pret:binsize:postt; yaxis = 1:size(rip_spkshist_multiZsort_Hp,1);
    imagesc(xaxis,yaxis,flipud(rip_spkshist_multiZsort_Hp));
    title(['Hp Multi Firing Z aligned to ripples-' num2str(binsize) 'ms bins'],...
        'FontSize',20,'Fontweight','normal');
    set(gca,'XLim',[-pret postt]);
    set(gca,'YLim',[0 size(rip_spkshist_multiZsort_Hp,1)])
    %xlabel('Time(ms)','FontSize',20,'Fontweight','normal');
    
    ypts = 0:1:size(rip_spkshist_multiZsort_Hp,1);
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
    
    subplot(3,1,3); hold on;
    xaxis = -pret:binsize:postt;
    plot(xaxis, mean(rip_spkshist_multiZsort_Hp),'Linewidth',3);
    set(gca,'XLim',[-pret postt]);
    xlabel('Time(ms)','FontSize',20,'Fontweight','normal');
    set(gca,'YLim',[min(mean(rip_spkshist_multiZsort_Hp)) max(mean(rip_spkshist_multiZsort_Hp))+0.5]);
    ypts = min(mean(rip_spkshist_multiZsort_Hp)):0.1:max(mean(rip_spkshist_multiZsort_Hp))+0.5;
    xpts = (0)*ones(size(ypts));
    % Plot Line at 0 ms - Onset of stimulation
    plot(xpts , ypts, 'k--','Linewidth',2);
    % Plot lines at 100ms and 2000ms
    xpts = (plotline1)*ones(size(ypts));
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
    xaxis = -pret:binsize:postt; yaxis = 1:size(rip_spkshist_multiZsort_PFC,1);
    imagesc(xaxis,yaxis,flipud(rip_spkshist_multiZsort_PFC));
    title(['PFC Multi Firing Z aligned to ripples-' num2str(binsize) 'ms bins'],...
        'FontSize',20,'Fontweight','normal');
    set(gca,'XLim',[-pret postt]);
    set(gca,'YLim',[0 size(rip_spkshist_multiZsort_PFC,1)])
    %xlabel('Time(ms)','FontSize',20,'Fontweight','normal');
    
    ypts = 0:1:size(rip_spkshist_multiZsort_PFC,1);
    xpts = (0)*ones(size(ypts));
    
    % Plot Line at 0 ms - Onset of stimulation
    plot(xpts , ypts, 'k--','Linewidth',3);
    % Plot lines at 100ms and 2000ms
    xpts = (plotline1)*ones(size(ypts));
    plot(xpts , ypts, 'k--','Linewidth',2);
    if ~isempty(plotline2)
        xpts = (plotline2)*ones(size(ypts));
        plot(xpts , ypts, 'k--','Linewidth',2);
    end
    
    subplot(3,1,3); hold on;
    xaxis = -pret:binsize:postt;
    plot(xaxis, mean(rip_spkshist_multiZsort_PFC),'Linewidth',3);
    set(gca,'XLim',[-pret postt]);
    xlabel('Time(ms)','FontSize',20,'Fontweight','normal');
    set(gca,'YLim',[min(mean(rip_spkshist_multiZsort_PFC)) max(mean(rip_spkshist_multiZsort_PFC))+0.5]);
    ypts = min(mean(rip_spkshist_multiZsort_PFC)):0.1:max(mean(rip_spkshist_multiZsort_PFC))+0.5;
    xpts = (0)*ones(size(ypts));
    % Plot Line at 0 ms - Onset of stimulation
    plot(xpts , ypts, 'k--','Linewidth',2);
    % Plot lines at 100ms and 2000ms
    xpts = (plotline1)*ones(size(ypts));
    plot(xpts , ypts, 'k--','Linewidth',1);
    
    if saveg1==1,
        figfile = [figdir,prefix,'_Day',num2str(day),'_Tet',num2str(tetu),'_PFCMultiMatrix'];
        print('-dpdf', figfile);
        print('-djpeg', figfile);
        saveas(gcf,figfile,'fig');
    end
    
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



