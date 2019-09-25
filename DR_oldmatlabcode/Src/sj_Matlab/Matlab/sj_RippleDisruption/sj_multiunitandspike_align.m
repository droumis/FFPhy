function [basefr0] = sj_multiunitandspike_align(prefix, day, epoch, tet, cell, mucluster, binsize, saveg1)

% Shantanu - May 2012
% From sj_multiunit_align, which is derived from sj_rippower_multiunit_calib.
% Change to getting times directly from epoch file. From given tetrode, plot mulitunit align, and spike align 
% for either given cell or all cells on tetrode.
% Also, adding ability to use MU cluster rather than the multi file.

% sj_multiunitandspike_align('JW7', 8, 1, 3, 1, [], 200, 0);
% sj_multiunitandspike_align('JW7', 7, 1, 3, 1, [], 200, 0);
% sj_multiunitandspike_align('JW7', 7, 1, 1, 1, [], 200, 0);
% 


% Only aligning multiunit activity to DIO times. 
% Used to enter epoch times manually. Implemented loading times file

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
    tet=1; % 
end
if nargin<5,
    cell=[]; % 
end
if nargin<6,
    mucluster=[]; % If no MU cluster is specified, then will lokk for Mu tag in spikes. If does not exist, will use MU file  
end
if nargin<7,
    binsize=100;
end
if nargin<8
    saveg1=0;
end

doraster_mu=0;
plotline1 = 2000; plotline2=[]; % Where to plot lines after stimulation
Screen = get(0,'Screensize');

% --------------- Parameters ---------------
pret=3200; postt=6200; %% Times to plot - For 2 sec
%pret=1600; postt=3600; %% Times to plot - For 1 sec
plot_ex=0; % For plotting example single-trial MU rate plots
dorip = 0; % SET TO 1 IF YOU ALSO WANT TO ALIGN TO RIPPLES
sd=4; %% SD for ripples

twoplot_idx=1;
if isempty(binsize)
    binsize = 100;  %% ms, for MU Fir Rate
end
binsize_plot=binsize; %% If you put binsize_plot=1000, then units are Nspikes/binsize, not inst. firing rate in Hz

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
    case 'JW7'
        rawdir = '/data25/sjadhav/HPExpt/JW7/';
        directoryname = '/data25/sjadhav/HPExpt/JW7_direct/';
end


clr = {'b','g','c','m','y','k','r'};

%% --------------------------------------------------
% Align MU and Single Cell Firing Rate to Stimulations and Ripples
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
userange = ranges(epoch+1,:);
usename = names{epoch+1}(end-15:end);
amp=70; amp1=amp;
nranges=1;


% Load dio file
%------------------
DIOfile = sprintf('%s/%sDIO%02d.mat', directoryname, prefix, day);
load(DIOfile);
stim = DIO{day}{epoch}{16};
stim_starttime = stim.pulsetimes(:,1)./10; %ms
stim_endtime = stim.pulsetimes(:,2)./10; %ms
stim_length = stim.pulselength;
stim_isi = stim.timesincelast(2:end)./10; %ms


% Load Spikes file
% ------------
spikefile = sprintf('%s/%sspikes%02d.mat', directoryname, prefix, day);
load(spikefile);

% Get Spike timestamps (divide by 10 for ms)
for i=1:length(spikes{day}{epoch}{tet})   
    cmd=sprintf('spike%d = spikes{day}{epoch}{tet}{%d}.data/10;',i,i); eval(cmd);
end

%% Set cell to use - If cell not specified, all cells combined
% -------------------------------------------------------------
if ~isempty(cell)
    cmd=sprintf('spikeu = spike%d;',cell); eval(cmd); 
 
else
    spikeu = [];
    for i=1:length(spikes{day}{epoch}{tet})   
        cmd=sprintf('spikeu = [spikeu;spike%d];',i); eval(cmd); 
    end
end


% MU - If cluster not specified, look for MU tag in spikes structure. If not found, use multi file. 
multiu=[];
if ~isempty(mucluster)
    multiu = spikes{day}{epoch}{tet}{mucluster}.data/10;
else
    for i=1:length(spikes{day}{epoch}{tet})
        cmd=sprintf('currtag = spikes{day}{epoch}{tet}{%d}.tag;',i); eval(cmd);
        if strcmp(currtag,'MU')
            multiu = spikes{day}{epoch}{tet}{i}.data/10;
            break
        end
    end      
end
% If MU is still empty, load multi file
% -------------------------------------
if isempty(multiu)
    multifile = sprintf('%s/%smulti%02d.mat', directoryname, prefix, day);
    load(multifile);
    cmd=sprintf('multiu = multi{day}{epoch}{%d}/10;',tet); eval(cmd);
end

% ALINGING TO STIMULATION
% ----------------------

% Set Counters
cnt=0; % for MU spikes,
cntrip=0; 
% Align MU spikes to stimulation
% -----------------------------
for i =1:length(stim_starttime)
    i;
    cnt=cnt+1;
    currstim = stim_starttime(i);
    currspks =  multiu(find( (multiu>=(currstim-pret)) & (multiu<=(currstim+postt)) ));
    currspks = currspks-(currstim-pret);
    histspks = histc(currspks,[0:binsize:pret+postt]);
    stim_spks_multi{cnt}=currspks;
    stim_spkshist_multi(cnt,:) = histspks;
end

% Align Cell spikes to stimulation
% -----------------------------
cnt=0;
for i =1:length(stim_starttime)
    i;
    cnt=cnt+1;
    currstim = stim_starttime(i);
    currspks =  spikeu(find( (spikeu>=(currstim-pret)) & (spikeu<=(currstim+postt)) ));
    currspks = currspks-(currstim-pret);
    histspks = histc(currspks,[0:binsize:pret+postt]);
    stim_spks_cell{cnt}=currspks;
    stim_spkshist_cell(cnt,:) = histspks;
end

% If asked to align to ripples
% ---------------------------
if dorip ==1
    % Load extracted ripple file
    ripfile = sprintf('%s/%sripplesep1%02d.mat', directoryname, prefix, day);
    load(ripfile);
    rip_starttime = 1000* ripples{day}{epoch}{tet}.starttime;   % in msec
    
    % Align spike to ripples    
    for i =2:length(rip_starttime)-1
        i;
        cntrip=cntrip+1;
        currrip = rip_starttime(i);
        currspks =  multiu(find( (multiu>=(currrip-pret)) & (multiu<=(currrip+postt)) ));
        currspks = currspks-(currrip-pret);
        histspks = histc(currspks,[0:binsize:pret+postt]);
        rip_spks{cntrip}=currspks;
        rip_spkshist(cntrip,:) = histspks;
    end
end

% Now Getting Range directly from times file
% % Set ranges for each stimulation amplitude Or Other Parameter
% % ------------------------------------------------------------
% switch prefix      % BY ANIMAL
%     
%     case 'JW7'        
%         
%         switch day  
%             case 2              
%                 nranges=1;
%                 [range1] = timetrans({'00:01:24', '00:07:11'},10000,2); amp1=80; %Laser power 80%, 5sec every 20 sec
%                 
%             case 1   
%                 nranges=3;
%                 [range1] = timetrans({'00:03:56', '00:06:14'},10000,2); amp1=50; %Laser power 50%, 2 sec every 10 sec  
%                 [range2] = timetrans({'00:13:30', '00:19:36'},10000,2); amp2=80; %Laser power 80%, 5 sec every 20 sec  
%                 [range3] = timetrans({'00:20:41', '00:23:47'},10000,2); amp3=80; %Laser power 80%, 2 sec every 10 sec  
%         end % end switch day
%         
% end % end switch prefix



%% Arrange in matrix form and bar graph form
% ------------------------------------------
Stimhist_matr_multi=mean(stim_spkshist_multi,1);
Stimhist_matr_norm_multi = [mean(stim_spkshist_multi,1)./max(mean(stim_spkshist_multi,1))];
stim_meanhist_multi = mean(stim_spkshist_multi,1);
stim_meanhistnorm_multi = stim_meanhist_multi./max(stim_meanhist_multi);

Stimhist_matr_cell=mean(stim_spkshist_cell,1);
Stimhist_matr_norm_cell = [mean(stim_spkshist_cell,1)./max(mean(stim_spkshist_cell,1))];
stim_meanhist_cell = mean(stim_spkshist_cell,1);
stim_meanhistnorm_cell = stim_meanhist_cell./max(stim_meanhist_cell);

% % Now, no division by Stimulation amplitude
% % ----------------------------------------------------------------------------
% % Each stimulus start time is converted to 0.1ms resolution to get index right
% stim_starttime_comp = stim_starttime*10; % Convert from ms to 10000 pts in sec
% Stimhist_matr=[]; Stimhist_matr_norm=[];
% for n=1:nranges
%     
%     cmd=sprintf('startu=range%d(1);',n); eval(cmd);
%     cmd=sprintf('endu=range%d(2);',n); eval(cmd);
%     % Index range by comparison of stim times to range start and end
%     currrange = find(stim_starttime_comp>=startu & stim_starttime_comp<=endu);
%     
%     temp = stim_spkshist( currrange,:);
%     
% %     if (nranges==1) || ( nranges>1 && n~=1) % To avoid removing stim artifact for 0 amplitude. Skip for optical stim
% %         
% %         if binsize==10

% %         end % end binsize
% %     end % end nranges
% %     
%      
%     cmd=sprintf('stimhistall%d = temp;',n); eval(cmd);
%     cmd=sprintf('stimhist%d = mean(temp,1);',n); eval(cmd);
%     cmd=sprintf('stimhistnorm%d = mean(temp,1)./max(mean(temp,1));',n); eval(cmd);
%     
%     Stimhist_matr = [Stimhist_matr; mean(temp,1)];
%     Stimhist_matr_norm = [Stimhist_matr_norm; mean(temp,1)./max(mean(temp,1))];
%     
% end % end nranges


%% ------------------------------------------------
% PLOT
% -------------------------------------------------


%% Matrix - MU
% -------

Stimhist_matr=Stimhist_matr_multi(:,1:end-1);
Stimhist_matr_norm=Stimhist_matr_norm_multi(:,1:end-1);
figure; hold on;
%redimscreen_figforppt1;
set(gcf,'Position',[0 Screen(4)*0.55 Screen(3)*0.2 Screen(4)*0.4])
imagesc(Stimhist_matr);
title(['Multiunit Firing aligned to stimulation-' num2str(binsize) 'ms bins'],...
    'FontSize',20,'Fontweight','normal');
ylabel('Stim Amp / Laser Power','FontSize',20,'Fontweight','normal');
xlabel('Time(ms)','FontSize',20,'Fontweight','normal');
axis([0 ((pret+postt)/binsize)+2 0.5 nranges+0.5]);
ypts = 0:1:nranges+1;
xpts = ((pret/binsize)-1)*ones(size(ypts));
% Plot Line at 0 ms - Onset of stimulation
plot(xpts , ypts, 'k--','Linewidth',4);
% Plot lines at 5000ms and 2000ms
xpts = ((pret+plotline1-binsize)/binsize)*ones(size(ypts));
plot(xpts , ypts, 'k--','Linewidth',3);
 if ~isempty(plotline2)
    xpts = ((pret+plotline2-binsize)/binsize)*ones(size(ypts));
    plot(xpts , ypts, 'k--','Linewidth',3);
 end
set(gca,'xtick',[-1,(pret/binsize)-2,(pret+postt)/binsize-1],'xticklabel',{num2str([-pret,0,postt]')},...
    'FontSize',20,'Fontweight','normal');

% Make Ylabel
for n=1:nranges,
    cmd=sprintf('amp = amp%d;',n); eval(cmd);
    String{n}=num2str(amp);    
end
set(gca,'ytick',[1:nranges],'yticklabel',String,...
    'FontSize',20,'Fontweight','normal');

%set(gca,'XLim',[19 79]); % -100 to 200 ms 

if saveg1==1,
    figfile = [figdir,prefix,'_Day',num2str(day),'_Tet',num2str(tetu),'_MultiMatrix'];
    print('-dpdf', figfile);
    print('-djpeg', figfile);
    saveas(gcf,figfile,'fig');
end


%% Plot Bar Graphs for MU
% ----------------------------------------------

figure; hold on;
%redimscreen_figforppt1;
%redimscreen_halfvert(1);
set(gcf,'Position',[0 Screen(4)*0.03 Screen(3)*0.2 Screen(4)*0.4])
for n=1:nranges
    cmd=sprintf('amp = amp%d;',n); eval(cmd);
    subplot(nranges,1,(nranges-n+1)); hold on;
    yplot = (1000/binsize_plot)*stim_meanhist_multi; %Multiunit fr in "binsize"ms bins
    
    % Or plot normalized rate
    %yplot = stimhistnorm;
    
    taxis = [-pret+binsize+binsize:binsize:postt+binsize+binsize];
    %taxis = taxis - pret/binsize;
    set(gca,'XLim',[-pret postt+binsize]);
    bar(taxis, yplot,'k');
    %axis tight;
    
    [ylimits] = get(gca,'YLim');
    y1=ylimits(1); y2=ylimits(2);
    basefr0 = mean(yplot(1:((pret/binsize)-1)));
    basefr0 = round(10*basefr0)./10; % Because roundn does not exist in Lab matlab
    std0 = std(yplot(1:(pret/2)/binsize)); std0 = roundn(std0,-1);
    plot(taxis, basefr0*ones(size(taxis)),'c--','Linewidth',2 );
    %plot(taxis, (basefr0+3*std0)*ones(size(taxis)),'c-.','Linewidth',2 );
    text(-2000,0.95*y2,['Mean fr ' num2str(basefr0) 'Hz'], 'FontSize',16,'Fontweight','normal');
    
    % Plot Line at 0ms
    ypts = 0:1:y2;
    xpts = 0*ones(size(ypts));
    plot(xpts , ypts, 'r-','Linewidth',3);
    % Plot lines at 5000ms and 2000ms
    xpts = plotline1*ones(size(ypts));
    plot(xpts , ypts, 'r--','Linewidth',3);
    if ~isempty(plotline2)
        xpts = plotline2*ones(size(ypts));
        plot(xpts , ypts, 'r--','Linewidth',3);
    end
    %set(gca,'XLim',[-100 200])
    % Set Ticks
%     set(gca,'ytick',[]);
%     if n~=1,
%         set(gca,'xtick',[]);
%     end
    % Set Labels
    title(['Multiunit Firing aligned to stimulation-' num2str(binsize) 'ms bins'],...
    'FontSize',20,'Fontweight','normal');
    ylabel('Firing Rate (Hz)','FontSize',20,'Fontweight','normal');   
    xlabel('Time(ms)','FontSize',20,'Fontweight','normal');
   
    % Save Graph if asked to
    if saveg1==1,
        figfile = [figdir,prefix,'_Day',num2str(day),'_Tet',num2str(tetu),'_Multi_Hist'];
        print('-dpdf', figfile);
        print('-djpeg', figfile);
        saveas(gcf,figfile,'fig');
    end
    
end




%% Matrix - Cell
% -------

Stimhist_matr=Stimhist_matr_cell(:,1:end-1);
Stimhist_matr_norm=Stimhist_matr_norm_cell(:,1:end-1);
figure; hold on;
%redimscreen_figforppt1;
set(gcf,'Position',[Screen(3)*0.205 Screen(4)*0.55 Screen(3)*0.2 Screen(4)*0.4]);
imagesc(Stimhist_matr);
title(['Cell ' num2str(cell) ' firing aligned to stimulation-' num2str(binsize) 'ms bins'],...
    'FontSize',20,'Fontweight','normal');
ylabel('Stim Amp / Laser Power','FontSize',20,'Fontweight','normal');
xlabel('Time(ms)','FontSize',20,'Fontweight','normal');
axis([0 ((pret+postt)/binsize)+2 0.5 nranges+0.5]);
ypts = 0:1:nranges+1;
xpts = ((pret/binsize)-1)*ones(size(ypts));
% Plot Line at 0 ms - Onset of stimulation
plot(xpts , ypts, 'k--','Linewidth',4);
% Plot lines at 5000ms and 2000ms
xpts = ((pret+plotline1-binsize)/binsize)*ones(size(ypts));
plot(xpts , ypts, 'k--','Linewidth',3);
 if ~isempty(plotline2)
    xpts = ((pret+plotline2-binsize)/binsize)*ones(size(ypts));
    plot(xpts , ypts, 'k--','Linewidth',3);
 end
set(gca,'xtick',[-1,(pret/binsize)-2,(pret+postt)/binsize-1],'xticklabel',{num2str([-pret,0,postt]')},...
    'FontSize',20,'Fontweight','normal');

% Make Ylabel
for n=1:nranges,
    cmd=sprintf('amp = amp%d;',n); eval(cmd);
    String{n}=num2str(amp);    
end
set(gca,'ytick',[1:nranges],'yticklabel',String,...
    'FontSize',20,'Fontweight','normal');

%set(gca,'XLim',[19 79]); % -100 to 200 ms 

if saveg1==1,
    figfile = [figdir,prefix,'_Day',num2str(day),'_Tet',num2str(tetu),'Cell',num2str(cell),'_Matrix'];
    print('-dpdf', figfile);
    print('-djpeg', figfile);
    saveas(gcf,figfile,'fig');
end







%% Plot Bar Graphs for Cell
% ----------------------------------------------
figure; hold on;
%redimscreen_figforppt1;
%redimscreen_halfvert(1);
set(gcf,'Position',[Screen(3)*0.205 Screen(4)*0.03 Screen(3)*0.2 Screen(4)*0.4])
for n=1:nranges
    cmd=sprintf('amp = amp%d;',n); eval(cmd);
    subplot(nranges,1,(nranges-n+1)); hold on;
    yplot = (1000/binsize_plot)*stim_meanhist_cell; %Multiunit fr in "binsize"ms bins
    
    % Or plot normalized rate
    %yplot = stimhistnorm;
    
    taxis = [-pret+binsize+binsize:binsize:postt+binsize+binsize];
    %taxis = taxis - pret/binsize;
    set(gca,'XLim',[-pret postt+binsize]);
    bar(taxis, yplot,'k');
    %axis tight;
    
    [ylimits] = get(gca,'YLim');
    y1=ylimits(1); y2=ylimits(2);
    basefr0 = mean(yplot(1:((pret/binsize)-1)));
    basefr0 = round(10*basefr0)./10; % Because roundn does not exist in Lab matlab
    std0 = std(yplot(1:(pret/2)/binsize)); std0 = roundn(std0,-1);
    plot(taxis, basefr0*ones(size(taxis)),'c--','Linewidth',2 );
    %plot(taxis, (basefr0+3*std0)*ones(size(taxis)),'c-.','Linewidth',2 );
    text(-2000,0.95*y2,['Mean fr ' num2str(basefr0) 'Hz'], 'FontSize',16,'Fontweight','normal')
    
    % Plot Line at 0ms
    ypts = 0:1:y2;
    xpts = 0*ones(size(ypts));
    plot(xpts , ypts, 'r-','Linewidth',3);
    % Plot lines at 5000ms and 2000ms
    xpts = plotline1*ones(size(ypts));
    plot(xpts , ypts, 'r--','Linewidth',3);
     if ~isempty(plotline2)
        xpts = plotline2*ones(size(ypts));
        plot(xpts , ypts, 'r--','Linewidth',3);
     end
    %set(gca,'XLim',[-100 200])
    % Set Ticks
%     set(gca,'ytick',[]);
%     if n~=1,
%         set(gca,'xtick',[]);
%     end
    % Set Labels
    title(['Cell ' num2str(cell) ' firing aligned to stimulation-' num2str(binsize) 'ms bins'],...
    'FontSize',20,'Fontweight','normal');
    ylabel('Firing Rate (Hz)','FontSize',20,'Fontweight','normal');   
    xlabel('Time(ms)','FontSize',20,'Fontweight','normal');
   
    % Save Graph if asked to
    if saveg1==1,
        figfile = [figdir,prefix,'_Day',num2str(day),'_Tet',num2str(tetu),'_Cell',num2str(cell),'_Hist'];
        print('-dpdf', figfile);
        print('-djpeg', figfile);
        saveas(gcf,figfile,'fig');
    end
    
end

  
%% Plot Raster for MU
% ----------------------------------------------
if doraster_mu==1
    figure; hold on;
    set(gcf,'Position',[Screen(3)*0.41 Screen(4)*0.55 Screen(3)*0.2 Screen(4)*0.4]);
    
    for tr=1:length(stim_spks_multi)
        
        currspkt = stim_spks_multi{tr}; % Spikes for current trial
        currspkt = currspkt-pret; % Set stimlation time to 0
        
        for i=1:length(currspkt),
            sj_plotraster(currspkt(i),tr,0.8);
        end
        
    end
    
    [ylimits] = get(gca,'YLim');
    y1=ylimits(1); y2=ylimits(2);
    % Plot Line at 0ms
    ypts = 0:1:y2;
    xpts = 0*ones(size(ypts));
    plot(xpts , ypts, 'r-','Linewidth',3);
    % Plot lines at 5000ms and 2000ms
    xpts = plotline1*ones(size(ypts));
    plot(xpts , ypts, 'r--','Linewidth',3);
     if ~isempty(plotline2)
        xpts = plotline2*ones(size(ypts));
        plot(xpts , ypts, 'r--','Linewidth',3);
     end
    % Set Labels
    title(['MultiUnit Raster'],'FontSize',20,'Fontweight','normal');
    ylabel('Trial Number','FontSize',20,'Fontweight','normal');
    xlabel('Time(ms)','FontSize',20,'Fontweight','normal');
    
    % Save Graph if asked to
    if saveg1==1,
        figfile = [figdir,prefix,'_Day',num2str(day),'_Tet',num2str(tetu),'_Cell',num2str(cell),'_Raster'];
        print('-dpdf', figfile);
        print('-djpeg', figfile);
        saveas(gcf,figfile,'fig');
    end
end



%% Plot Raster for Cell - Using tick marks: Can also use sparse matrix strategy
% ----------------------------------------------
figure; hold on;
set(gcf,'Position',[Screen(3)*0.41 Screen(4)*0.03 Screen(3)*0.2 Screen(4)*0.4]);
for tr=1:length(stim_spks_cell)
    
    currspkt = stim_spks_cell{tr}; % Spikes for current trial
    currspkt = currspkt-pret; % Set stimlation time to 0 
    
    for i=1:length(currspkt),
        sj_plotraster(currspkt(i),tr,0.8);
    end
    
end

[ylimits] = get(gca,'YLim');
y1=ylimits(1); y2=ylimits(2);
% Plot Line at 0ms
ypts = 0:1:y2;
xpts = 0*ones(size(ypts));
plot(xpts , ypts, 'r-','Linewidth',3);
% Plot lines at 5000ms and 2000ms
xpts = plotline1*ones(size(ypts));
plot(xpts , ypts, 'r--','Linewidth',3);
 if ~isempty(plotline2)
    xpts = plotline2*ones(size(ypts));
    plot(xpts , ypts, 'r--','Linewidth',3);
 end
% Set Labels
title(['Cell ' num2str(cell) ' Raster'],'FontSize',20,'Fontweight','normal');
ylabel('Trial Number','FontSize',20,'Fontweight','normal');
xlabel('Time(ms)','FontSize',20,'Fontweight','normal');

% Save Graph if asked to
if saveg1==1,
    figfile = [figdir,prefix,'_Day',num2str(day),'_Tet',num2str(tetu),'_Cell',num2str(cell),'_Raster'];
    print('-dpdf', figfile);
    print('-djpeg', figfile);
    saveas(gcf,figfile,'fig');
end







%figfile = ['Figures/JW7_d7t1c1_inh_hist100_1sec'];saveas(gcf,figfile,'fig');print('-djpeg', figfile);print('-dpdf', figfile);

































%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%% Plot example single-trial firing rate and ripple power aligned to stimulation
% ------------------------------------------------------------------------

if plot_ex==1
    
    figure; hold on;
    orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
    redimscreen_figforppt1;
    
    [ns]=[1,twoplot_idx]; % 0uA and 2nd amplitude
    for i=1:length(ns),
        n=ns(i);
        cmd=sprintf('stimhistall = stimhistall%d;',n); eval(cmd);
        cmd=sprintf('amp = amp%d;',n); eval(cmd);
        
        % Get 5 random indices from trials in stimhistall
        vec=1:size(stimhistall,1);
        r=randperm(length(vec));
        tr=vec(r(1:5));
        for p=1:5
            subplot(5,length(ns),2*(p-1)+i); hold on;
            
            stimhist=stimhistall(tr(p),:);
            yplot = (1000/binsize_plot)*stimhist; %Multiunit fr in "binsize"ms bins
            taxis = [-pret:binsize:postt];
            bar(taxis, yplot,'r');
            set(gca,'XLim',[-pret-binsize postt+binsize]);
            if n==1
                [ylimits] = get(gca,'YLim');
                y1=ylimits(1); y2=ylimits(2);
                basefr0 = mean(yplot(1:(pret/2)/binsize));
                basefr0 = round(10*basefr0)./10; % Because roundn does not exist in Lab matlab
                std0 = std(yplot(1:(pret/2)/binsize));
                std0 = roundn(std0,-1);
                plot(taxis, (basefr0+3*std0)*ones(size(taxis)),'c-.','Linewidth',2 );
                text(-pret+binsize,0.8*y2,['Fir Rate ' num2str(basefr0) 'Hz'], 'FontSize',14,'Fontweight','normal');
            end
            % Plot Line at 0ms
            ypts = 0:1:y2;
            xpts = 0*ones(size(ypts));
            % Plot lines at 100ms and 200ms
            plot(xpts , ypts, 'k-.','Linewidth',2);
            % Plot lines at 100ms and 200ms
            xpts = 100*ones(size(ypts));
            plot(xpts , ypts, 'k:','Linewidth',3);
            xpts = 200*ones(size(ypts));
            plot(xpts , ypts, 'k:','Linewidth',3);
            % Set Yaxis based on 0uA
            set(gca,'YLim',[y1 y2]);   % Same axis - All Graphs, based on 0uA
            % Plot Baseline FR and STD
            plot(taxis, basefr0*ones(size(taxis)),'k--','Linewidth',2 );
            
            if n~=1,
                curr_basefr = mean(yplot(1:(pret/2)/binsize));
                curr_basefr = round(10*curr_basefr)./10; % Because roundn does not exist in Lab matlab
                cstd=std(yplot(1:(pret/2)/binsize));
                cstd=roundn(cstd,-1);
                plot(taxis, curr_basefr*ones(size(taxis)),'b--','Linewidth',2 );
                plot(taxis, (curr_basefr+3*cstd)*ones(size(taxis)),'g-.','Linewidth',2 );
                text(-pret+binsize,0.8*y2,['Fir Rate ' num2str(curr_basefr) 'Hz'],...
                    'FontSize',14,'Fontweight','normal','Color','b');
            end
            
            
            %set(gca,'xtick',[]);
            if p==1
                title([num2str(amp) 'uA; ' num2str(binsize) 'ms bins; Tet ' num2str(tetu)],'FontSize',24,'Fontweight','normal');
                %ylabel([ num2str(amp) 'uA'],'FontSize',14,'Fontweight','bold');
            end
            
            if p==5,
                xlabel('Time(ms)','FontSize',24,'Fontweight','normal');
            end
            
            %     if i==1
            %         title(['Stimulation Calibration'],'FontSize',12,'Fontweight','bold');
            %     end
            
            
        end
        
    end
    
end


%% Plot Multiunit rate around ripples
% -------------------------------------

if dorip==1
    %plot(taxis, 2*ripnostim_pre(i,:),['r-'],'Linewidth',2,'Markersize',6);
    
    ypts = 0:1.1*max(yplot);
    xpts = 200*ones(size(ypts));
    plot(xpts , ypts, 'k--','Linewidth',2);
    
    title(['Sleep Multiunit Firing around ripples-Day' num2str(day)],...
        'FontSize',24,'Fontweight','normal');
    %axis([0 800 -800 600]);
    ylabel('Instantaeous Multiunit Firing Rate','FontSize',24,'Fontweight','normal');
    xlabel('Time(ms)','FontSize',20,'Fontweight','normal');
    
    %text( 4, 2450,['DetRate(4,6,7) =' num2str(round(DetRateper*100)/100) '%'],'FontSize', 24, 'FontWeight','bold');
    
    if saveg1==1,
        orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
        saveas(gcf,['MultiFRrip_' num2str(day) '_tet' num2str(tet) '_SEPstd' num2str(sd)],'fig');
        saveas(gcf,['MultiFrrip_' num2str(day) '_tet' num2str(tet) '_SEPstd' num2str(sd)],'jpg');
    end
    
end



