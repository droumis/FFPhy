function [basefr0] = sj_multiunit_align(prefix, day, epoch, binsize, manytet, maintet, tets, saveg1)

% Shantanu - May 2012
% From sj_rippower_multiunit_calib.
% Only aligning multiunit activity to DIO times. 
% Enter times manually. Implement loading times file later

% From sj_rippower_stimampl3, which split into sj_rippower_multiunit_run and sj_rippower_multiunit_calib
% Only for calibration now, unlike sj_rippower_stimampl3  - needs manual input of time ranges

% Plot Ripple Power Calibration Curve - Change from version 1:
% Be careful what you input! If Stimulator is on, get rid of artifact (ripple_nostim files)
% Options from Multiunit Data also in here (sj_multiunitdr_rip_stimampl1.m)
% eg

% sj_rippower_multiunit_calib('RCb', 15,1,5,1,12,[3 4 9 10 11 12],0); % days 15,16
% sj_rippower_multiunit_calib('REd', 18,1,5,1,10,[3 4 5 6 10 11],0); % days 17,18
% sj_rippower_stimampl3('REd', 18,1,5,1,10,[3 4 5 6 10 11],0);
% sj_rippower_stimampl3('REd', 18,1,5,0,10,[],1);


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
    binsize=5;
end
if nargin<5,
    manytet=0; % Set to 1 if you want combine MU firing across tetrodes, 0 otherwise
end
if nargin<6,
    maintet=4; % If you want to look at MU firing on single tetrode
end
if nargin<7,
    tets=[4 5]; % If you want to look at MU firing on multiple tetrodes combined
end
if nargin<8
    saveg1=0;
end

% --------------- Parameters ---------------
pret=2000; postt=6000; %% Times to plot
plot_ex=0; % For plotting example single-trial MU rate plots
dorip = 0; % SET TO 1 IF YOU ALSO WANT TO ALIGN TO RIPPLES
sd=4; %% SD for ripples

twoplot_idx=1;
if isempty(binsize)
    binsize = 5;  %% ms, for MU Fir Rate
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
        directoryname = '/data25/sjadhav/HPExpt/JW7_direct/';
end


clr = {'b','g','c','m','y','k','r'};

%% --------------------------------------------------
% Align MU Firing Rate to Stimulations and Ripples
% -------------------------------------------------

% Set Counters
cnt=0; % for MU spikes,
cntrip=0; 

% Load dio file
%------------------
DIOfile = sprintf('%s/%sDIO%02d.mat', directoryname, prefix, day);
load(DIOfile);
stim = DIO{day}{epoch}{16};
stim_starttime = stim.pulsetimes(:,1)./10; %ms
stim_endtime = stim.pulsetimes(:,2)./10; %ms
stim_length = stim.pulselength;
stim_isi = stim.timesincelast(2:end)./10; %ms


% Load MU file
% ------------
multifile = sprintf('%s/%smulti%02d.mat', directoryname, prefix, day);
load(multifile);

% Get Multiunit timestamps (divide by 10 for ms)
cmd=sprintf('multi%d = multi{day}{epoch}{%d}/10;',maintet,maintet); eval(cmd);
for t=1:length(tets),
    currtet=tets(t);
    cmd=sprintf('multi%d = multi{day}{epoch}{%d}/10;',currtet,currtet); eval(cmd);
end


%% Set which multiunit firing rate to use
% -------------------------------------------------------------
if manytet==0
    cmd=sprintf('multi_this = multi%d;',maintet); eval(cmd); % eg %multiu=[multi4];
    multiu=[multi_this];
    tetu=maintet;
else
    multiu=[];
    for t=1:length(tets),
        currtet=tets(t);
        cmd=sprintf('curr_multi = multi%d;',currtet); eval(cmd);
        multiu=[multiu;curr_multi];
        % eg. multiu=[multi3; multi4; multi5; multi10; multi11; multi12];
    end
    tetu=tets;
end

% Align spikes to stimulation
% -----------------------------
for i =1:length(stim_starttime)
    i;
    cnt=cnt+1;
    currstim = stim_starttime(i);
    currspks =  multiu(find( (multiu>=(currstim-pret)) & (multiu<=(currstim+postt)) ));
    currspks = currspks-(currstim-pret);
    histspks = histc(currspks,[0:binsize:pret+postt]);
    stim_spks{cnt}=currspks;
    stim_spkshist(cnt,:) = histspks;
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


% Set ranges for each stimulation amplitude Or Other Parameter
% ------------------------------------------------------------
switch prefix      % BY ANIMAL
    
    case 'JW7'        
        
        switch day  
            case 2              
                nranges=1;
                [range1] = timetrans({'00:01:24', '00:07:11'},10000,2); amp1=80; %Laser power 80%, 5sec every 20 sec
                
            case 1   
                nranges=3;
                [range1] = timetrans({'00:03:56', '00:06:14'},10000,2); amp1=50; %Laser power 50%, 2 sec every 10 sec  
                [range2] = timetrans({'00:13:30', '00:19:36'},10000,2); amp2=80; %Laser power 80%, 5 sec every 20 sec  
                [range3] = timetrans({'00:20:41', '00:23:47'},10000,2); amp3=80; %Laser power 80%, 2 sec every 10 sec  
        end % end switch day
        
end % end switch prefix



%% Divide by Stimulation amplitude - Arrange in matrix form and bar graph form
% ----------------------------------------------------------------------------
% Each stimulus start time is converted to 0.1ms resolution to get index right
stim_starttime_comp = stim_starttime*10; % Convert from ms to 10000 pts in sec

Stimhist_matr=[]; Stimhist_matr_norm=[];
for n=1:nranges
    
    cmd=sprintf('startu=range%d(1);',n); eval(cmd);
    cmd=sprintf('endu=range%d(2);',n); eval(cmd);
    % Index range by comparison of stim times to range start and end
    currrange = find(stim_starttime_comp>=startu & stim_starttime_comp<=endu);
    
    temp = stim_spkshist( currrange,:);
    
%     if (nranges==1) || ( nranges>1 && n~=1) % To avoid removing stim artifact for 0 amplitude. Skip for optical stim
%         
%         if binsize==10
%             temp(:,(pret/binsize)) = 0;
%             temp(:,(pret/binsize)+1) = 0;
%             %temp(:,(pret/binsize)+2) = 0;        
%         end 
%         if binsize==5  
%             temp(:,(pret/binsize)) = 0;
%             temp(:,(pret/binsize)+1) = 0;         
%             temp(:,(pret/binsize)+2) = 0;
%         end % end binsize=5
%     end % end nranges
%     
%     % align
%     temp(:,(pret/binsize)-1)=temp(:,(pret/binsize));
%     temp(:,2:(pret/binsize))=temp(:,1:(pret/binsize)-1);
%     temp(:,2) = temp(:,1);    
     
    cmd=sprintf('stimhistall%d = temp;',n); eval(cmd);
    cmd=sprintf('stimhist%d = mean(temp,1);',n); eval(cmd);
    cmd=sprintf('stimhistnorm%d = mean(temp,1)./max(mean(temp,1));',n); eval(cmd);
    
    Stimhist_matr = [Stimhist_matr; mean(temp,1)];
    Stimhist_matr_norm = [Stimhist_matr_norm; mean(temp,1)./max(mean(temp,1))];
    
end % end nranges


%% ------------------------------------------------
% PLOT
% -------------------------------------------------


%% Matrix
% -------

Stimhist_matr=Stimhist_matr(:,1:end-1);
Stimhist_matr_norm=Stimhist_matr_norm(:,1:end-1);
figure; hold on;
redimscreen_figforppt1;
imagesc(Stimhist_matr);
title(['Multiunit Firing aligned to stimulation-' num2str(binsize) 'ms bins'],...
    'FontSize',24,'Fontweight','normal');
ylabel('Stim Amp / Laser Power','FontSize',24,'Fontweight','normal');
xlabel('Time(ms)','FontSize',20,'Fontweight','normal');
axis([0 ((pret+postt)/binsize)+2 0.5 nranges+0.5]);
ypts = 0:1:nranges+1;
xpts = ((pret/binsize)-1)*ones(size(ypts));
% Plot Line at 0 ms - Onset of stimulation
plot(xpts , ypts, 'k--','Linewidth',2);
% Plot lines at 5000ms and 2000ms
xpts = ((pret+5000-binsize)/binsize)*ones(size(ypts));
plot(xpts , ypts, 'k:','Linewidth',3);
xpts = ((pret+2000-binsize)/binsize)*ones(size(ypts));
plot(xpts , ypts, 'k:','Linewidth',3);

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
    if length(tetu)==1
        figfile = [figdir,prefix,'_Day',num2str(day),'_Tet',num2str(tetu),'_CalibMatrix'];
    else
        figfile = [figdir,prefix,'_Day',num2str(day),'_AllTets_CalibMatrix_UnNor'];
    end
    print('-dpdf', figfile);
    print('-djpeg', figfile);
    saveas(gcf,figfile,'fig');
end



%% Plot Bar Graphs for all stimulations on one plot
% ----------------------------------------------

figure; hold on;
%redimscreen_figforppt1;
redimscreen_halfvert(1);
for n=1:nranges
    cmd=sprintf('stimhist = stimhist%d;',n); eval(cmd);
    cmd=sprintf('stimhistnorm = stimhistnorm%d;',n); eval(cmd);
    cmd=sprintf('amp = amp%d;',n); eval(cmd);
    subplot(nranges,1,(nranges-n+1)); hold on;
    yplot = (1000/binsize_plot)*stimhist; %Multiunit fr in "binsize"ms bins
    
    % Or plot normalized rate
    %yplot = stimhistnorm;
    
    taxis = [-pret+binsize+binsize:binsize:postt+binsize+binsize];
    %taxis = taxis - pret/binsize;
    set(gca,'XLim',[-pret postt+binsize]);
    bar(taxis, yplot,'k');
    axis tight;
    % Get Axis if 0uA and get Baseline Fr0 from (pret/2)/binsize bins, eg. 100ms if pret is 200ms
    if n==1
        [ylimits] = get(gca,'YLim');
        y1=ylimits(1); y2=ylimits(2);
        basefr0 = mean(yplot(1:(pret/2)/binsize));
        basefr0 = round(10*basefr0)./10; % Because roundn does not exist in Lab matlab
        std0 = std(yplot(1:(pret/2)/binsize)); std0 = roundn(std0,-1);
        plot(taxis, basefr0*ones(size(taxis)),'r--','Linewidth',2 );
        %plot(taxis, (basefr0+3*std0)*ones(size(taxis)),'c-.','Linewidth',2 );
        text(-90,0.8*y2,['Fir Rate ' num2str(basefr0) 'Hz'], 'FontSize',10,'Fontweight','normal')
    end
    % Plot Line at 0ms
    ypts = 0:1:y2;
    xpts = 0*ones(size(ypts));
    plot(xpts , ypts, 'r-.','Linewidth',2);
    
    % Plot lines at 5000ms and 2000ms
    xpts = 5000*ones(size(ypts));
    plot(xpts , ypts, 'r:','Linewidth',3);
    xpts = 2000*ones(size(ypts));
    plot(xpts , ypts, 'r:','Linewidth',3);
    
    % Set Yaxis based on 0uA
    %set(gca,'YLim',[y1 y2]);   % Same axis - All Graphs, based on 0uA
    
    % Plot Baseline0 FR
    % If you want to plot baseline firing rate for 0uA on all plot, use this
    % plot(taxis, basefr0*ones(size(taxis)),'k--','Linewidth',2 );
    
    if n~=1,
        curr_basefr = mean(yplot(1:(pret/2)/binsize));
        curr_basefr = round(10*curr_basefr)./10; % Because roundn might not exist in Lab matlab
        cstd = std(yplot(1:(pret/2)/binsize)); cstd = roundn(cstd,-1);
        plot(taxis, curr_basefr*ones(size(taxis)),'b--','Linewidth',2 );
        %plot(taxis, (curr_basefr+3*cstd)*ones(size(taxis)),'g-.','Linewidth',2 );
        text(-90,0.8*y2,['Fir Rate ' num2str(curr_basefr) 'Hz'],'FontSize',10,'Fontweight','normal','Color','b');
        %axis off
    end
    
    %set(gca,'XLim',[-100 200])
    % Set Ticks
    set(gca,'ytick',[]);
    if n~=1,
        set(gca,'xtick',[]);
    end
    % Set Labels
    ylabel([num2str(amp) ' uA/power'],'FontSize',12,'Fontweight','normal');
    if n==1
        %ylabel('MU FR');
        xlabel('Time(ms)','FontSize',12,'Fontweight','normal');
    end
    % Set Title
    %     if n==nranges,
    %         title(['Multiunit Firing aligned to stimulation in ' num2str(binsize) 'ms bins - Tets ' num2str(tetu)],...
    %             'FontSize',12,'Fontweight','normal');
    %     end
    % Save Graph if asked to
    if saveg1==1,
        if length(tetu)==1
            figfile = [figdir,prefix,'_Day',num2str(day),'_Tet',num2str(tetu),'_CalibHist'];
        else
            figfile = [figdir,prefix,'_Day',num2str(day),'_AllTets_CalibHist_UnNor'];
        end
        print('-dpdf', figfile);
        print('-djpeg', figfile);
        saveas(gcf,figfile,'fig');
    end
    
end


%% Plot Bar graph for 0uA and 2nd chosen amplitude
% --------------------------------------------------

figure; hold on;
%redimscreen_widevert(0);
redimscreen_2horsubplots;
%redimscreen_figforppt1;

[ns]=[1,twoplot_idx]; % 0uA and 2nd amplitude
for i=1:length(ns),
    n=ns(i);
    cmd=sprintf('stimhist = stimhist%d;',n); eval(cmd);
    cmd=sprintf('stimhistnorm = stimhistnorm%d;',n); eval(cmd);
    cmd=sprintf('amp = amp%d;',n); eval(cmd);
    subplot(1,length(ns),i); hold on;
    
    yplot = (1000/binsize_plot)*stimhist; %Multiunit fr in "binsize"ms bins
    % Or plot normalized rate
    %yplot = stimhistnorm;
    
    taxis = [-pret+binsize+binsize:binsize:postt+binsize+binsize];
    bar(taxis, yplot,'k');
    set(gca,'XLim',[-pret postt+binsize]);
    if n==1
        [ylimits] = get(gca,'YLim');
        y1=ylimits(1); y2=ylimits(2);
        basefr0 = mean(yplot(1:(pret/2)/binsize));
        basefr0 = round(10*basefr0)./10; % Because roundn does not exist in Lab matlab
        std0 = std(yplot(1:(pret/2)/binsize));
        std0 = roundn(std0,-1);
        % Plot Baseline FR and STD
        plot(taxis, basefr0*ones(size(taxis)),'r--','Linewidth',2 );
        plot(taxis, (basefr0+3*std0)*ones(size(taxis)),'c-.','Linewidth',2 );
        text(-pret+binsize,0.8*y2,['Fir Rate ' num2str(basefr0) 'Hz'], 'FontSize',16,'Fontweight','normal');
    end
    % Plot Line at 0ms
    ypts = 0:1:y2;
    xpts = 0*ones(size(ypts));
    % Plot lines at 1000ms and 2000ms
    plot(xpts , ypts, 'r-.','Linewidth',2);
    % Plot lines at 5000ms and 2000ms
    xpts = 5000*ones(size(ypts));
    plot(xpts , ypts, 'r:','Linewidth',3);
    xpts = 2000*ones(size(ypts));
    plot(xpts , ypts, 'r:','Linewidth',3);
    % Set Yaxis based on 0uA
    set(gca,'YLim',[y1 y2]);   % Same axis - All Graphs, based on 0uA
    % Plot Baseline FR and STD
    %plot(taxis, basefr0*ones(size(taxis)),'k--','Linewidth',2 );
    %plot(taxis, (basefr0+3*std0)*ones(size(taxis)),'c-.','Linewidth',2 );
    
    if n~=1,
        curr_basefr = mean(yplot(1:(pret/2)/binsize));
        curr_basefr = round(10*curr_basefr)./10; % Because roundn does not exist in Lab matlab
        cstd=std(yplot(1:(pret/2)/binsize));
        cstd=roundn(cstd,-1);
        plot(taxis, curr_basefr*ones(size(taxis)),'b--','Linewidth',2 );
        %plot(taxis, (curr_basefr+3*cstd)*ones(size(taxis)),'g-.','Linewidth',2 );
        text(-pret+binsize,0.8*y2,['Fir Rate ' num2str(curr_basefr) 'Hz'],...
            'FontSize',16,'Fontweight','normal','Color','b');
    end
    
    
    %set(gca,'xtick',[]);
    %set(gca,'XLim',[-100 200])
    ylabel([ num2str(amp) 'uA'],'FontSize',20,'Fontweight','normal');
    if i==length(ns)
        xlabel('Time(ms)','FontSize',20,'Fontweight','normal');
    end
    if i==1
        title(['Stimulation Calibration- ' num2str(binsize) 'ms bins; Tet ' num2str(tetu)],'FontSize',24,'Fontweight','normal');
    end
end











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



