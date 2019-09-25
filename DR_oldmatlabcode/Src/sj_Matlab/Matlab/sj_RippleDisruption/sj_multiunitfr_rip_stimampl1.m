function [basefr0] = sj_multiunitfr_rip_stimampl1(prefix, day, binsize, manytet, maintet, tets, saveg1)
%% Plot Stimulation Calibration Curve from Data Taken Especially for this purpose
% eg
% sj_multiunitfr_rip_stimampl1('REb', 27,5,1,4,[3 4 8 10 11 12]);
% sj_multiunitfr_rip_stimampl1('RCa', 21,5,0,4);
% sj_multiunitfr_rip_stimampl1('REb', 20,5,0,4);

if nargin<1,
    keyboard
    error('Please enter Expt Prefix and Day No!');
end
if nargin<2,
    keyboard
    error('Please enter Day No!');
end
if nargin<3,
    binsize=5;
end
if nargin<4,
    manytet=0; % Set to 1 if you want combine MU firing across tetrodes, 0 otherwise
end
if nargin<5,
    maintet=4; % If you want to look at MU firing on single tetrode
end
if nargin<6,
    tets=[4 5]; % If you want to look at MU firing on multiple tetrodes combined
end
if nargin<7
    saveg1=0;
end

twoplot_idx=1; % for the plot comparing 0uA with 2nd stimulation amplitude, choose which index

plot_ex=0; % For plotting example single-trial MU rate plots

epoch = [4]; %% Epoch - Usually 1 for calibration

set(0,'defaultaxesfontsize',14);set(0,'defaultaxesfontweight','bold'); set(0,'defaultaxeslinewidth',2);



%% -----------------------------------------
% SET DATA
% -------------------------------------------

switch prefix
    case 'REb'
        directoryname = '/data25/sjadhav/RippleInterruption/REb_direct/StimAmpl';
    case 'RCa'
        directoryname = '/data25/sjadhav/RippleInterruption/RCa_direct/StimAmpl';
    case 'RE1'
        directoryname = '/data25/sjadhav/RippleInterruption/RE1_direct/StimAmpl';
    case 'REc'
        directoryname = '/data25/sjadhav/RippleInterruption/REc_direct/StimAmpl';
         if day>10
            directoryname = '/data25/sjadhav/RippleInterruption/REc_direct/StimAmpl';
        else
            directoryname = '/data25/sjadhav/RippleInterruption/REc_direct';
        end
end

dorip = 0; % SET TO 1 IF YOU ALSO WANT TO ALIGN TO RIPPLES

sd=4; %% SD for ripples

pret=200; postt=250; %% Times to plot
%binsize = 5;  %% ms, for MU Fir Rate
binsize_plot=binsize; %% If you put binsize_plot=1000, then units are Nspikes/binsize, not inst. firing rate in Hz

%% For LFP
eeg_pre = []; eegnostim_pre = []; rip_pre=[]; ripnostim_pre=[];
eeg_post = []; eegnostim_post = []; rip_post=[]; ripnostim_post=[];
eeg_run = []; eegnostim_run = []; rip_run=[]; ripnostim_run=[];
dio_pre = [];
s_pre = []; s_run = []; s_post = [];
t_pre = []; t_post = []; t_run=[];
clr = {'b','g','c','m','y','k','r'};

%% --------------------------------------------------
% Align MU Firing Rate to Stimulations and Ripples
% -------------------------------------------------

% Set Counters
cnt=0; cntrip=0;

% Load extracted ripple file
if dorip ==1
    ripfile = sprintf('%s/%sripples%02dstd%02d.mat', directoryname, prefix, day, sd);
    load(ripfile);
    rip_starttime = 1000* ripples{day}{epoch}{tet}.starttime;   % in msec
end

% Load dio file
DIOfile = sprintf('%s/%sDIO%02d.mat', directoryname, prefix, day);
load(DIOfile);
stim = DIO{day}{epoch}{15};
stim_starttime = stim.pulsetimes(:,1)./10; %ms
stim_endtime = stim.pulsetimes(:,2)./10; %ms
stim_length = stim.pulselength;
stim_isi = stim.timesincelast(2:end)./10; %ms

% Load MU file
multifile = sprintf('%s/%smulti%02d.mat', directoryname, prefix, day);
load(multifile);

% Get Multiunit timestamps (divide by 10 for ms)
cmd=sprintf('multi%d = multi{day}{epoch}{%d}/10;',maintet,maintet); eval(cmd);
for t=1:length(tets),
    currtet=tets(t);
    cmd=sprintf('multi%d = multi{day}{epoch}{%d}/10;',currtet,currtet); eval(cmd);
end
%multi3 = multi{day}{epoch}{3}/10 ; multi4 = multi{day}{epoch}{4}/10 ;

%% Set which multiunit firing rate to use
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

% Align spike to ripples if asked to
if dorip==1
    for i =5:length(rip_starttime)-10
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


%% Set ranges for each stimulation amplitude Or Other Parameter

switch prefix      % BY ANIMAL
    
    
    
    case 'REc'
        
        switch day        % BY DAY WITHIN ANIMAL
            
            case 15
                % 15ripdistest1 - Stim Off; Rip Detectn on 5 Tets
                % 7,8,9,10,4; SD=4; Nthr=1 and Nthr =2;
                nranges=2;
                [range1] = timetrans({'00:21:52', '00:24:52'},10000,2); amp1=1;
                [range2] = timetrans({'00:25:00', '00:28:00'},10000,2); amp2=2;
                
            case 16
                % 16ripdistest3 - Stim Off; Rip Detecn on Tet 10; Nthr=1;
                nranges=1;
                [range1] = timetrans({'00:05:22', '00:08:22'},10000,2); amp1=1; %0
                
                    case 18
                % 18lintest1 - Stim Off; Rip Detecn on 5 Tets 7,8,9,10,4;
                % Nthr=2; Lots of cable noise. 5 epochs, use only 1st one
                nranges=1;
                [range1] = timetrans({'00:47:11', '01:02:17'},10000,2); amp1=1; %0
                
            case 19
                % 19ripdistest4 - 900uA Stim ; Rip Detecn on 5 Tets 6,7,8,9,10;
                % Nthr=3; Cable noise, but rat mostly sleeping. 
                nranges=1;
                [range1] = timetrans({'00:24:49', '00:27:49'},10000,2); amp1=1; %0
                
            case 20
                % 20ripdistest5 Stim ; Rip Detecn on 4 Tets 6,8,9,10;
                % Nthr=2; Cable noise, but rat mostly sleeping.
                nranges=4;
                [range1] = timetrans({'00:48:15', '00:49:26'},10000,2); amp1=0; %0
                [range2] = timetrans({'00:45:37', '00:48:37'},10000,2); amp2=300;
                [range3] = timetrans({'00:42:20', '00:45:20'},10000,2); amp3=500;
                [range4] = timetrans({'00:39:05', '00:42:05'},10000,2); amp4=700;
                
            case 01
                % 01_112310 Stim ; Rip Detecn on 4 Tets 6,8,9,10;
                % Nthr=2; Cable noise
                nranges=1;
                [range1] = timetrans({'00:24:15', '00:39:16'},10000,2); amp1=900;
                
            case 02
                % 02_112410 Stim ; Rip Detecn on 4 Tets 6,8,9,10;
                % Nthr=2; 
                nranges=1;
                %[range1] = timetrans({'00:22:30', '00:37:31'},10000,2); amp1=1500; % epoch 2
                [range1] = timetrans({'01:15:00', '01:30:01'},10000,2); amp1=2000; % epoch 4
                
        end % end switch day
        
        
    case 'REb'
        
        switch day        % BY DAY WITHIN ANIMAL
            
            case 16
                % 16ripdistest1
                nranges=7; % 0uA to 120 uA
                [range1] = timetrans({'00:20:45', '00:22:50'},10000,2); amp1=0; %0
                [range2] = timetrans({'00:39:00', '00:42:00'},10000,2); amp2=20; %20
                [range3] = timetrans({'00:35:50', '00:38:50'},10000,2); amp3=40;%40
                [range4] = timetrans({'00:32:40', '00:35:40'},10000,2); amp4=60;%60
                [range5] = timetrans({'00:29:29', '00:32:29'},10000,2); amp5=80;%80
                [range6] = timetrans({'00:22:50', '00:25:50'},10000,2); amp6=100;%100
                [range7] = timetrans({'00:25:58', '00:28:58'},10000,2); amp7=120;%120
                
            case 17
                % 17riptest
                nranges=4;
                [range1] = timetrans({'00:01:37', '00:02:40'},10000,2); amp1=0;
                [range2] = timetrans({'00:02:48', '00:04:02'},10000,2); amp2=60;
                [range3] = timetrans({'00:04:14', '00:05:14'},10000,2); amp3=80;
                [range4] = timetrans({'00:05:23', '00:06:33'},10000,2); amp4=100;
                
            case 19
                %19ripdistest2
                nranges=5; % 0uA, 80, 100, 130 150
                [range1] = timetrans({'00:50:51', '00:52:51'},10000,2); amp1=0;
                [range2] = timetrans({'00:52:51', '00:54:51'},10000,2); amp2=80;
                [range3] = timetrans({'00:44:24', '00:46:24'},10000,2); amp3=100;
                [range4] = timetrans({'00:46:40', '00:48:40'},10000,2); amp4=130;
                [range5] = timetrans({'00:48:50', '00:50:50'},10000,2); amp5=150;
                
            case 20
                % 20ripdistest3
                nranges=5; % 0uA, 150, 300, 500 700
                [range1] = timetrans({'00:11:39', '00:13:45'},10000,2); amp1=0; %0
                [range2] = timetrans({'00:09:35', '00:11:35'},10000,2); amp2=150;%150
                [range3] = timetrans({'00:07:16', '00:09:16'},10000,2); amp3=300;%300
                [range4] = timetrans({'00:02:35', '00:04:35'},10000,2); amp4=500;%500
                [range5] = timetrans({'00:04:55', '00:06:56'},10000,2); amp5=700;%700
                
            case 21
                % 21ripdistest4
                nranges=4; % 0uA, 150, 300, 500
                [range1] = timetrans({'00:08:04', '00:10:20'},10000,2); amp1=0;%0
                [range2] = timetrans({'00:01:26', '00:03:30'},10000,2); amp2=150;%150
                [range3] = timetrans({'00:03:42', '00:05:46'},10000,2); amp3=300;%300
                [range4] = timetrans({'00:05:58', '00:08:02'},10000,2); amp4=500;%500
                
        end % end switch day
        
    case 'RCa'
        
        switch day
            
            case 20
                % 20Riptesttime
                nranges=4;
                [range1] = timetrans({'00:01:50', '00:03:50'},10000,2);
                [range2] = timetrans({'00:03:52', '00:05:52'},10000,2);
                [range3] = timetrans({'00:05:55', '00:07:55'},10000,2);
                [range4] = timetrans({'00:07:57', '00:10:00'},10000,2);
                
            case 21
                % 21Ripdis
                nranges=1;
                [range1] = timetrans({'00:02:05', '00:07:10'},10000,2); amp1=250; % 300us
                
                
            case 23
                %23RCaRipamptest_051910
                nranges=5;
                [range1] = timetrans({'00:29:25', '00:31:25'},10000,2);
                [range2] = timetrans({'00:31:26', '00:33:31'},10000,2);
                [range3] = timetrans({'00:33:42', '00:35:42'},10000,2);
                [range4] = timetrans({'00:35:50', '00:38:05'},10000,2);
                [range5] = timetrans({'00:38:13', '00:40:13'},10000,2);
                
            case 26
                %26RCaRipamptestn_051910
                nranges=7;
                [range1] = timetrans({'00:07:42', '00:09:42'},10000,2);
                [range2] = timetrans({'00:09:44', '00:11:44'},10000,2);
                [range3] = timetrans({'00:11:54', '00:13:57'},10000,2);
                [range4] = timetrans({'00:15:03', '00:17:03'},10000,2);
                [range5] = timetrans({'00:17:11', '00:19:11'},10000,2);
                [range6] = timetrans({'00:19:18', '00:21:18'},10000,2);
                [range7] = timetrans({'00:21:25', '00:23:30'},10000,2);
                
            case 27
                %27RCaRipamptestn_052110
                nranges=7;
                [range1] = timetrans({'00:03:25', '00:05:25'},10000,2);
                [range2] = timetrans({'00:05:27', '00:07:27'},10000,2);
                [range3] = timetrans({'00:07:28', '00:09:28'},10000,2);
                [range4] = timetrans({'00:09:30', '00:11:35'},10000,2);
                [range6] = timetrans({'00:11:45', '00:13:45'},10000,2);
                [range5] = timetrans({'00:13:50', '00:15:50'},10000,2);
                [range7] = timetrans({'00:15:58', '00:18:09'},10000,2);
                
        end
        
    case 'RE1'
        
        switch days
            
            case 15
                % 15StimAmpltest
                nranges=8;
                %TO DO
                
        end
        
end % end switch prefix

%% Divide by Stimulation amplitude - Arrange in matrix form and bar graph form

% Each stimulus start time is converted to 0.1ms resolution to get index right
stim_starttime_comp = stim_starttime*10; % Convert from ms to 10000 pts in sec

Stimhist_matr=[];
for n=1:nranges
    
    cmd=sprintf('startu=range%d(1);',n); eval(cmd);
    cmd=sprintf('endu=range%d(2);',n); eval(cmd);
    % Index range by comparison of stim times to range start and end
    currrange = find(stim_starttime_comp>=startu & stim_starttime_comp<=endu);
    
    temp = stim_spkshist( currrange,:);
    temp_matrix= stim_spkshist( currrange,:); % Separate for matrix to optimize color plot
    
    if n~=1
        if binsize==10
            temp_matrix(:,(pret/binsize)+1) = 0;
            temp_matrix(:,(pret/binsize)+2) = 0;
        end
        if binsize==5
            temp_matrix(:,(pret/binsize)+1) = 0;
            temp_matrix(:,(pret/binsize)+2) = 0;
            temp_matrix(:,(pret/binsize)+3) = 0;
            
%             temp(:,(pret/binsize)+1) = 0;
%             temp(:,(pret/binsize)+2) = 0;
%             temp(:,(pret/binsize)+3) = 0;
            
        end
    end
    
    
    
    cmd=sprintf('stimhistall%d = temp;',n); eval(cmd);
    cmd=sprintf('stimhist%d = mean(temp,1);',n); eval(cmd);
    
    Stimhist_matr = [Stimhist_matr; mean(temp_matrix,1)];
    
end


%% ------------------------------------------------
% PLOT
% -------------------------------------------------

%% Matrix
% -------
Stimhist_matr=Stimhist_matr(:,1:end-1);
figure; hold on;
orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
imagesc(Stimhist_matr);
title(['Multiunit Firing aligned to stimulation-' num2str(binsize) 'ms bins'],...
    'FontSize',16,'Fontweight','bold');
ylabel('Stimulation Amplitude (uA)','FontSize',16,'Fontweight','bold');
xlabel('Time(ms)','FontSize',16,'Fontweight','bold');
axis([0 ((pret+postt)/binsize)+2 0.5 nranges+0.5]);
ypts = 0:1:nranges+1;
xpts = ((pret/binsize)+1)*ones(size(ypts));
% Plot Line at 0 ms - Onset of stimulation
plot(xpts , ypts, 'k--','Linewidth',2);
% Plot lines at 100ms and 200ms
xpts = ((pret+100)/binsize)*ones(size(ypts));
plot(xpts , ypts, 'k:','Linewidth',3);
xpts = ((pret+200)/binsize)*ones(size(ypts));
plot(xpts , ypts, 'k:','Linewidth',3);

set(gca,'xtick',[0,(pret/binsize)+1,(pret+postt)/binsize],'xticklabel',{num2str([-pret,0,postt]')},...
    'FontSize',14,'Fontweight','bold');

% Make Ylabel
for n=1:nranges,
    cmd=sprintf('amp = amp%d;',n); eval(cmd);
    String{n}=num2str(amp);
    
end

set(gca,'ytick',[1:nranges],'yticklabel',String,...
    'FontSize',14,'Fontweight','bold');

if saveg1==1,
    orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
    saveas(gcf,['Image_MultiFR_StimAmpl_' num2str(day) '_tet' num2str(tet) '_SEPstd' num2str(std)],'fig');
    saveas(gcf,['Image_MultiFrStimAmpl_' num2str(day) '_tet' num2str(tet) '_SEPstd' num2str(std)],'jpg');
end



%% Plot Bar Graphs for all stimulations on one plot
% ----------------------------------------------
figure; hold on;
redimscreen_halfvert(0);
orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
for n=1:nranges
    cmd=sprintf('stimhist = stimhist%d;',n); eval(cmd);
    cmd=sprintf('amp = amp%d;',n); eval(cmd);
    subplot(nranges,1,(nranges-n+1)); hold on;
    yplot = (1000/binsize_plot)*stimhist; %Multiunit fr in "binsize"ms bins
    taxis = [-pret:binsize:postt];
    % taxis = taxis - pret/binsize;
    set(gca,'XLim',[-pret-binsize postt+binsize]);
    bar(taxis, yplot,'r');
    % Get Axis if 0uA and get Baseline Fr0 from (pret/2)/binsize bins, eg. 100ms if pret is 200ms
    if n==1
        [ylimits] = get(gca,'YLim');
        y1=ylimits(1); y2=ylimits(2);
        basefr0 = mean(yplot(1:(pret/2)/binsize));
        basefr0 = round(10*basefr0)./10; % Because roundn does not exist in Lab matlab
        std0 = std(yplot(1:(pret/2)/binsize)); std0 = roundn(std0,-1);
        plot(taxis, (basefr0+3*std0)*ones(size(taxis)),'c-.','Linewidth',2 );
        text(-pret+binsize,0.8*y2,['Fir Rate ' num2str(basefr0) 'Hz'], 'FontSize',14,'Fontweight','bold')
    end
    % Plot Line at 0ms
    ypts = 0:1:y2;
    xpts = 0*ones(size(ypts));
    plot(xpts , ypts, 'k-.','Linewidth',2);
    % Plot lines at 100ms and 200ms
    xpts = 100*ones(size(ypts));
    plot(xpts , ypts, 'k:','Linewidth',3);
    xpts = 200*ones(size(ypts));
    plot(xpts , ypts, 'k:','Linewidth',3);
    % Set Yaxis based on 0uA
    set(gca,'YLim',[y1 y2]);   % Same axis - All Graphs, based on 0uA
    % Plot Baseline FR
    plot(taxis, basefr0*ones(size(taxis)),'k--','Linewidth',2 );
    if n~=1,
        curr_basefr = mean(yplot(1:(pret/2)/binsize));
        curr_basefr = round(10*curr_basefr)./10; % Because roundn does not exist in Lab matlab
        cstd = std(yplot(1:(pret/2)/binsize)); cstd = roundn(cstd,-1);
        plot(taxis, curr_basefr*ones(size(taxis)),'b--','Linewidth',2 );
        %plot(taxis, (curr_basefr+3*cstd)*ones(size(taxis)),'g-.','Linewidth',2 );
        text(-pret+binsize,0.8*y2,['Fir Rate ' num2str(curr_basefr) 'Hz'],...
            'FontSize',14,'Fontweight','bold','Color','b');
    end
    % Set Ticks
    set(gca,'ytick',[]);
    if n~=1,
        set(gca,'xtick',[]);
    end
    % Set Labels
    ylabel([num2str(amp) ' uA'],'FontSize',16,'Fontweight','bold');
    if n==1
        %ylabel('MU FR');
        xlabel('Time(ms)');
    end
    % Set Title
    if n==nranges,
        title(['Multiunit Firing aligned to stimulation in ' num2str(binsize) 'ms bins - Tets ' num2str(tetu)],...
            'FontSize',16,'Fontweight','bold');
    end
    % Save Graph if asked to
    if saveg1==1,
        orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
        saveas(gcf,['MultiFR_StimAmpl_' num2str(day) '_tet' num2str(tet) '_SEPstd' num2str(std)],'fig');
        saveas(gcf,['MultiFrStimAmpl_' num2str(day) '_tet' num2str(tet) '_SEPstd' num2str(std)],'jpg');
    end
    
end


%% Plot Bar graph for 0uA and 2nd chosen amplitude
% --------------------------------------------------
figure; hold on;
orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
redimscreen_widevert(0);

[ns]=[1,twoplot_idx]; % 0uA and 2nd amplitude
for i=1:length(ns),
    n=ns(i);
    cmd=sprintf('stimhist = stimhist%d;',n); eval(cmd);
    cmd=sprintf('amp = amp%d;',n); eval(cmd);
    subplot(1,length(ns),i); hold on;
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
        text(-pret+binsize,0.8*y2,['Fir Rate ' num2str(basefr0) 'Hz'], 'FontSize',14,'Fontweight','bold');
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
    plot(taxis, (basefr0+3*std0)*ones(size(taxis)),'c-.','Linewidth',2 );
    if n~=1,
        curr_basefr = mean(yplot(1:(pret/2)/binsize));
        curr_basefr = round(10*curr_basefr)./10; % Because roundn does not exist in Lab matlab
        cstd=std(yplot(1:(pret/2)/binsize));
        cstd=roundn(cstd,-1);
        plot(taxis, curr_basefr*ones(size(taxis)),'b--','Linewidth',2 );
        plot(taxis, (curr_basefr+3*cstd)*ones(size(taxis)),'g-.','Linewidth',2 );
        text(-pret+binsize,0.8*y2,['Fir Rate ' num2str(curr_basefr) 'Hz'],...
            'FontSize',14,'Fontweight','bold','Color','b');
    end
    
    
    %set(gca,'xtick',[]);
    ylabel([ num2str(amp) 'uA'],'FontSize',16,'Fontweight','bold');
    if i==length(ns)
        xlabel('Time(ms)','FontSize',16,'Fontweight','bold');
    end
    if i==1
        title(['Stimulation Calibration- ' num2str(binsize) 'ms bins; Tet ' num2str(tetu)],'FontSize',14,'Fontweight','bold');
    end
end


%% Plot example single-trial firing rate and ripple power aligned to stimulation
% ------------------------------------------------------------------------

if plot_ex==1
    
    figure; hold on;
    orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
    redimscreen_widevert(0);
    
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
                text(-pret+binsize,0.8*y2,['Fir Rate ' num2str(basefr0) 'Hz'], 'FontSize',10,'Fontweight','bold');
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
                    'FontSize',10,'Fontweight','bold','Color','b');
            end
            
            
            %set(gca,'xtick',[]);
            if p==1
                title([num2str(amp) 'uA; ' num2str(binsize) 'ms bins; Tet ' num2str(tetu)],'FontSize',14,'Fontweight','bold');
                %ylabel([ num2str(amp) 'uA'],'FontSize',14,'Fontweight','bold');
            end
            
            if p==5,
                xlabel('Time(ms)','FontSize',14,'Fontweight','bold');
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
        'FontSize',24,'Fontweight','bold');
    %axis([0 800 -800 600]);
    ylabel('Instantaeous Multiunit Firing Rate');
    xlabel('Time(ms)');
    
    %text( 4, 2450,['DetRate(4,6,7) =' num2str(round(DetRateper*100)/100) '%'],'FontSize', 24, 'FontWeight','bold');
    
    if saveg1==1,
        orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
        saveas(gcf,['MultiFRrip_' num2str(day) '_tet' num2str(tet) '_SEPstd' num2str(sd)],'fig');
        saveas(gcf,['MultiFrrip_' num2str(day) '_tet' num2str(tet) '_SEPstd' num2str(sd)],'jpg');
    end
    
end



