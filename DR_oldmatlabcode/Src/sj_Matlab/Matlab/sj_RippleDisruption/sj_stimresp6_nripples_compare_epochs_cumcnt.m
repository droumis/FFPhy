
function [ep] = sj_stimresp6_nripples_compare_epochs_cumcnt (animdirect,prefix,days,allepochs,tet,std,figopt1,saveg1,figopt2,saveg2, savedata)

%% PLOT CUMULATIVE COUNT OF SIZE AND LENGTH OF RIPPLES IN DIFFERENT EPOCHS ON GIVEN DAY OR ACROSS DAYS
% eg
% sj_stimresp6_nripples_compare_epochs_cumcnt('/data25/sjadhav/RippleInterruption/RE1_direct','RE1',2,[1 2 3 4 5],5,3,1,0,1,0,0);
% std=standard deviation

%%
if nargin<7,
    figopt1 = 0; % Plot individual day figures
end
if nargin<8,
    saveg1 = 0;
end
if nargin<9,
    figopt1 = 0; % Plot individual day figures
end
if nargin<10,
    saveg1 = 0;
end
if nargin<11,
    savedata = 0;
end

%% Fixed parameters
directoryname = animdirect;
clr = {'b','g','c','m','y','k','r'};

set(0,'defaultaxesfontsize',20);set(0,'defaultaxesfontweight','normal');
set(0,'defaultaxeslinewidth',2);

%% Loop over days
allday_n_pre = []; allday_n_run = []; allday_n_post = []; allday_n_run2=[]; allday_n_post2=[];
allday_s_pre = []; allday_s_run = []; allday_s_post = []; allday_s_run2=[]; allday_s_post2=[];
allday_t_pre = []; allday_t_post = []; allday_t_run = []; allday_t_run2=[]; allday_t_post2=[];

for d = 1:length(days)
    
    day = days(d);
    n_pre = []; n_run = []; n_post = []; n_run2=[]; n_post2=[]; dio_pre = [];
    s_pre = []; s_run = []; s_post = []; s_run2=[]; s_post2=[];
    t_pre = []; t_post = []; t_run=[]; t_run2=[]; t_post2=[];
    
    ripfile = sprintf('%s/%sripples%02dstd%02d.mat', directoryname, prefix, day, std);
    load(ripfile);
    
    %%Load dio file
    DIOfile = sprintf('%s/%sDIO%02d.mat', directoryname, prefix, day);
    load(DIOfile);
    
    for ep = 1:length(allepochs)
        
        ep;
        
        epoch = allepochs(ep);
        nrip = length(ripples{day}{epoch}{tet}.startind);
        all_trip = ripples{day}{epoch}{tet}.endtime - ripples{day}{epoch}{tet}.starttime;
        srip = ripples{day}{epoch}{tet}.maxthresh;
        
        stim = DIO{day}{epoch}{15};
        ndio = length(stim.pulselength);
        
        if ep==1
            n_pre = [n_pre; nrip]; dio_pre=[dio_pre; ndio];
            t_pre = [t_pre; 1000*all_trip];
            s_pre = [s_pre; srip];
        end
        
        if ep==2
            n_run = [n_run; nrip];
            s_run = [s_run; srip];
            t_run = [t_run; 1000*all_trip];
        end
        
        if ep==3
            n_post = [n_post; nrip];
            t_post = [t_post; 1000*all_trip];
            s_post = [s_post; srip];
        end
        
        if ep==4
            n_run2 = [n_run2; nrip];
            s_run2 = [s_run2; srip];
            t_run2 = [t_run2; 1000*all_trip];
        end
        
        if ep==5
            n_post2 = [n_post2; nrip];
            t_post2 = [t_post2; 1000*all_trip];
            s_post2 = [s_post2; srip];
        end
        
        
    end
    
    %% Save across Days
    day_s_pre{d} = s_pre; day_t_pre{d} = t_pre;
    day_s_run{d} = s_run; day_t_run{d} = t_run;
    day_s_post{d} = s_post; day_t_post{d} = t_post;
    
    allday_s_pre = [allday_s_pre; s_pre]; allday_s_run = [allday_s_run; s_run]; allday_s_post = [allday_s_post; s_post];
    allday_t_pre = [allday_t_pre; t_pre]; allday_t_post = [allday_t_post; t_post]; allday_t_run = [allday_t_run; t_run];
    
    if length(allepochs)>3
        day_s_run2{d} = s_run2; day_t_run2{d} = t_run2;
        day_s_post2{d} = s_post2; day_t_post2{d} = t_post2;
        allday_s_run2 = [allday_s_run2; s_run2]; allday_s_post2 = [allday_s_post2; s_post2];
        allday_t_run2 = [allday_t_run2; t_run2]; allday_t_post2 = [allday_t_post2; t_post2];
    end
    
    %% Plot for current day
    if figopt1==1
        figure; hold on;
        redimscreen_figforppt1;
        orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
        
        hpre = histc(s_pre,[3:0.25:20]); chpre = cumsum(hpre); nchpre = chpre./max(chpre);
        hrun = histc(s_run,[3:0.25:20]); chrun = cumsum(hrun);  nchrun = chrun./max(chrun);
        hpost = histc(s_post,[3:0.25:20]); chpost = cumsum(hpost);  nchpost = chpost./max(chpost);
        
        plot([0,3, 3:0.25:20], [0;0; nchpre], ['k.-'], 'Linewidth',4, 'Markersize',6);
        plot([0,3, 3:0.25:20], [0;0; nchrun], ['r.-'], 'Linewidth',4, 'Markersize',6);
        plot([0,3, 3:0.25:20], [0;0; nchpost], ['c.-'], 'Linewidth',4, 'Markersize',6);
        if length(allepochs)>3
            hrun2 = histc(s_run2,[3:0.25:20]); chrun2 = cumsum(hrun2);  nchrun2 = chrun2./max(chrun2);
            hpost2 = histc(s_post2,[3:0.25:20]); chpost2 = cumsum(hpost2);  nchpost2 = chpost2./max(chpost2);
            plot([0,3, 3:0.25:20], [0;0; nchrun2], ['m.-'], 'Linewidth',4, 'Markersize',6);
            plot([0,3, 3:0.25:20], [0;0; nchpost2], ['b.-'], 'Linewidth',4, 'Markersize',6);
        end
        
        title(['Ripple Size in Diff Epochs (sd:' num2str(std) ') on Tet ' num2str(tet) '- Day' num2str(day)],...
            'FontSize',24,'Fontweight','normal');
        axis([0 21 -0.05 1.05]);
        xlabel('Ripple size (std)','FontSize',20,'Fontweight','normal');
        ylabel('Cumulative Proportion','FontSize',24,'Fontweight','normal');
        
        legend('Sleep1','Run1','Sleep2','Run2','Sleep3','Location','SouthEast')
        
        %text( 4, 2450,['DetRate(4,6,7) =' num2str(round(DetRateper*100)/100) '%'],'FontSize', 24, 'FontWeight','bold');
        
        if saveg1==1,
            orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
            saveas(gcf,[prefix 'RippleCumProp_tet' num2str(tet) '_Day' num2str(day) '_std' num2str(std)],'fig');
            saveas(gcf,[prefix 'RippleCumProp_tet' num2str(tet) '_Day' num2str(day) '_std' num2str(std)],'jpg');
        end
        
        
        % %%%%%% Plot Durn %%%%%%%%
        
        figure; hold on;
        redimscreen_figforppt1;
        orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
        
        hpre = histc(t_pre,[20:5:250]); chpre = cumsum(hpre); nchpre = chpre./max(chpre);
        hrun = histc(t_run,[20:5:250]); chrun = cumsum(hrun);  nchrun = chrun./max(chrun);
        hpost = histc(t_post,[20:5:250]); chpost = cumsum(hpost);  nchpost = chpost./max(chpost);
        
        plot([0, 20:5:250], [0;nchpre], ['k.-'], 'Linewidth',4, 'Markersize',6);
        plot([0, 20:5:250], [0;nchrun], ['r.-'], 'Linewidth',4, 'Markersize',6);
        plot([0, 20:5:250], [0;nchpost], ['c.-'], 'Linewidth',4, 'Markersize',6);
        if length(allepochs)>3
            hrun2 = histc(t_run2,[20:5:250]); chrun2 = cumsum(hrun2);  nchrun2 = chrun2./max(chrun2);
            hpost2 = histc(t_post2,[20:5:250]); chpost2 = cumsum(hpost2);  nchpost2 = chpost2./max(chpost2);
            plot([0, 20:5:250], [0;nchrun2], ['m.-'], 'Linewidth',4, 'Markersize',6);
            plot([0, 20:5:250], [0;nchpost2], ['b.-'], 'Linewidth',4, 'Markersize',6);
        end
        
        title(['Ripple Duration in Diff Epochs (sd:' num2str(std) ') on Tet ' num2str(tet) '- Day' num2str(day)],...
            'FontSize',24,'Fontweight','normal');
        axis([0 260 -0.05 1.05]);
        xlabel('Ripple Duration (ms)','FontSize',20,'Fontweight','normal');
        ylabel('Cumulative Proportion', 'FontSize',24,'Fontweight','normal');
        
        legend('Sleep1','Run1','Sleep2','Run2','Sleep3','Location','SouthEast')
        
        %text( 4, 2450,['DetRate(4,6,7) =' num2str(round(DetRateper*100)/100) '%'],'FontSize', 24, 'FontWeight','bold');
        
        if saveg1==1,
            orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
            saveas(gcf,[prefix 'RippleDurnCumProp_tet' num2str(tet) '_Day' num2str(day) '_std' num2str(std)],'fig');
            saveas(gcf,[prefix 'RippleDurnCumProp_tet' num2str(tet) '_Day' num2str(day) '_std' num2str(std)],'jpg');
        end
        
    end
    
    
end


%% Plot across days

if figopt2==1
    
    figure; hold on;
    redimscreen_figforppt1;
    orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
    
    hpre = histc(allday_s_pre,[3:0.25:20]); chpre = cumsum(hpre); nchpre = chpre./max(chpre);
    hrun = histc(allday_s_run,[3:0.25:20]); chrun = cumsum(hrun);  nchrun = chrun./max(chrun);
    hpost = histc(allday_s_post,[3:0.25:20]); chpost = cumsum(hpost);  nchpost = chpost./max(chpost);
    
    plot([0,3, 3:0.25:20], [0;0; nchpre], ['k.-'], 'Linewidth',4, 'Markersize',6);
    plot([0,3, 3:0.25:20], [0;0; nchrun], ['r.-'], 'Linewidth',4, 'Markersize',6);
    plot([0,3, 3:0.25:20], [0;0; nchpost], ['c.-'], 'Linewidth',4, 'Markersize',6);
    if length(allepochs)>3
        hrun2 = histc(allday_s_run2,[3:0.25:20]); chrun2 = cumsum(hrun2);  nchrun2 = chrun2./max(chrun2);
        hpost2 = histc(allday_s_post2,[3:0.25:20]); chpost2 = cumsum(hpost2);  nchpost2 = chpost2./max(chpost2);
        plot([0,3, 3:0.25:20], [0;0; nchrun2], ['m.-'], 'Linewidth',4, 'Markersize',6);
        plot([0,3, 3:0.25:20], [0;0; nchpost2], ['b.-'], 'Linewidth',4, 'Markersize',6);
    end
    
    title(['Ripple Size in Diff Epochs (sd:' num2str(std) ') on Tet ' num2str(tet) ' - Days ' num2str(min(days)) ' to ' max(num2str(days))],...
        'FontSize',24,'Fontweight','normal');
    axis([0 21 -0.05 1.05]);
    xlabel('Ripple size (std)','FontSize',20,'Fontweight','normal');
    ylabel('Cumulative Proportion','FontSize',24,'Fontweight','normal');
    
    legend('Sleep1','Run1','Sleep2','Run2','Sleep3','Location','SouthEast')
    
    %text( 4, 2450,['DetRate(4,6,7) =' num2str(round(DetRateper*100)/100) '%'],'FontSize', 24, 'FontWeight','bold');
    
    if saveg2==1,
        orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
        saveas(gcf,[prefix '_RippleCumProp_tet' num2str(tet) '_Days' min(num2str(days)) 'to' max(num2str(days)) '_std' num2str(std)],'fig');
        saveas(gcf,[prefix '_RippleCumProp_tet' num2str(tet) '_Days' min(num2str(days)) 'to' max(num2str(days)) '_std' num2str(std)],'jpg');
    end
    
    
    %
    %
    % %%%%%% Plot Durn %%%%%%%%
    
    figure; hold on;
    redimscreen_figforppt1;
    orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
    
    hpre = histc(allday_t_pre,[20:5:250]); chpre = cumsum(hpre); nchpre = chpre./max(chpre);
    hrun = histc(allday_t_run,[20:5:250]); chrun = cumsum(hrun);  nchrun = chrun./max(chrun);
    hpost = histc(allday_t_post,[20:5:250]); chpost = cumsum(hpost);  nchpost = chpost./max(chpost);
    
    
    plot([0, 20:5:250], [0;nchpre], ['k.-'], 'Linewidth',4, 'Markersize',6);
    plot([0, 20:5:250], [0;nchrun], ['r.-'], 'Linewidth',4, 'Markersize',6);
    plot([0, 20:5:250], [0;nchpost], ['c.-'], 'Linewidth',4, 'Markersize',6);
    
    if length(allepochs)>3
        hrun2 = histc(allday_t_run2,[20:5:250]); chrun2 = cumsum(hrun2);  nchrun2 = chrun2./max(chrun2);
        hpost2 = histc(allday_t_post2,[20:5:250]); chpost2 = cumsum(hpost2);  nchpost2 = chpost2./max(chpost2);
        plot([0, 20:5:250], [0;nchrun2], ['m.-'], 'Linewidth',4, 'Markersize',6);
        plot([0, 20:5:250], [0;nchpost2], ['b.-'], 'Linewidth',4, 'Markersize',6);
    end
    
    title(['Ripple Duration in Diff Epochs (sd:' num2str(std) ') on Tet ' num2str(tet) ' - Days ' num2str(min(days)) ' to ' max(num2str(days))],...
        'FontSize',24,'Fontweight','normal');
    axis([0 260 0 1.1]);
    xlabel('Ripple Duration (ms)','FontSize',20,'Fontweight','normal');
    ylabel('Cumulative Proportion','FontSize',24,'Fontweight','normal');
    
    legend('Sleep1','Run1','Sleep2','Run2','Sleep3','Location','SouthEast')
    
    %text( 4, 2450,['DetRate(4,6,7) =' num2str(round(DetRateper*100)/100) '%'],'FontSize', 24, 'FontWeight','bold');
    
    if saveg2==1,
        orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
        saveas(gcf,[prefix '_RippleDurnCumProp_tet' num2str(tet) '_Days' num2str(min(days)) 'to' max(num2str(days)) '_std' num2str(std)],'fig');
        saveas(gcf,[prefix '_RippleDurnCumProp_tet' num2str(tet) '_Days' num2str(min(days)) 'to' max(num2str(days)) '_std' num2str(std)],'jpg');
    end
    
end


if savedata==1
    savefile = sprintf('%s/ProcessedData/%sRippleStats_Days%01dto%01d_Tet%02d.mat', animdirect, prefix, min(days), max(days), tet);
    save(savefile);
end

% %%%%%%%%%%%%%%% %%%%%%%%


