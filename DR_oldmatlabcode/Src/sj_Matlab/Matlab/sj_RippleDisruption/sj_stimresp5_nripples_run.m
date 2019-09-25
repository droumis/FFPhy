
function [ep] = sj_stimresp5_nripples_run (animdirect,prefix,days,allepochs,tet,std,saveg1)

%%%% PLOT NO of ripples and stimulations in given run epoch/s across days 
% eg
% sj_stimresp5_nripples_run('/data25/sjadhav/RippleInterruption/RE1_direct','RE1',1:7,[2 4],5,3);
%std=standard deviation

%%
if nargin<7,
    saveg1 = 0;
end

%% Getfiles for day

if (animdirect(end) == '/')
   animdirect = animdirect(1:end-1);
end
cd(animdirect);
clr = {'b','g','c','m','y','k','r'};

% directoryname = '/data25/sjadhav/SJStimC_direct';
% prefix = 'sjc';
% directoryname = '/data25/sjadhav/RE1_direct';
% prefix = 'RE1';
% days = [1:7];  allepochs = [2 4];
% tet=5;
% std = 3;

n_run1 = []; n_run2 = []; dio_run1 = [];
t_run1 = []; t_run2 = []; dio_run2 = [];

for d = 1:length(days)
    
    day = days(d);
    
    %% Load Ripple FIle
    ripfile = sprintf('%s/%sripples%02dstd%02d.mat', animdirect, prefix, day, std);
    load(ripfile);
    
     %%Load dio file
    DIOfile = sprintf('%s/%sDIO%02d.mat', animdirect, prefix, day);
    load(DIOfile);
    
    %%% Loop over epochs
    for ep = 1:length(allepochs)
        
        epoch = allepochs(ep);
        nrip = length(ripples{day}{epoch}{tet}.startind);
        all_trip = ripples{day}{epoch}{tet}.endtime - ripples{day}{epoch}{tet}.starttime;
        
        stim = DIO{day}{epoch}{15};
        ndio = length(stim.pulselength);
    
        if ep==1
            n_run1 = [n_run1; nrip]; dio_run1=[dio_run1; ndio];
            t_run1 = [t_run1; 1000*mean(all_trip)]; 
            length_ep1(d) = (ripples{day}{epoch}{tet}.timerange(2) - ripples{day}{epoch}{tet}.timerange(1));
            
        end
        if ep==2            
            n_run2 = [n_run2; nrip]; dio_run2=[dio_run2; ndio];
            t_run2 = [t_run2; 1000*mean(all_trip)];
            length_ep2(d) = (ripples{day}{epoch}{tet}.timerange(2) - ripples{day}{epoch}{tet}.timerange(1));
        end
        
    end
    
end

%npre5 = n_pre(5); n_pre(5) = []; tpre5 = t_pre(5); t_pre(5) = []; diopre5=dio_pre(5); 
%npost5 = n_post(5); n_post(5) = []; tpost5 = t_post(5); t_post(5) = [];


%%%%%% 1st Plot - Now Defunct Plot Nripples  %%%%%%%%
% figure; hold on;
% redimscreen_halfvert(0);
% orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
% plot(days,n_pre, ['ko--'], 'Linewidth',4, 'Markersize',12);
% plot(days,n_pre+dio_pre, ['rs-'], 'Linewidth',2, 'Markersize',12);
% title(['No of ripples (SEP sd:' num2str(std) ') on Tet ' num2str(tet) ' in Epoch' num2str(epoch)],...
%     'FontSize',24,'Fontweight','bold');
% %axis([0.8 2.2 0 1400]);
% xlabel('Day');
% ylabel('No of ripples (red = stimln)');
% DetRate = dio_pre([1:7])./(n_pre([1:7]) + dio_pre([1:7]));
% DetRateper = mean(DetRate*100);
% %[h,p,ci] = ttest2(n_pre, n_post, 0.05, 'left' );
% text( 4, 2450,['DetRate(Days 1:7) =' num2str(round(DetRateper*100)/100) '%'],'FontSize', 24, 'FontWeight','bold');
% if saveg1==1,
%     orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
%     saveas(gcf,['RippleN_Run_tet' num2str(tet) '_SEPstd' num2str(std)],'fig');
%     saveas(gcf,['RippleN_Run_tet' num2str(tet) '_SEPstd' num2str(std)],'jpg');
% end



%%%% Plot Ripple/Stimulation Rate for each epoch individually and combined if it exists %%%%

%%% For Each Epoch %%%%

for ep = 1:length(allepochs)
    
    epoch = allepochs(ep);
    
    figure; hold on;
    redimscreen_halfvert(0);
    orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
    
    cmd = sprintf('dio_run = dio_run%d;', ep); eval(cmd);
    cmd = sprintf('n_run = n_run%d;', ep); eval(cmd);
    cmd = sprintf('length_ep = length_ep%d;', ep); eval(cmd);
    
    stim_rate = dio_run./length_ep';
    allriprate = (dio_run + n_run) ./length_ep';
    
    plot(days,allriprate, ['co--'], 'Linewidth',4, 'Markersize',12);
    plot(days,stim_rate, ['rs-'], 'Linewidth',2, 'Markersize',12);
    
    title(['Ripple and Stimln Rate (SEP sd:' num2str(std) ') on Tet ' num2str(tet) ' in Epoch' num2str(epoch)],...
        'FontSize',24,'Fontweight','bold');
    axis([0 max(days)+1 0 4]);
    xlabel('Day');
    ylabel('Ripple / Stimuln Rate(Hz)');
    
    DetRate = dio_run([days])./(n_run([days]) + dio_run([days]));
    DetRateper = mean(DetRate*100);
    
    %[h,p,ci] = ttest2(n_run1, n_post, 0.05, 'left' );
    text( 4, 3.5,['DetRate(Days 1:' num2str(max(days)) ') =' num2str(round(DetRateper*100)/100) '%'],'FontSize', 24, 'FontWeight','bold');
    
    if saveg1==1,
        orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
        saveas(gcf,['RippleStimRate_Epoch' num2str(epoch) '_tet' num2str(tet) '_SEPstd' num2str(std)],'fig');
        saveas(gcf,['RippleStimRate_Epoch' num2str(epoch) '_tet' num2str(tet) '_SEPstd' num2str(std)],'jpg');
    end
    
end

if length(allepochs)>1,
    
    cnt_dio_run=0; cnt_n_run=0; total_length=0;
    
    for ep = 1:length(allepochs)
        epoch = allepochs(ep);
        
        cmd = sprintf('dio_run = dio_run%d;', ep); eval(cmd);
        cmd = sprintf('n_run = n_run%d;', ep); eval(cmd);
        cmd = sprintf('length_ep = length_ep%d;', ep); eval(cmd);
        
        cnt_dio_run = cnt_dio_run + dio_run;
        cnt_n_run = cnt_n_run + n_run;
        total_length = total_length + length_ep;       
    
    end
    
    figure; hold on;
    redimscreen_halfvert(0);
    orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
    
    stim_rate = cnt_dio_run./total_length';
    allriprate = (cnt_dio_run + cnt_n_run) ./total_length';
    
    plot(days,allriprate, ['co--'], 'Linewidth',4, 'Markersize',12);
    plot(days,stim_rate, ['rs-'], 'Linewidth',2, 'Markersize',12);
    
     title(['Ripple and Stimln Rate (SEP sd:' num2str(std) ') on Tet ' num2str(tet) ' in all Run Epochs'],...
        'FontSize',24,'Fontweight','bold');
    axis([0 max(days)+1 0 4]);
    xlabel('Day');
    ylabel('Ripple / Stimuln Rate(Hz)');
    
    DetRate = cnt_dio_run([days])./(cnt_n_run([days]) + cnt_dio_run([days]));
    DetRateper = mean(DetRate*100); 
    text( 4, 3.5,['DetRate(Days 1:' num2str(max(days)) ') =' num2str(round(DetRateper*100)/100) '%'],'FontSize', 24, 'FontWeight','bold');
    
    if saveg1==1,
        orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
        saveas(gcf,['RippleStimRate_AllRunEpochs_tet' num2str(tet) '_SEPstd' num2str(std)],'fig');
        saveas(gcf,['RippleStimRate_AllRunEpochs_tet' num2str(tet) '_SEPstd' num2str(std)],'jpg');
    end
    
end

    





