
function [ep] = sj_stimresp6_check_ripplesfromextract (animdirect,prefix,day,allepochs,tet,std,saveg1)

%%%% PLOT EXAMPLE EXTRACTED RIPPLES ON GIVEN DAY IN GIVEN EPOCHS 
% eg
% sj_stimresp6_check_ripplesfromextract('/data25/sjadhav/RippleInterruption/RE1_direct','RE1',3,[1 2 3 4 5],5,3);
%std=standard deviation

if nargin<7,
    saveg1 = 0;
end

%%% Plot example extracted ripples in given epochs ###

% directoryname = '/data25/sjadhav/SJStimC_direct';
% prefix = 'sjc';
% directoryname = '/data25/sjadhav/RE1_direct';
% prefix = 'RE1';
% days = 4;  allepochs = [1 2 3];
% tet=5;
% std = 3;

clr = {'b','g','c','m','y','k','r'};

eeg_pre = []; eegnostim_pre = []; rip_pre=[]; ripnostim_pre=[];
eeg_post = []; eegnostim_post = []; rip_post=[]; ripnostim_post=[];
eeg_run = []; eegnostim_run = []; rip_run=[]; ripnostim_run=[];
eeg_run2 = []; eegnostim_run2 = []; rip_run2=[]; ripnostim_run2=[];
eeg_post2 = []; eegnostim_post2 = []; rip_post2=[]; ripnostim_post2=[];

dio_pre = [];
s_pre = []; s_run = []; s_post = []; s_run2=[]; s_post2=[];

t_pre = []; t_post = []; t_run=[]; t_run2=[]; t_post2=[];

%day = days(d);
%d=1;
directoryname = animdirect;

%% Load extracted ripple file

ripfile = sprintf('%s/%sripples%02dstd%02d.mat', directoryname, prefix, day, std);
load(ripfile);

%%Load dio file
DIOfile = sprintf('%s/%sDIO%02d.mat', directoryname, prefix, day);
load(DIOfile);

for ep = 1:length(allepochs)
    
    epoch = allepochs(ep);
    rip_midtime = ripples{day}{epoch}{tet}.midtime;
    
    %%Load eeg files
    EEGfile = sprintf('%s/EEG/%seeg%02d-%01d-%02d.mat', directoryname, prefix, day,epoch,tet);
    load(EEGfile);
    e = eeg{day}{epoch}{tet};
    t = geteegtimes(e);
    eind = lookup(rip_midtime, t);
    
    EEGnostimfile = sprintf('%s/EEG/%seegnostim%02d-%01d-%02d.mat', directoryname, prefix, day,epoch,tet);
    load(EEGnostimfile);
    enostim = eeg{day}{epoch}{tet}.data;
    
    %%Load ripple files
    ripfile = sprintf('%s/EEG/%sripple%02d-%01d-%02d.mat', directoryname, prefix, day,epoch,tet);
    load(ripfile);
    ripamp = ripple{day}{epoch}{tet}.data(:,1);
    ripenv = ripple{day}{epoch}{tet}.data(:,3);
    
    ripnostimfile = sprintf('%s/EEG/%sripplenostim%02d-%01d-%02d.mat', directoryname, prefix, day,epoch,tet);
    load(ripnostimfile);
    ripampnostim = ripple{day}{epoch}{tet}.data(:,1);
    ripenvnostim = ripple{day}{epoch}{tet}.data(:,3);
    
    
    %%%%%%%%%%%%%%%%
    
    
    cnt=0;
    etmp=[]; enostimtmp=[]; ripamptmp=[]; ripampnostimtmp=[];
    
    for i=1:5
        cnt=cnt+1;
        etmp = e.data(eind(i+5)-600:eind(i+5)+600-1);
        enostimtmp = enostim(eind(i+5)-600:eind(i+5)+600-1);
        ripamptmp = ripamp(eind(i+5)-600:eind(i+5)+600-1);
        ripampnostimtmp = ripampnostim(eind(i+5)-600:eind(i+5)+600-1);
        
        %     nrip = length(ripples{day}{epoch}{tet}.startind);
        %     all_trip = ripples{day}{epoch}{tet}.endtime - ripples{day}{epoch}{tet}.starttime;
        %     srip = ripples{day}{epoch}{tet}.maxthresh;
        %     stim = DIO{day}{epoch}{15};
        %     ndio = length(stim.pulselength);
        
        if epoch==1
            eeg_pre = [eeg_pre; etmp'];
            eegnostim_pre = [eegnostim_pre; enostimtmp'];
            rip_pre = [rip_pre; ripamptmp'];
            ripnostim_pre = [ripnostim_pre; ripampnostimtmp'];
        end
        
        if epoch==2
            eeg_run = [eeg_run; etmp'];
            eegnostim_run = [eegnostim_run; enostimtmp'];
            rip_run = [rip_run; ripamptmp'];
            ripnostim_run = [ripnostim_run; ripampnostimtmp'];
        end
        
        if epoch==3
            eeg_post = [eeg_post; etmp'];
            eegnostim_post = [eegnostim_post; enostimtmp'];
            rip_post = [rip_post; ripamptmp'];
            ripnostim_post = [ripnostim_post; ripampnostimtmp'];
        end
        
        if epoch==4
            eeg_run2 = [eeg_run2; etmp'];
            eegnostim_run2 = [eegnostim_run2; enostimtmp'];
            rip_run2 = [rip_run2; ripamptmp'];
            ripnostim_run2 = [ripnostim_run2; ripampnostimtmp'];
        end
        
        if epoch==5
            eeg_post2 = [eeg_post2; etmp'];
            eegnostim_post2 = [eegnostim_post2; enostimtmp'];
            rip_post2 = [rip_post2; ripamptmp'];
            ripnostim_post2 = [ripnostim_post2; ripampnostimtmp'];
        end
    end
end


%%%%%% Plot Example DataSets  %%%%%%%%

tmp=eeg_pre(1,:);
taxis = [1:size(eeg_pre,2)]*1000/e.samprate;
taxis = taxis - 600*1000/e.samprate;

% for ep = 1:length(allepochs)
%     epoch = allepochs(ep);
% end



%%%% Epoch 1 - Sleep 1: %%%%
for i=5,
    figure; hold on;
    redimscreen_land;
    orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
    plot(taxis, eegnostim_pre(i,:),['k-'],'Linewidth',2,'Markersize',6);
    plot(taxis, 2*ripnostim_pre(i,:),['r-'],'Linewidth',2,'Markersize',6);
    
    title(['Pre-Sleep1:Ex. Extracted ripple (Day:' num2str(day) ', SEP sd:' num2str(std) ') on Tet ' num2str(tet)],...
        'FontSize',24,'Fontweight','bold');
    axis([-400 400 -800 600]);
    ylabel('uV / X2uV');
    xlabel('Time(ms)');
    
    %text( 4, 2450,['DetRate(4,6,7) =' num2str(round(DetRateper*100)/100) '%'],'FontSize', 24, 'FontWeight','bold');
    
    if saveg1==1,
        orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
        saveas(gcf,[prefix 'RippleExSleep1_day_' num2str(day) '_tet' num2str(tet) '_std' num2str(std)],'fig');
        saveas(gcf,[prefix 'RippleExSleep1_day_' num2str(day) '_tet' num2str(tet) '_std' num2str(std)],'jpg');
    end
end

%%%% Epoch 2 - Run 1: %%%%
if length(allepochs)>1
    for i=2:5
        figure; hold on;
        redimscreen_land;
        orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
        plot(taxis, eeg_run(i,:),['k-'],'Linewidth',2,'Markersize',6);
        plot(taxis, 2*rip_run(i,:),['r-'],'Linewidth',2,'Markersize',6);
        
        title(['Run1:Ex. Extracted ripple (Day:' num2str(day) ', SEP sd:' num2str(std) ') on Tet ' num2str(tet)],...
            'FontSize',24,'Fontweight','bold');
        axis([-400 400 -800 600]);
        ylabel('uV / X2uV');
        xlabel('Time(ms)');
        
        %text( 4, 2450,['DetRate(4,6,7) =' num2str(round(DetRateper*100)/100) '%'],'FontSize', 24, 'FontWeight','bold');
        
        if saveg1==1,
            orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
            saveas(gcf,[prefix 'RippleEx3Run1_day_' num2str(day) '_tet' num2str(tet) '_SEPstd' num2str(std)],'fig');
            saveas(gcf,[prefix 'RippleEx3Run1_day_' num2str(day) '_tet' num2str(tet) '_SEPstd' num2str(std)],'jpg');
        end
    end
end


%%%% Epoch 3 - Sleep 2: %%%%
if length(allepochs)>2
    for i=2
        figure; hold on;
        redimscreen_land;
        orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
        plot(taxis, eegnostim_post(i,:),['k-'],'Linewidth',2,'Markersize',6);
        plot(taxis, 2*ripnostim_post(i,:),['r-'],'Linewidth',2,'Markersize',6);
        
        title(['Post-Sleep2:Ex. Extracted ripple (Day:' num2str(day) ', SEP sd:' num2str(std) ') on Tet ' num2str(tet)],...
            'FontSize',24,'Fontweight','bold');
        axis([-400 400 -800 600]);
        ylabel('uV / X2uV');
        xlabel('Time(ms)');
        
        %text( 4, 2450,['DetRate(4,6,7) =' num2str(round(DetRateper*100)/100) '%'],'FontSize', 24, 'FontWeight','bold');
        
        if saveg1==1,
            orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
            saveas(gcf,[prefix 'RippleExSleep2_day_' num2str(day) '_tet' num2str(tet) '_SEPstd' num2str(std)],'fig');
            saveas(gcf,[prefix 'RippleExSleep2_day_' num2str(day) '_tet' num2str(tet) '_SEPstd' num2str(std)],'jpg');
        end
    end
end

%%%% Epoch 4 - Run 2: %%%%
if length(allepochs)>3
    for i=2:5
        figure; hold on;
        redimscreen_land;
        orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
        plot(taxis, eeg_run2(i,:),['k-'],'Linewidth',2,'Markersize',6);
        plot(taxis, 2*rip_run2(i,:),['r-'],'Linewidth',2,'Markersize',6);
        
        title(['Run2:Ex. Extracted ripple (Day:' num2str(day) ', SEP sd:' num2str(std) ') on Tet ' num2str(tet)],...
            'FontSize',24,'Fontweight','bold');
        axis([-400 400 -800 600]);
        ylabel('uV / X2uV');
        xlabel('Time(ms)');
        
        %text( 4, 2450,['DetRate(4,6,7) =' num2str(round(DetRateper*100)/100) '%'],'FontSize', 24, 'FontWeight','bold');
        
        if saveg1==1,
            orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
            saveas(gcf,[prefix 'RippleExRun2_day_' num2str(day) '_tet' num2str(tet) '_SEPstd' num2str(std)],'fig');
            saveas(gcf,[prefix 'RippleExRun2_day_' num2str(day) '_tet' num2str(tet) '_SEPstd' num2str(std)],'jpg');
        end
    end
end

%%%% Epoch 5 - Sleep 3: %%%%
if length(allepochs)>4
    for i=2
        figure; hold on;
        redimscreen_land;
        orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
        plot(taxis, eegnostim_post2(i,:),['k-'],'Linewidth',2,'Markersize',6);
        plot(taxis, 2*ripnostim_post2(i,:),['r-'],'Linewidth',2,'Markersize',6);
        
        title(['Post-Sleep3:Ex. Extracted ripple (Day:' num2str(day) ', SEP sd:' num2str(std) ') on Tet ' num2str(tet)],...
            'FontSize',24,'Fontweight','bold');
        axis([-400 400 -800 600]);
        ylabel('uV / X2uV');
        xlabel('Time(ms)');
        
        %text( 4, 2450,['DetRate(4,6,7) =' num2str(round(DetRateper*100)/100) '%'],'FontSize', 24, 'FontWeight','bold');
        
        if saveg1==1,
            orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
            saveas(gcf,[prefix 'RippleExSleep3_day_' num2str(day) '_tet' num2str(tet) '_SEPstd' num2str(std)],'fig');
            saveas(gcf,[prefix 'RippleExSleep3_day_' num2str(day) '_tet' num2str(tet) '_SEPstd' num2str(std)],'jpg');
        end
    end
end

taxis;





% %%%%%%%%%%%%%%% %%%%%%%%


