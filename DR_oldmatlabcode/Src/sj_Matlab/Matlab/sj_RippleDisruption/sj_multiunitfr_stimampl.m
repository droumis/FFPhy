
% Use latest function sj_multiunitfr_rip_stimampl1.m / sj_rippower_stimampl3

saveg1=0;

%%% Plot Stimulation Calibration Curve from Data Taken Especially for this purpose ###

%directoryname = '/data25/sjadhav/RippleInterruption/SJStimC_direct';
%prefix = 'RE1';
directoryname = '/data25/sjadhav/RippleInterruption/RE1_direct';
prefix = 'RE1';
days = [15];
allepochs = [1];
clr = {'b','g','c','m','y','k','r'};
tet=5;
std = 3;
eeg_pre = []; eegnostim_pre = []; rip_pre=[]; ripnostim_pre=[];
eeg_post = []; eegnostim_post = []; rip_post=[]; ripnostim_post=[];
eeg_run = []; eegnostim_run = []; rip_run=[]; ripnostim_run=[];
dio_pre = [];
s_pre = []; s_run = []; s_post = [];
t_pre = []; t_post = []; t_run=[];

binsize = 10;  %% ms, for MU Fir Rate
binsize_plot=10; %% If you put binsize_plot=1000, then units are Nspikes/binsize, not inst. firing rate in Hz

cnt=0; cntrip=0;
for d = 1:length(days)
    day = days(d);
    day
    %% Load extracted ripple file
    
    ripfile = sprintf('%s/%sripples%02dstd%02d.mat', directoryname, prefix, day, std);
    load(ripfile);
    
    
    for ep = 1:length(allepochs)
        
        epoch = 1;
        %%Load dio file
        DIOfile = sprintf('%s/%sDIO%02d.mat', directoryname, prefix, day);
        load(DIOfile);
        stim = DIO{day}{epoch}{15};
        stim_starttime = stim.pulsetimes(:,1)./10; %ms
        stim_endtime = stim.pulsetimes(:,2)./10; %ms
        stim_length = stim.pulselength;
        stim_isi = stim.timesincelast(2:end)./10; %ms
        
        epoch = allepochs(ep);
        rip_starttime = 1000* ripples{day}{epoch}{tet}.starttime;   % in msec
        multifile = sprintf('%s/%smulti%02d.mat', directoryname, prefix, day);
        load(multifile);
        multi1 = multi{day}{epoch}{1}/10 ;   %in ms
        multi5 = multi{day}{epoch}{5}/10 ;
        multi7 = multi{day}{epoch}{7}/10 ;
%         
       
        for i =1:length(stim_starttime)
            i;
            cnt=cnt+1;
            currstim = stim_starttime(i);
            currspks =  multi5(find( (multi5>=(currstim-200)) & (multi5<=(currstim+400)) ));
            currspks = currspks-(currstim-200);
            histspks = histc(currspks,[0:binsize:600]);
            stim_spks{cnt}=currspks;
            stim_spkshist(cnt,:) = histspks;
        end
        
         for i =5:length(rip_starttime)-10
            i;
            cntrip=cntrip+1;
            currrip = rip_starttime(i);
            currspks =  multi5(find( (multi5>=(currrip-200)) & (multi5<=(currrip+400)) ));
            currspks = currspks-(currrip-200);
            histspks = histc(currspks,[0:binsize:600]);
            rip_spks{cntrip}=currspks;
            rip_spkshist(cntrip,:) = histspks;
         end
         
         %%%% Divide stimulations by amplitude: 8 types in all - 10uA to
         %%%% 80uA
         
         per = size(stim_spkshist,1)/8;  %304/8= 38 pulses per Stim-Ampl
         
         
         Stimhist_matr=[];
         for n=1:8
             
             temp = stim_spkshist( per*(n-1)+1:per*n,:);
             %temp(:,21)= temp(:,20);
             temp(:,(200/binsize)+1) = 0;
             if binsize==5
                 temp(:,(200/binsize)+2) = 0;
                 temp(:,(200/binsize)+3) = 0;
             end
             cmd=sprintf('stimhistall%d = temp;',n); eval(cmd);
             cmd=sprintf('stimhist%d = mean(temp,1);',n); eval(cmd);
             
             Stimhist_matr = [Stimhist_matr; mean(temp,1)];
             
         end

         
        
    end
    
end


%%%%%% Plot Diff Ampl Stim Resp %%%%%

%%% Matrix
Stimhist_matr=Stimhist_matr(:,1:end-1);
figure; hold on;
orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
imagesc(Stimhist_matr);
title(['Multiunit Firing around stimulation'],...
    'FontSize',24,'Fontweight','bold');
ylabel('Stim Ampl');
xlabel('Time(ms)');
axis([0 (600/binsize)+2 0.5 8.5]);
ypts = 0:1:9;
xpts = ((200/binsize)+1)*ones(size(ypts));
plot(xpts , ypts, 'r--','Linewidth',2);

set(gca,'xtick',[0:20:60],'xticklabel',{num2str([-200,0,200,400]')});


if saveg1==1,
    orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
    saveas(gcf,['Image_MultiFR_StimAmpl_' num2str(day) '_tet' num2str(tet) '_SEPstd' num2str(std)],'fig');
    saveas(gcf,['Image_MultiFrStimAmpl_' num2str(day) '_tet' num2str(tet) '_SEPstd' num2str(std)],'jpg');
end



%%% Each Bar Graph
figure; hold on;
redimscreen_widevert(0);
orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
for n=1:8
    cmd=sprintf('stimhist = stimhist%d;',n); eval(cmd);
    subplot(8,1,(8-n+1)); hold on;
    yplot = (1000/binsize_plot)*stimhist; %Multiunit fr in "binsize"ms bins
    taxis = [0:binsize:600];
    %plot(taxis, yplot,['r.-'],'Linewidth',2,'Markersize',12);
    bar(taxis, yplot,'r');
    
    ypts = 0:1.1*max(yplot);
    xpts = 200*ones(size(ypts));
    plot(xpts , ypts, 'k--','Linewidth',2);
        
    set(gca,'ytick',[]);
    
    if n~=1,
        set(gca,'xtick',[]);
    end
    
    if n==1,
        
        %ylabel('MU FR');
        xlabel('Time(ms)');
    end
    
    if n==8,
        title(['Multiunit Firing around stimulation-Ampl 10-80uA'],...
        'FontSize',16,'Fontweight','bold');
    end
    
    if saveg1==1,
    orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
    saveas(gcf,['MultiFR_StimAmpl_' num2str(day) '_tet' num2str(tet) '_SEPstd' num2str(std)],'fig');
    saveas(gcf,['MultiFrStimAmpl_' num2str(day) '_tet' num2str(tet) '_SEPstd' num2str(std)],'jpg');
    end

end




%%%%%% Plot Multiunit rate around stimulations  %%%%%%%%

figure; hold on;
redimscreen_land;
orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');

yplot = (1000/binsize_plot)*mean(stim_spkshist,1); %Multiunit fr in binsize ms bins
taxis = [0:binsize:600];
plot(taxis, yplot,['r.-'],'Linewidth',2,'Markersize',12);
bar(taxis, yplot,'r');
%plot(taxis, 2*ripnostim_pre(i,:),['r-'],'Linewidth',2,'Markersize',6);

ypts = 0:1.1*max(yplot);
xpts = 200*ones(size(ypts));
plot(xpts , ypts, 'k--','Linewidth',2);

title(['Multiunit Firing around stimulation-AllStimAmpl'],...
    'FontSize',24,'Fontweight','bold');
%axis([0 800 -800 600]);
ylabel('Instantaeous Multiunit Firing Rate');
xlabel('Time(ms)');

%text( 4, 2450,['DetRate(4,6,7) =' num2str(round(DetRateper*100)/100) '%'],'FontSize', 24, 'FontWeight','bold');

if saveg1==2,
    orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
    saveas(gcf,['MultiFRstim_' num2str(day) '_tet' num2str(tet) '_SEPstd' num2str(std)],'fig');
    saveas(gcf,['MultiFrstim_' num2str(day) '_tet' num2str(tet) '_SEPstd' num2str(std)],'jpg');
end




%%%%%% Plot Multiunit rate around ripples  %%%%%%%%

figure; hold on;
redimscreen_land;
orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');

yplot = (1000/binsize_plot)*mean(rip_spkshist,1); %Multiunit fr in binsize ms bins
taxis = [0:binsize:600];
plot(taxis, yplot,['r.-'],'Linewidth',2,'Markersize',12);
bar(taxis, yplot,'r');
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
    saveas(gcf,['MultiFRrip_' num2str(day) '_tet' num2str(tet) '_SEPstd' num2str(std)],'fig');
    saveas(gcf,['MultiFrrip_' num2str(day) '_tet' num2str(tet) '_SEPstd' num2str(std)],'jpg');
end










%%%% Sleep 2: %%%%
% for i=5
%     figure; hold on;
%     redimscreen_land;
%     orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
%     plot(taxis, eegnostim_post(i,:),['k-'],'Linewidth',2,'Markersize',6);
%     plot(taxis, 2*ripnostim_post(i,:),['r-'],'Linewidth',2,'Markersize',6);
%     
%     title(['Post-Sleep:Ex. Extracted ripple (Day:' num2str(day) ', SEP sd:' num2str(std) ') on Tet ' num2str(tet)],...
%         'FontSize',24,'Fontweight','bold');
%     axis([0 800 -900 600]);
%     ylabel('uV / X2uV');
%     xlabel('Time(ms)');
%     
%     %text( 4, 2450,['DetRate(4,6,7) =' num2str(round(DetRateper*100)/100) '%'],'FontSize', 24, 'FontWeight','bold');
%     
%     if saveg1==1,
%         orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
%         saveas(gcf,['RippleExSleep2_day_' num2str(day) '_tet' num2str(tet) '_SEPstd' num2str(std)],'fig');
%         saveas(gcf,['RippleExSleep2_day_' num2str(day) '_tet' num2str(tet) '_SEPstd' num2str(std)],'jpg');
%     end
% end
% 
% %%%% Run1: %%%%
% for i=5
%     figure; hold on;
%     redimscreen_land;
%     orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
%     plot(taxis, eegnostim_run(i,:),['k-'],'Linewidth',2,'Markersize',6);
%     plot(taxis, 2*ripnostim_run(i,:),['r-'],'Linewidth',2,'Markersize',6);
%     
%     title(['Run:Ex. Extracted ripple (Day:' num2str(day) ', SEP sd:' num2str(std) ') on Tet ' num2str(tet)],...
%         'FontSize',24,'Fontweight','bold');
%     axis([0 800 -900 900]);
%     ylabel('uV / X2uV');
%     xlabel('Time(ms)');
%     
%     %text( 4, 2450,['DetRate(4,6,7) =' num2str(round(DetRateper*100)/100) '%'],'FontSize', 24, 'FontWeight','bold');
%     
%     if saveg1==1,
%         orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
%         saveas(gcf,['RippleEx3Run_day_' num2str(day) '_tet' num2str(tet) '_SEPstd' num2str(std)],'fig');
%         saveas(gcf,['RippleEx3Run_day_' num2str(day) '_tet' num2str(tet) '_SEPstd' num2str(std)],'jpg');
%     end
% end



% %%%%%%%%%%%%%%% %%%%%%%%


