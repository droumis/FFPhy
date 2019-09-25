
saveg1=0;

%%% Plot Stimulation Calibration Curve from Data Taken Especially for this purpose ###

%directoryname = '/data25/sjadhav/RippleInterruption/SJStimC_direct';
%prefix = 'RE1';
directoryname = '/data25/sjadhav/RippleInterruption/REb_direct/StimAmpl';
prefix = 'REb';
days = [20];
allepochs = [1];
clr = {'b','g','c','m','y','k','r'};
tet=4;
std=4;
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
    %ripfile = sprintf('%s/%sripples%02dstd%02d.mat', directoryname, prefix, day, std);
    %load(ripfile);
    
    
    for ep = 1:length(allepochs)
        
        epoch = 1;
        %rip_starttime = 1000* ripples{day}{epoch}{tet}.starttime;   % in msec
        
        %%Load dio file
        DIOfile = sprintf('%s/%sDIO%02d.mat', directoryname, prefix, day);
        load(DIOfile);
        stim = DIO{day}{epoch}{15};
        stim_starttime = stim.pulsetimes(:,1)./10; %ms
        stim_endtime = stim.pulsetimes(:,2)./10; %ms
        stim_length = stim.pulselength;
        stim_isi = stim.timesincelast(2:end)./10; %ms
        
        epoch = allepochs(ep);
       
        multifile = sprintf('%s/%smulti%02d.mat', directoryname, prefix, day);
        load(multifile);
        multi2 = multi{day}{epoch}{2}/10 ;   %in ms
        multi3 = multi{day}{epoch}{3}/10 ;   %in ms
        multi4 = multi{day}{epoch}{4}/10 ;
        multi5 = multi{day}{epoch}{5}/10 ;
        multi8 = multi{day}{epoch}{8}/10 ;
        multi9 = multi{day}{epoch}{9}/10 ;
        
        multiu=[multi8];
        pret=200; postt=250;
%         
       
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
        
%          for i =5:length(rip_starttime)-10
%             i;
%             cntrip=cntrip+1;
%             currrip = rip_starttime(i);
%             currspks =  multiu(find( (multiu>=(currrip-pret)) & (multiu<=(currrip+postt)) ));
%             currspks = currspks-(currrip-pret);
%             histspks = histc(currspks,[0:binsize:pret+postt]);
%             rip_spks{cntrip}=currspks;
%             rip_spkshist(cntrip,:) = histspks;
%          end
         





         %%%% Divide stimulations by amplitude
         
          % 16ripdistest1
%          nranges=7; % 0uA to 120 uA
%         
%          [range1] = timetrans({'00:20:45', '00:22:50'},10000,2); %0
%          [range2] = timetrans({'00:39:00', '00:42:00'},10000,2); %20
%          [range3] = timetrans({'00:35:50', '00:38:50'},10000,2); %40
%          [range4] = timetrans({'00:32:40', '00:35:40'},10000,2); %60
%          [range5] = timetrans({'00:29:29', '00:32:29'},10000,2); %80
%          [range6] = timetrans({'00:22:50', '00:25:50'},10000,2); %100
%          [range7] = timetrans({'00:25:58', '00:28:58'},10000,2); %120
         
           % 19ripdistest2
%          nranges=5; % 0uA, 80, 100, 130 150 
%         
%          [range1] = timetrans({'00:50:51', '00:52:51'},10000,2); 
%          [range2] = timetrans({'00:52:51', '00:54:51'},10000,2); 
%          [range3] = timetrans({'00:44:24', '00:46:24'},10000,2); 
%          [range4] = timetrans({'00:46:40', '00:48:40'},10000,2); 
%          [range5] = timetrans({'00:48:50', '00:50:50'},10000,2); 
         
                % 20ripdistest3
         nranges=5; % 0uA, 150, 300, 500 700 
        
         [range1] = timetrans({'00:11:39', '00:13:45'},10000,2); %0
         [range2] = timetrans({'00:09:35', '00:11:35'},10000,2); %150
         [range3] = timetrans({'00:07:16', '00:09:16'},10000,2); %300
         [range4] = timetrans({'00:02:35', '00:04:35'},10000,2); %500
         [range5] = timetrans({'00:04:55', '00:06:56'},10000,2); %700
        
          
                % 21ripdistest4
%          nranges=4; % 0uA, 150, 300, 500 
%         
%          [range1] = timetrans({'00:08:04', '00:10:20'},10000,2); %0
%          [range2] = timetrans({'00:01:26', '00:03:30'},10000,2); %150
%          [range3] = timetrans({'00:03:42', '00:05:46'},10000,2); %300
%          [range4] = timetrans({'00:05:58', '00:08:02'},10000,2); %500
        
         
         % 17riptest
%          nranges=4;
%          [range1] = timetrans({'00:01:37', '00:02:40'},10000,2);
%          [range2] = timetrans({'00:02:48', '00:04:02'},10000,2);
%          [range3] = timetrans({'00:04:14', '00:05:14'},10000,2);
%          [range4] = timetrans({'00:05:23', '00:06:33'},10000,2);
%          
         
         
           % 20Riptesttime
%          nranges=4;
%          [range1] = timetrans({'00:01:50', '00:03:50'},10000,2);
%          [range2] = timetrans({'00:03:52', '00:05:52'},10000,2);
%          [range3] = timetrans({'00:05:55', '00:07:55'},10000,2);
%          [range4] = timetrans({'00:07:57', '00:10:00'},10000,2);

           %23RCaRipamptest_051910
%            nranges=5;
%            [range1] = timetrans({'00:29:25', '00:31:25'},10000,2);
%            [range2] = timetrans({'00:31:26', '00:33:31'},10000,2);
%            [range3] = timetrans({'00:33:42', '00:35:42'},10000,2);
%            [range4] = timetrans({'00:35:50', '00:38:05'},10000,2);
%            [range5] = timetrans({'00:38:13', '00:40:13'},10000,2);
%            
           %26RCaRipamptestn_051910
%            nranges=7;
%            [range1] = timetrans({'00:07:42', '00:09:42'},10000,2);
%            [range2] = timetrans({'00:09:44', '00:11:44'},10000,2);
%            [range3] = timetrans({'00:11:54', '00:13:57'},10000,2);
%            [range4] = timetrans({'00:15:03', '00:17:03'},10000,2);
%            [range5] = timetrans({'00:17:11', '00:19:11'},10000,2);
%            [range6] = timetrans({'00:19:18', '00:21:18'},10000,2);
%            [range7] = timetrans({'00:21:25', '00:23:30'},10000,2);

            %27RCaRipamptestn_052110
%              nranges=7;
%              [range1] = timetrans({'00:03:25', '00:05:25'},10000,2);
%              [range2] = timetrans({'00:05:27', '00:07:27'},10000,2);
%              [range3] = timetrans({'00:07:28', '00:09:28'},10000,2);
%              [range4] = timetrans({'00:09:30', '00:11:35'},10000,2);
%              [range6] = timetrans({'00:11:45', '00:13:45'},10000,2);
%              [range5] = timetrans({'00:13:50', '00:15:50'},10000,2);
%              [range7] = timetrans({'00:15:58', '00:18:09'},10000,2);



         
         stim_starttime_comp = stim_starttime*10; % Convert from ms to 10000 pts in sec
               
         %per = floor(size(stim_spkshist,1)/4);  %155/4= 38 pulses per Stim-Ampl
                
         Stimhist_matr=[];
         for n=1:nranges
             
             cmd=sprintf('startu=range%d(1);',n); eval(cmd);
             cmd=sprintf('endu=range%d(2);',n); eval(cmd);
             currrange = find(stim_starttime_comp>=startu & stim_starttime_comp<=endu);
             
             temp = stim_spkshist( currrange,:);
             %temp(:,21)= temp(:,20);
             if n~=1
                 temp(:,(pret/binsize)+1) = 0;
                 temp(:,(pret/binsize)+2) = 0;
             end
             
             if binsize==5
                 temp(:,(pret/binsize)+2) = 0;
                 %temp(:,(pret/binsize)+3) = 0;
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
axis([0 ((pret+postt)/binsize)+2 0.5 nranges+0.5]);
ypts = 0:1:nranges+1;
xpts = ((pret/binsize)+1)*ones(size(ypts));
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
for n=1:nranges
    cmd=sprintf('stimhist = stimhist%d;',n); eval(cmd);
    subplot(nranges,1,(nranges-n+1)); hold on;
    yplot = (1000/binsize_plot)*stimhist; %Multiunit fr in "binsize"ms bins
    taxis = [0:binsize:pret+postt];
    %plot(taxis, yplot,['r.-'],'Linewidth',2,'Markersize',12);
    bar(taxis, yplot,'r');
    
    ypts = 0:1.1*max(yplot);
    xpts = pret*ones(size(ypts));
    plot(xpts , ypts, 'k--','Linewidth',2);
        
    set(gca,'ytick',[]);
    
    if n~=1,
        set(gca,'xtick',[]);
    end
    
    if n==1,
        
        %ylabel('MU FR');
        xlabel('Time(ms)');
    end
    
    if n==nranges,
        title(['Multiunit Firing around stimulation-Ampl'],...    %250uA 100-400us
        'FontSize',16,'Fontweight','bold');
    end
    
    if saveg1==1,
    orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
    saveas(gcf,['MultiFR_StimAmpl_' num2str(day) '_tet' num2str(tet) '_SEPstd' num2str(std)],'fig');
    saveas(gcf,['MultiFrStimAmpl_' num2str(day) '_tet' num2str(tet) '_SEPstd' num2str(std)],'jpg');
    end

end




%%%%%% Plot Multiunit rate around stimulations  %%%%%%%%

% figure; hold on;
% redimscreen_land;
% orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
% 
% yplot = (1000/binsize_plot)*mean(stim_spkshist,1); %Multiunit fr in binsize ms bins
% taxis = [0:binsize:pret+postt];
% plot(taxis, yplot,['r.-'],'Linewidth',2,'Markersize',12);
% bar(taxis, yplot,'r');
% %plot(taxis, 2*ripnostim_pre(i,:),['r-'],'Linewidth',2,'Markersize',6);
% 
% ypts = 0:1.1*max(yplot);
% xpts = pret*ones(size(ypts));
% plot(xpts , ypts, 'k--','Linewidth',2);
% 
% title(['Multiunit Firing around stimulation-AllStimAmpl'],...
%     'FontSize',24,'Fontweight','bold');
% %axis([0 800 -800 600]);
% ylabel('Instantaeous Multiunit Firing Rate');
% xlabel('Time(ms)');
% 
% %text( 4, 2450,['DetRate(4,6,7) =' num2str(round(DetRateper*100)/100) '%'],'FontSize', 24, 'FontWeight','bold');
% 
% if saveg1==2,
%     orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
%     saveas(gcf,['MultiFRstim_' num2str(day) '_tet' num2str(tet) '_SEPstd' num2str(std)],'fig');
%     saveas(gcf,['MultiFrstim_' num2str(day) '_tet' num2str(tet) '_SEPstd' num2str(std)],'jpg');
% end


%%% Particular amplitude
figure; hold on;
orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
n=1;
cmd=sprintf('stimhist = stimhist%d;',n); eval(cmd);

yplot = (1000/binsize_plot)*stimhist; %Multiunit fr in "binsize"ms bins
taxis = [0:binsize:pret+postt];
%plot(taxis, yplot,['r.-'],'Linewidth',2,'Markersize',12);
bar(taxis, yplot,'r');

ypts = 0:1.1*max(yplot);
xpts = pret*ones(size(ypts));
plot(xpts , ypts, 'k--','Linewidth',2);


%set(gca,'xtick',[]);
ylabel('MU FR');
xlabel('Time(ms)');
title(['Calibration-REb-0uA'],...  
        'FontSize',16,'Fontweight','bold');


    
    
    
    %%%%%%%%%%
    
%%%%%% Plot Multiunit rate around ripples  %%%%%%%%

% %plot(taxis, 2*ripnostim_pre(i,:),['r-'],'Linewidth',2,'Markersize',6);
% 
% ypts = 0:1.1*max(yplot);
% xpts = 200*ones(size(ypts));
% plot(xpts , ypts, 'k--','Linewidth',2);
% 
% title(['Sleep Multiunit Firing around ripples-Day' num2str(day)],...
%     'FontSize',24,'Fontweight','bold');
% %axis([0 800 -800 600]);
% ylabel('Instantaeous Multiunit Firing Rate');
% xlabel('Time(ms)');
% 
% %text( 4, 2450,['DetRate(4,6,7) =' num2str(round(DetRateper*100)/100) '%'],'FontSize', 24, 'FontWeight','bold');
% 
% if saveg1==1,
%     orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
%     saveas(gcf,['MultiFRrip_' num2str(day) '_tet' num2str(tet) '_SEPstd' num2str(std)],'fig');
%     saveas(gcf,['MultiFrrip_' num2str(day) '_tet' num2str(tet) '_SEPstd' num2str(std)],'jpg');
% end







