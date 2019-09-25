
saveg1=0;

%%% Plot Stimulation Calibration Curve from Data Taken Especially for this purpose ###

%directoryname = '/data25/sjadhav/RippleInterruption/SJStimC_direct';
%prefix = 'RE1';

directoryname = '/data25/sjadhav/RippleInterruption/REb_direct/StimAmpl';
%directoryname = '/data25/sjadhav/RippleInterruption/REb_direct';

prefix = 'REb';
days = [17];
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

nstim = 1 % No of stimulation types

cnt=0; cntrip=0;
for d = 1:length(days)
    day = days(d);
    day
    %% Load extracted ripple file
%     ripfile = sprintf('%s/%sripples%02dstd%02d.mat', directoryname, prefix, day, std);
%     load(ripfile);
    
    
    
    for ep = 1:length(allepochs)
        
        epoch = allepochs(ep);
%         rip_starttime = 1000* ripples{day}{epoch}{tet}.starttime;   % in msec
        
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
        multi8 = multi{day}{epoch}{8}/10 ;   %in ms
        multi3 = multi{day}{epoch}{3}/10 ;   %in ms
        multi4 = multi{day}{epoch}{4}/10 ;
        multi12 = multi{day}{epoch}{12}/10 ;
        multi5 = multi{day}{epoch}{5}/10 ;
        multi10 = multi{day}{epoch}{10}/10 ;
        
        
        
        multiu=multi8;
        %multiu=[multi4; multi12; multi5];
        pret=200; postt=400;
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
         
         %%%% Divide stimulations by amplitude: 6 types in all -
         %%%% 0:20uA:100uA 
   
         
         per = floor(size(stim_spkshist,1)/nstim);  %eg. 155/4= 38 pulses per Stim-Ampl
                
         Stimhist_matr=[];
         for n=1:nstim
             
             temp = stim_spkshist( per*(n-1)+1:per*n,:);
             %temp(:,21)= temp(:,20);
             temp(:,(pret/binsize)+1) = 0;
             if binsize==5
                 temp(:,(pret/binsize)+2) = 0;
                 temp(:,(pret/binsize)+3) = 0;
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
axis([0 ((pret+postt)/binsize)+2 0.5 nstim+0.5]);
ypts = 0:1:nstim+1;
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
for n=1:nstim
    cmd=sprintf('stimhist = stimhist%d;',n); eval(cmd);
    subplot(nstim,1,(nstim-n+1)); hold on;
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
    
    if n==nstim,
        title(['Multiunit Firing around stimulation-Ampl 100-400uA 200us'],...
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
taxis = [0:binsize:pret+postt];
plot(taxis, yplot,['r.-'],'Linewidth',2,'Markersize',12);
bar(taxis, yplot,'r');
%plot(taxis, 2*ripnostim_pre(i,:),['r-'],'Linewidth',2,'Markersize',6);

ypts = 0:1.1*max(yplot);
xpts = pret*ones(size(ypts));
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

% figure; hold on;
% redimscreen_land;
% orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
% 
% yplot = (1000/binsize_plot)*mean(rip_spkshist,1); %Multiunit fr in binsize ms bins
% taxis = [0:binsize:600];
% plot(taxis, yplot,['r.-'],'Linewidth',2,'Markersize',12);
% bar(taxis, yplot,'r');
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
% 






