
saveg1=0;

%%% Plot Stimulation  ###

%directoryname = '/data25/sjadhav/RippleInterruption/SJStimC_direct';
%prefix = 'RE1';
directoryname = '/data25/sjadhav/RippleInterruption/RCa_direct';
prefix = 'RCa';
days = [1];
allepochs = [2 4];
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
    ripfile = sprintf('%s/%sripples%02dstd%02d.mat', directoryname, prefix, day, std);
    load(ripfile);
    
    
    
    for ep = 1:length(allepochs)
        
        epoch = allepochs(ep);
        rip_starttime = 1000* ripples{day}{epoch}{tet}.starttime;   % in msec
        
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
        multi4 = multi{day}{epoch}{4}/10 ;
        multi9 = multi{day}{epoch}{9}/10 ;
%         
        multiu=multi4;
        pret=200; postt=400;
       
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
    
end


%%%%%% Plot Multiunit rate around stimulations  %%%%%%%%

figure; hold on;
redimscreen_land;
orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');

yplot = (1000/binsize_plot)*mean(stim_spkshist,1); %Multiunit fr in binsize ms bins
taxis = [0:binsize:pret+postt];
taxis = taxis - pret;
plot(taxis, yplot,['r.-'],'Linewidth',2,'Markersize',12);
bar(taxis, yplot,'r');
%plot(taxis, 2*ripnostim_pre(i,:),['r-'],'Linewidth',2,'Markersize',6);

ypts = 0:1.1*max(yplot);
xpts = 0*ones(size(ypts));
plot(xpts , ypts, 'k--','Linewidth',2);

title(['Multiunit Firing around stimulation-AllStimAmpl'],...
    'FontSize',24,'Fontweight','bold');
%axis([0 800 -800 600]);
ylabel('Instantaeous Multiunit Firing Rate');
xlabel('Time(ms)');
set(gca,'XLim',[-(pret+binsize) postt])

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
taxis = [0:binsize:pret+postt];
taxis = taxis - pret;
%plot(taxis, yplot,['r.-'],'Linewidth',2,'Markersize',12);
bar(taxis, yplot,'r');
%plot(taxis, 2*ripnostim_pre(i,:),['r-'],'Linewidth',2,'Markersize',6);

ypts = 0:1.1*max(yplot);
xpts = 0*ones(size(ypts));
plot(xpts , ypts, 'k--','Linewidth',2);

title(['Multiunit Firing around ripples during run - Day' num2str(day)],...
    'FontSize',20,'Fontweight','bold');
%axis([0 800 -800 600]);
ylabel('Instantaeous Multiunit Firing Rate','FontSize',16,'Fontweight','bold');
xlabel('Time(ms)','FontSize',16,'Fontweight','bold');
set(gca,'XLim',[-125.5 305.5])


%text( 4, 2450,['DetRate(4,6,7) =' num2str(round(DetRateper*100)/100) '%'],'FontSize', 24, 'FontWeight','bold');

if saveg1==1,
    orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
    saveas(gcf,['MultiFRrip_' num2str(day) '_tet' num2str(tet) '_SEPstd' num2str(std)],'fig');
    saveas(gcf,['MultiFrrip_' num2str(day) '_tet' num2str(tet) '_SEPstd' num2str(std)],'jpg');
end







