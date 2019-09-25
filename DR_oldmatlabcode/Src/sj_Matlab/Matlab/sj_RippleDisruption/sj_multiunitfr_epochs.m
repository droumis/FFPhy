function [ep] = sj_multiunitfr_epochs (animdirect,prefix,days,allepochs,type,tet,std,figopt1,saveg1,figopt2,saveg2, savedata)

%% Generate and Plot Still and Moving Proportions during run OR sleep epochs for all days
% sj_multiunitfr_epochs('/data25/sjadhav/RippleInterruption/RE1_direct','RE1',2,[2 4],'Run',5,3,1,0,1,0,0);
% sj_multiunitfr_epochs('/data25/sjadhav/RippleInterruption/RE1_direct','RE1',3,[1 3 5],'Sleep',5,3,1,0,1,0,0);
%%
if nargin<8,
    figopt1 = 0; % Plot individual day figures
end
if nargin<9,
    saveg1 = 0;
end
if nargin<10,
    figopt1 = 0; % Plot individual day figures
end
if nargin<11,
    saveg1 = 0;
end
if nargin<12,
    savedata = 0;
end

%% Fixed parameters
directoryname = animdirect;
clr = {'b','g','c','m','y','k','r'};
thrsvel = 5;  % <5cm per sec is still during run epochs (<3cm for sleep epochs)
binsize = 10; %ms
binsize_plot=10; %% If you put binsize_plot=1000, then units are Nspikes/binsize, not inst. firing rate in Hz


%% Loop over days
alldays_stim_spkshist_mat =[];
alldays_rip_spkshist_mat =[];
for d = 1:length(days)
    
    day = days(d);
    day_stim_spkshist = [];
    day_rip_spkshist = [];
    
    % Load extracted ripple file
    ripfile = sprintf('%s/%sripples%02dstd%02d.mat', directoryname, prefix, day, std);
    load(ripfile);
    
    % Load dio file
    DIOfile = sprintf('%s/%sDIO%02d.mat', directoryname, prefix, day);
    load(DIOfile);
    
    % Load multiunit file
    multifile = sprintf('%s/%smulti%02d.mat', directoryname, prefix, day);
    load(multifile);
    
    for ep = 1:length(allepochs)
        
        epoch = allepochs(ep);
        cnt=0; cntrip=0;
        
        % Stimulations
        stim = DIO{day}{epoch}{15};
        stim_starttime = stim.pulsetimes(:,1)./10; %ms
        stim_endtime = stim.pulsetimes(:,2)./10; %ms
        stim_length = stim.pulselength;
        stim_isi = stim.timesincelast(2:end)./10; %ms
        
        % Ripples
        rip_starttime = 1000* ripples{day}{epoch}{tet}.starttime;   % in ms
        multi5 = multi{day}{epoch}{tet}/10; % in ms
        
        for i = 2:length(stim_starttime)-1
            i;
            cnt=cnt+1;
            currstim = stim_starttime(i);
            currspks =  multi5(find( (multi5>=(currstim-200)) & (multi5<=(currstim+400)) ));
            currspks = currspks-(currstim-200);
            histspks = histc(currspks,[0:binsize:600]);
            stim_spks{cnt} = currspks;
            stim_spkshist(cnt,:) = histspks;
        end
        allepochs_stim_spks{d}{ep} = stim_spks;
        allepochs_stim_spkshist{d}{ep} = stim_spkshist;
        day_stim_spkshist = [day_stim_spkshist; stim_spkshist];
        
        for i = 2:length(rip_starttime)-1
            i;
            cnt=cnt+1;
            currrip = rip_starttime(i);
            currspks =  multi5(find( (multi5>=(currrip-200)) & (multi5<=(currrip+400)) ));
            currspks = currspks-(currrip-200);
            histspks = histc(currspks,[0:binsize:600]);
            rip_spks{cnt} = currspks;
            rip_spkshist(cnt,:) = histspks;
        end
        allepochs_rip_spks{d}{ep} = rip_spks;
        allepochs_rip_spkshist{d}{ep} = rip_spkshist;
        day_rip_spkshist = [day_rip_spkshist; rip_spkshist];
        
    end % end epoch
    
    alldays_stim_spkshist{d} = day_stim_spkshist;
    alldays_rip_spkshist{d} = day_rip_spkshist;
    alldays_stim_spkshist_mat = [alldays_stim_spkshist_mat; day_stim_spkshist];
    alldays_rip_spkshist_mat = [alldays_rip_spkshist_mat; day_rip_spkshist];
    
    % Plot for day
    if figopt1==1
        
        figure(day); hold on;
        redimscreen_2versubplots;
        orient(gcf,'landscape'); hold on; set(gcf, 'PaperPositionMode', 'auto');
        
        %subplot(2,1,1); hold on;
        yplot = (1000/binsize_plot)*mean(day_stim_spkshist,1); %Multiunit fr (Hz) in binsize ms bins
        
        % Remove Stimln artifacts
        yplot((200/binsize)+1)=0;
        yplot(50)= yplot(50)/2;
        yplot([53 57 60])= yplot([53 57 60])./1.2;
        
        
        taxis = [0:binsize:600];
        taxis = taxis - 200;
        %plot(taxis, yplot,['r.-'],'Linewidth',2,'Markersize',12);
        bar(taxis, yplot,'r');
        %plot(taxis, 2*ripnostim_pre(i,:),['r-'],'Linewidth',2,'Markersize',6);
        
        ypts = 0:1.1*max(yplot); xpts = 0*ones(size(ypts));
        plot(xpts , ypts, 'k--','Linewidth',2);
        title(['MUStim-' type '-Day' num2str(day)],'FontSize',20,'Fontweight','bold');
        ylabel('Instantaeous MU Firing Rate'); %xlabel('Time(ms)');
        xlabel('Time(ms)');
        set(gca,'XLim',[-210 400])
        
        %subplot(2,1,2); hold on;
        figure(100+day); hold on;
         redimscreen_2versubplots;
        orient(gcf,'landscape'); hold on; set(gcf, 'PaperPositionMode', 'auto');
        yplot = (1000/binsize_plot)*mean(day_rip_spkshist,1); %Multiunit fr (Hz) in binsize ms bins
        taxis = [0:binsize:600];
        taxis = taxis - 200;
        %plot(taxis, yplot,['r.-'],'Linewidth',2,'Markersize',12);
        bar(taxis, yplot,'r');
        %plot(taxis, 2*ripnostim_pre(i,:),['r-'],'Linewidth',2,'Markersize',6);
        
        ypts = 0:1.1*max(yplot); xpts = 0*ones(size(ypts));
        plot(xpts , ypts, 'k--','Linewidth',2);
        title(['MU Firing around ripples during sleep'],'FontSize',20,'Fontweight','bold');
        ylabel('Instantaeous MU Firing Rate','FontSize',16,'Fontweight','bold');
        xlabel('Time(ms)','FontSize',16,'Fontweight','bold');
        set(gca,'XLim',[-210 400])
        
        if saveg1==1,
            orient(gcf,'landscape'); hold on; set(gcf, 'PaperPositionMode', 'auto');
            saveas(gcf,[prefix 'MultiFRstimrip_' type '_' num2str(day) '_tet' num2str(tet) '_SEPstd' num2str(std)],'fig');
            saveas(gcf,[prefix 'MultiFrstimrip_' type '_' num2str(day) '_tet' num2str(tet) '_SEPstd' num2str(std)],'jpg');
        end
        
    end  % end figure
    
end % end day


%%%%%% Plot Across Days  %%%%%%%%


if figopt2==1
    
    figure; hold on;
    redimscreen_2versubplots;
    orient(gcf,'landscape'); hold on; set(gcf, 'PaperPositionMode', 'auto');
    
    subplot(2,1,1); hold on;
    yplot = (1000/binsize_plot)*mean(alldays_stim_spkshist_mat,1); %Multiunit fr (Hz) in binsize ms bins
    taxis = [0:binsize:600];
    taxis = taxis - 200;
    plot(taxis, yplot,['r.-'],'Linewidth',2,'Markersize',12);
    bar(taxis, yplot,'r');
    %plot(taxis, 2*ripnostim_pre(i,:),['r-'],'Linewidth',2,'Markersize',6);
    
    ypts = 0:1.1*max(yplot); xpts = 0*ones(size(ypts));
    plot(xpts , ypts, 'k--','Linewidth',2);
    title(['MUStim-' type '-Days' num2str(min(days)) 'to' num2str(max(days))],'FontSize',20,'Fontweight','bold');
    ylabel('Instantaeous MU Firing Rate'); %xlabel('Time(ms)');
    
    subplot(2,1,2); hold on;
    yplot = (1000/binsize_plot)*mean(alldays_rip_spkshist_mat,1); %Multiunit fr (Hz) in binsize ms bins
    taxis = [0:binsize:600];
    taxis = taxis - 200;
    plot(taxis, yplot,['r.-'],'Linewidth',2,'Markersize',12);
    bar(taxis, yplot,'r');
    %plot(taxis, 2*ripnostim_pre(i,:),['r-'],'Linewidth',2,'Markersize',6);
    
    ypts = 0:1.1*max(yplot); xpts = 0*ones(size(ypts));
    plot(xpts , ypts, 'k--','Linewidth',2);
    title(['MURip'],'FontSize',20,'Fontweight','bold');
    ylabel('Instantaeous MU Firing Rate');
    xlabel('Time(ms)');
    
    if saveg2==1,
        orient(gcf,'landscape'); hold on; set(gcf, 'PaperPositionMode', 'auto');
        saveas(gcf,[prefix 'MultiFRstimrip_' type 'Avg_' num2str(min(days)) 'to' num2str(max(days)) '_tet' num2str(tet) '_SEPstd' num2str(std)],'fig');
        saveas(gcf,[prefix 'MultiFrstimrip_' type 'Avg_' num2str(min(days)) 'to' num2str(max(days)) '_tet' num2str(tet) '_SEPstd' num2str(std)],'jpg');
    end
    
end  % end figure

if savedata==1
    savefile = sprintf('%s/ProcessedData/%sMultiFRstimrip_Days%01dto%01d_Tet%02d.mat', animdirect, prefix, min(days), max(days), tet);
    save(savefile);
end







