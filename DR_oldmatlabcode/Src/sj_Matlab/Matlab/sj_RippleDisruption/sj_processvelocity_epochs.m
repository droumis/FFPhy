
function [ep] = sj_processvelocity_epochs (animdirect,prefix,days,allepochs,type,figopt1,saveg1,savedata)

%% Generate and Plot Still and Moving Proportions during run OR sleep epochs for all days
% sj_processvelocity_epochs('/data25/sjadhav/RippleInterruption/RE1_direct','RE1',1:9,[2 4],'Run',1,0,0);
% sj_processvelocity_epochs('/data25/sjadhav/RippleInterruption/RE1_direct','RE1',1:9,[1 3 5],'Sleep',1,0,0);
%%
if nargin<6,
    figopt1 = 0;
end
if nargin<7,
    saveg1 = 0;
end
if nargin<8,
    savedata = 0;
end

%% Fixed parameters
directoryname = animdirect;
clr = {'b','g','c','m','y','k','r'};
thrsvel = 5;  % <5cm per sec is still during run epochs (<3cm for sleep epochs)
Fs=30; %video sampling rate
tsamp=1/30; % in sec
velfiltlth = Fs; % Filter for smoothing velocity - Length 1 sec=Fs points, Std Dev = Lth/4


%% Loop over days

for d = 1:length(days)
    
    day = days(d);
    currdayvel = [];
    % Pos file
    posfile = sprintf('%s/%spos%02d.mat', directoryname, prefix, day);
    load(posfile);
    
    % Loop over epochs
    for ep = 1:length(allepochs)
        
        epoch = allepochs(ep);
        
        currvel = pos{day}{epoch}.data(:,5);
        % Look at instantaneous smoothed velocity
        filt = gaussian(velfiltlth/4,velfiltlth);
        smooth_currvel = smoothvect(currvel,filt);
        % Get velocity over a second
        tempvel = smooth_currvel(1:Fs*floor(length(currvel)/Fs));
        x = mean(reshape(tempvel,Fs,floor(length(tempvel)/Fs)),1);
        
        %% Calculate for epochs
        currtime =  pos{day}{epoch}.data(:,1);
        epochtime(d,ep) = currtime(end)-currtime(1); % in sec
        remove_ind = find(currvel < 0);
        meanvel_epoch(d,ep) = mean(x);
        errvel_epoch(d,ep) = std(x);
        
        % Separate into still vs. moving
        still_ind_ep=find(x<=thrsvel);
        move_ind_ep=find(x>thrsvel);
        still_time_epoch(d,ep)=length(still_ind_ep)*1;  %% Each bin is now 1 second
        move_time_epoch(d,ep)=length(move_ind_ep)*1;
        meanvelmove_epoch(d,ep) = mean(x(move_ind_ep));
        errvelmove_epoch(d,ep) = std(x(move_ind_ep));
        
        % Save
        allvel_epoch{d}{ep} = x;
        currdayvel = [currdayvel, x];
        
    end
    
    % Collapse for a day also
    
    totaltime(d) = sum(epochtime(d,:)); % in sec
    meanvel(d) = mean(currdayvel);
    errvel(d) = std(currdayvel);
    allvel{d} = currdayvel;
    
    % Separate into still vs. moving
    still_ind=find(currdayvel<=thrsvel);
    move_ind=find(currdayvel>thrsvel);
    still_time(d)=length(still_ind)*1;
    move_time(d)=length(move_ind)*1;
    meanvelmove(d) = mean(currdayvel(move_ind));
    errvelmove(d) = std(currdayvel(move_ind));
    
end


%% Plotting  %%%%%%%%

if figopt1 ==1
    
    figure; hold on;
    redimscreen_halfvert(0);
    orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
    
    %plot(meanvel, ['r.'], 'Linewidth',2, 'Markersize',12);
    errorbar(days,meanvel,errvel, 'r.-', 'Linewidth',2, 'Markersize',12);
    title(['Mean Velocity (per sec) during ' type ' epochs'],...
        'FontSize',24,'Fontweight','bold');
    axis([min(days)-0.5 max(days)+0.5 0 30]);
    xlabel('Days');
    ylabel('Mean Vel (cm/sec)');
    
    %[h,p,ci] = ttest2(n_pre, n_post, 0.05, 'left' );
    text( 2, 25,['Mean Velocity =' num2str(round(mean(meanvel)*10)/10) '+-' num2str(round(sem(meanvel)*10)/10)],'FontSize', 24, 'FontWeight','bold');
    
    if saveg1==1,
        orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
        saveas(gcf,[type '_MeanVel'],'fig');
        saveas(gcf,[type '_MeanVel'],'jpg');
    end
    
    
    %%%% Plot Velocity for Moving EPochs%%%%%
    
    figure; hold on;
    redimscreen_halfvert(0);
    orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
    
    %plot(meanvel, ['r.'], 'Linewidth',2, 'Markersize',12);
    errorbar(days,meanvelmove,errvelmove, 'r.-', 'Linewidth',2, 'Markersize',12);
    
    title(['Mean Velocity (per sec) during MOVE'],...
        'FontSize',24,'Fontweight','bold');
    axis([min(days)-0.5 max(days)+0.5 0 40]);
    xlabel('Days');
    ylabel('Mean Vel (cm/sec)');
    
    %[h,p,ci] = ttest2(n_pre, n_post, 0.05, 'left' );
    text( 2, 35,['Mean Velocity =' num2str(round(mean(meanvelmove)*10)/10) '+-' num2str(round(sem(meanvelmove)*10)/10)],'FontSize', 24, 'FontWeight','bold');
    
    if saveg1==1,
        orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
        saveas(gcf,[type '_MeanVel_MOVE'],'fig');
        saveas(gcf,[type '_MeanVel_MOVE'],'jpg');
    end
    
    
    %%%% Plot Still vs Moving times %%%%
    
    still_move = [(still_time')./totaltime' (move_time')./totaltime'];
    
    figure; hold on;
    redimscreen_halfvert(0);
    orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
    
    bar(still_move,'stack')
    title(['Still and Moving Times for ' type ' Epochs'],...
        'FontSize',24,'Fontweight','bold');
    axis([min(days)-0.5 max(days)+0.5 0 1.1]);
    xlabel('Day');
    ylabel('Fraction of total time (sec)');
    
    %[h,p,ci] = ttest2(n_pre, n_post, 0.05, 'left' );
    %text( 4, 2.450,['DetRate(4,6,7) =' num2str(round(DetRateper*100)/100) '%'],'FontSize', 24, 'FontWeight','bold');
    
    if saveg1==1,
        orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
        saveas(gcf,[type '_StillMoveTime'],'fig');
        saveas(gcf,[type '_StillMoveTime'],'jpg');
    end
    
    
    figure; hold on;
    redimscreen_halfvert(0);
    orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
    
    bar(totaltime);
    title(['Total Time for velocity for ' type ' Epochs'],...
        'FontSize',24,'Fontweight','bold');
    %axis([0.8 2.2 0 1400]);
    xlabel('Day');
    ylabel('Time (sec)');
    
end

%%
if savedata ==1
    savefile = sprintf('%s/ProcessedData/%s_%s_Velocity.mat', animdirect, prefix,type);
    save(savefile);
end





