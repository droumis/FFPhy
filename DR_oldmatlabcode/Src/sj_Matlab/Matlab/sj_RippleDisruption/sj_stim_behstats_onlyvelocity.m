
function [ep] = sj_stim_behstats_onlyvelocity (animdirect,prefix,days,allepochs,figopt1,saveg1,figopt2, saveg2, savedata)

%% PLOT Position and Velocity of Stimulations in run epochs
% eg
% sj_stim_behstats_onlyvelocity('/data25/sjadhav/RippleInterruption/RCb_direct','RCb',[1:8],[2 4],1,0,0);
% sj_stim_behstats_onlyvelocity('/data25/sjadhav/RippleInterruption/REe_direct','REe',[1:8],[2 4],1,0,0);
% sj_stim_behstats_onlyvelocity('/data25/sjadhav/RippleInterruption/REb_direct/StimAmpl','REb',17,[1],1,0,0);
% sj_stim_behstats_onlyvelocity('/data25/sjadhav/RippleInterruption/SJStimC_direct','sjc',1:3,[2 4],1,0,0);

% figopt1: Make plots for posn and vel mapping of stimulations
% saveg1: Save graphs
% figopt2, saveg2: Plots across days
% savedata: Save mat file with data

%%
if nargin<5,
    figopt1 = 0;
end
if nargin<6,
    saveg1 = 0;
end
if nargin<7,
    figopt2 = 0;
end
if nargin<8,
    saveg2 = 0;
end
if nargin<9,
    savedata = 0;
end


%% Fixed parameters
Fs=29.97; %video sampling rate
tsamp=1/Fs; % in sec

% Use Calebs filter
defaultfilter = 'velocitydayprocess_filter.mat';
eval(['load ', defaultfilter]);
L = length(velocityfilter.kernel);


% Variable parameters
thrsvel=5;  %<5cm per sec is still on track (3 for sleep session in sj_stimresp2_withvel.m)
pastsecs = 1; % For looking at history of Beh, eg. velocity: used in stimresp2
thrsdist = 20; %in cm, for defining well and intersection boundary
armthrsdist = 20;
dur = 15; % Standard epoch duration in mins

%% Initialize
directoryname = animdirect;
if (animdirect(end) == '/')
    animdirect = animdirect(1:end-1);
end
cd(animdirect);
clr = {'b','g','c','m','y','k','b','g','c'};
wellclr = {'r','m','g'};


if prefix=='sjc',
    fix=1;
else
    fix=0;
end

% ------------------------------
% Figure and Font Sizes

forppr = 1;
% If yes, everything set to redimscreen_figforppr1
% If not, everything set to redimscreen_figforppt1

figdir = '/data25/sjadhav/RippleInterruption/Figures/01AugSep11_RippleDisFigs/StimulationStats/';
summdir = figdir;
set(0,'defaultaxesfontweight','normal'); set(0,'defaultaxeslinewidth',2);

if forppr==1
    set(0,'defaultaxesfontsize',16);
    tfont = 18; % title font
    xfont = 16;
    yfont = 16;
else
    set(0,'defaultaxesfontsize',24);
    tfont = 28;
    xfont = 20;
    yfont = 20;
end

clr = {'b','r','g','c','m','y','k','r'};
savefig1=0;
% -------------------------------------


%% Loop over days and load files
ratecnt=0;
for d = 1:length(days)
    
    day = days(d);
    % Day n
    DIOfile = sprintf('%s/%sDIO%02d.mat', directoryname, prefix, day);
    load(DIOfile);
    
    % Pos file
    posfile = sprintf('%s/%spos%02d.mat', directoryname, prefix, day);
    load(posfile);
    
    % LinPos file
    linposfile = sprintf('%s/%slinpos%02d.mat', directoryname, prefix, day);
    load(linposfile);
    
    % Figure for each day - epochs in subplot (1 fig for position, 1 fig for velocity)
    if figopt1==1
        figure(day); hold on; redimscreen_2horsubplots;
        
        orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
        
    end
    
    %     if (figopt1==1 || figopt2==1)
    %         figure(11); hold on;
    %         orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
    %     end
    
    % Epoch 2,4... - Run 1,2,...
    % Loop over epochs
    for ep = 1:length(allepochs)
        
        % Get Stimln information
        epoch = allepochs(ep);
        stim = DIO{day}{epoch}{15};
        if isempty(stim),
            stim = DIO{day}{epoch}{16};
        end
        stim_starttime = stim.pulsetimes(:,1);
        stim_endtime = stim.pulsetimes(:,2);
        stim_length = stim.pulselength;
        stim_isi = stim.timesincelast(2:end);
        pt = stim.pulsetimes ./ 10000;
        
        %% Get Well and Intersection positions
        wellpos = linpos{day}{epoch}.wellSegmentInfo.wellCoord;
        
        %%% Get Position and Velocity
        currtime = pos{day}{epoch}.data(:,1);
        totaltime = length(currtime)/Fs;
        
        if size(pos{day}{epoch}.data,2)>5 % already smoothed position and filtered velocity
            currpos = pos{day}{epoch}.data(:,6:7);
            currvel = pos{day}{epoch}.data(:,9);
            smooth_currvel=currvel;
            
        else
            currpos = pos{day}{epoch}.data(:,2:3);
            currvel = pos{day}{epoch}.data(:,5);
            % Look at instantaneous smoothed velocity - Smooth velocity
            smooth_currvel = [filtfilt(velocityfilter.kernel, 1, currvel)];
            %smooth_currvel = smoothvect(currvel,filt); 2nd option for convolution
        end
        
        % Get indices for stimulation time in the pos vector
        posind = lookup(pt(:,1), currtime);
        % Skip first few and last few if needed
        posind=posind(1:end);
        pos_stim = currpos(posind,:);
        
        % look at instantaneous velocity at stimulation times
        velindk = posind;
        for i=1:length(velindk)
            vel_stim(i) = min(smooth_currvel(velindk(i):velindk(i)));
        end
        
        %figure; hold on;
        %plot(smooth_currvel); plot(velindk,vel_stim,'ko');
        
        % Divide into moving stimulations and still stimulations
        moving_stimidx = find(vel_stim>thrsvel); moving_stimvel = vel_stim(find(vel_stim>thrsvel));
        still_stimidx = find(vel_stim<=thrsvel); still_stimvel = vel_stim(find(vel_stim<=thrsvel));
        
        % Divide all velocity into moving and stim
        moving_vel =  smooth_currvel(find( smooth_currvel>thrsvel));
        still_vel =  smooth_currvel(find( smooth_currvel<=thrsvel));
        time_moving = length(moving_vel)/Fs;
        time_still = length(still_vel)/Fs;
        
%         % Save for day
%         day_totalstim(d,ep) = length(pos_stim);
%         day_totaltime(d,ep) = totaltime;
%         day_velstim{d}{ep} = vel_stim;
%         day_vel{d}{ep} = smooth_currvel;
        
        ratecnt=ratecnt+1;
        stimrate(ratecnt)=length(pos_stim)./(dur*60);   % This is per epoch
        
        %%    Figure for Current Day and Epoch
        
        if figopt1==1
            figure(day); hold on;
            % Plot stimulation spreads
            subplot(1,length(allepochs),ep); hold on;
            plot(pos_stim(:,1),pos_stim(:,2),'.')
            
            %axis([-5 150 -5 150]);
            xlabel('X-position (mm)');
            ylabel('Y-position (mm)');
            title(['Posn at Stim: Day' num2str(day) ' - Ep' num2str(allepochs(ep)) '; Nstim: ' num2str(length(pos_stim)) '; Rate: ' num2str(roundn(stimrate(ratecnt))) ' Hz'],'FontSize',20,'Fontweight','normal');
            
        end
        
        xhist = [0:1:floor(max(vel_stim))];
        yhist = hist(vel_stim,xhist);
        
        % Plot velocity distribution at stimulation
        if (figopt2==1)
            figure(day*10+1); hold on; redimscreen;
            subplot(length(allepochs),2,2*(ep-1)+1); hold on;
            plot(xhist,yhist,[clr{d} '.-'],'MarkerSize',16,'LineWidth',2);
            %plot([0:1:round(max(vel_stim))],hist(smooth_currvel,[0:1:round(max(vel_stim))]),['r.-'],'MarkerSize',8,'LineWidth',2);
            plot([0:1:50],hist(smooth_currvel,[0:1:50]),['r.-'],'MarkerSize',8,'LineWidth',2);
            plot(thrsvel*ones(size([0:10:max(hist(smooth_currvel,[0:1:50]))])),[0:10:max(hist(smooth_currvel,[0:1:50]))],'k--','Linewidth',1);
            
            % if d==length(days)
            xlabel('Speed (cm/sec)');
            ylabel(['No of Stimulations']);
            title([prefix ': Vel Distr at Stim: Day ' num2str(day) ' - Ep' num2str(ep)],'FontSize',20,'FontWeight','normal');
            %title([prefix ': Vel Distr at Stim: Day ' num2str(min(days)) ' to ' num2str(max(days)) ' - Ep' num2str(ep)],'FontSize',20,'Fontweight','normal');
            text(20, 0.9*max(hist(smooth_currvel,[0:1:50])),['(Avg Stim Rate = ' num2str(roundn(stimrate(ratecnt))) ' Hz)'],...
                'FontSize',16,'FontWeight','normal');
            % end
            
            subplot(length(allepochs),2,2*(ep-1)+2); hold on;
            norm_yhist = yhist./max(yhist);
            vel_hist = hist(smooth_currvel,[0:1:50]);
            norm_vel_hist = vel_hist ./ max(vel_hist);
            plot(xhist,norm_yhist,[clr{d} '.-'],'MarkerSize',16,'LineWidth',2);
            %plot([0:1:round(max(vel_stim))],hist(smooth_currvel,[0:1:round(max(vel_stim))]),['r.-'],'MarkerSize',8,'LineWidth',2);
            plot([0:1:50],norm_vel_hist,['r.-'],'MarkerSize',8,'LineWidth',2);
            plot(thrsvel*ones(size([0:0.1:1])),[0:0.1:1],'k--','Linewidth',1);
            
            frac_below5 = sum(norm_yhist(1:5)) ./ sum(norm_yhist);
            frac_below5 = roundn(frac_below5,-3);
            per_below5 = frac_below5*100;
            
            frac_below10 = sum(norm_yhist(1:10)) ./ sum(norm_yhist);
            frac_below10 = roundn(frac_below10,-3);
            per_below10 = frac_below10*100;
            
            %if d==length(days)
            xlabel('Speed (cm/sec)');
            ylabel(['Norm No of Stimulations']);
            title([prefix ': Vel Distr at Stim: Day ' num2str(day) ' - Ep' num2str(ep)],'FontSize',20,'Fontweight','normal');
            text(15, 0.9,['(Avg Stim Rate = ' num2str(roundn(stimrate(ratecnt))) ' Hz)'],...
                'FontSize',16,'FontWeight','normal');
            text(15, 0.7,['NStim below 5 cm/sec = ' num2str(per_below5) '%'],'FontSize',14,'Fontweight','normal');
            text(15, 0.55,['NStim below 10cm/sec = ' num2str(per_below10) '%'],'FontSize',14,'Fontweight','normal');
            
            %end
        end
        
        
        %%%%%%%%%%%%%%%% Save for day-epoch %%%%%%%%%%%%%%%%%%%
        dayep_totalstim(d,ep) = length(pos_stim);
        dayep_totaltime(d,ep) = totaltime; % in sec
        %dayep_stimrate(d,ep) = length(vel_stim) / totaltime; % in Hz
        dayep_stimrate(d,ep) = length(pos_stim)./(dur*60);; % in Hz
        
        dayep_velstim{d}{ep} = vel_stim;
        dayep_velstim_moving{d}{ep} = moving_stimvel;
        dayep_velstim_still{d}{ep} = still_stimvel;
        dayep_time_moving(d,ep) = time_moving; % in sec
        dayep_time_still(d,ep) = time_still; % in sec
        
        dayep_vel{d}{ep} = smooth_currvel;
        dayep_vel_moving{d}{ep} = moving_vel;
        dayep_vel_still{d}{ep} = still_vel;
        
    end  %% end epoch
    
    if saveg1==1,
        figure(day); hold on;
        orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
        saveas(gcf,[prefix '_StimlnBehStats' '_Day' num2str(day)],'fig');
        saveas(gcf,[prefix '_StimlnBehStats' '_Day' num2str(day)],'jpg');
    end
    
    
    %%%%%%%%%%%%%%%% Save for Day%%%%%%%%%%%%%%%%%%%
    day_totalstim(d) = mean(dayep_totalstim(d,:),2);
    day_totaltime(d) = mean(dayep_totaltime(d,:),2); % in sec
    day_stimrate(d) = mean( dayep_stimrate(d,:),2); % in Hz
    
    day_velstim{d} = [dayep_velstim{d}{1} dayep_velstim{d}{2}];
    day_velstim_moving{d} = [dayep_velstim_moving{d}{1} dayep_velstim_moving{d}{2}];
    day_velstim_still{d} = [dayep_velstim_still{d}{1} dayep_velstim_still{d}{2}];
    day_time_moving(d) = mean(dayep_time_moving(d,:),2); % in sec
    day_time_still(d) = mean(dayep_time_still(d,:),2); % in sec
    
    day_vel{d} = [dayep_vel{d}{1}; dayep_vel{d}{2}];
    day_vel_moving{d} = [dayep_vel_moving{d}{1}; dayep_vel_moving{d}{2}];
    day_vel_still{d} = [dayep_vel_still{d}{1}; dayep_vel_still{d}{2}];
    
end %% end day


day_stimrate = mean(reshape(stimrate,2,8),1);
mean_stimrate = mean(day_stimrate);


%% Save Vel figure

if saveg2==1,
    figure(11); hold on; redimscreen_figforppt1;
    saveas(gcf,[prefix '_StimlnBehVelStats_AcrossDays'],'fig');
    saveas(gcf,[prefix '_StimlnBehVelStats_AcrossDays'],'jpg');
end


%% Plot across day comparisons
% if figopt2==1
%     figure(21); hold on;
%     redimscreen_2versubplots
%     orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
%
%     wellintvsarm = []; well1=[]; well2 = []; well3=[];
%     for d = 1:length(days)
%
%         day = days(d);
%
%         % Well and intersection vs. arm across days
%         currday  = [sum(sum(nstimwell(d,:,:),3)) + sum(sum(nstimint(d,:,:),3)), sum(sum(nstimarm(d,:,:),3))];
%         wellintvsarm = [wellintvsarm; currday];
%
%         % Wells across days
%         well1 = [well1, sum(nstimwell(d,:,1))];
%         well2 = [well2, sum(nstimwell(d,:,2))];
%         well3 = [well3, sum(nstimwell(d,:,3))];
%
%     end
%
%     % Well and intersection vs. arm
%     subplot(2,1,1); hold on;
%     bar(wellintvsarm,'group');
%     title([prefix ': Arm stimlns vs. Well+Interstn stimlns']);
%     legend('All Wells+Interstns','All Arms','Location','NorthEastOutside');
%     xlabel(['Days']); ylabel(['No of stimulations']);
%     %set(gca,'xtick',[1 2 3],'xticklabel',{'2';'1';'3'});
%     %set(gca,'xtick',[days],'xticklabel',{num2str(days)'});
%
%     % Wells across days
%     subplot(2,1,2); hold on;
%     Y = [well2; well1; well3];
%     bar(Y,'group');
%     title([prefix ': No. of well stimulns across days']);
%     %legend('All Wells+Interstns','All Arms','Location','NorthEastOutside');
%     xlabel(['Wells']); ylabel(['No of stimulations']);
%     set(gca,'xtick',[1 2 3],'xticklabel',{'2';'1';'3'});
% end
%
% if saveg2==1,
%     figure(21); hold on;
%     orient(gcf,'landscape'); hold on; set(gcf, 'PaperPositionMode', 'auto');
%     saveas(gcf,[prefix '_StimlnBehStats_AcrossDays'],'fig');
%     saveas(gcf,[prefix '_StimlnBehStats_AcrossDays'],'jpg');
% end
%
% if savedata==1,
%
%     savefile = sprintf('%s/ProcessedData/%s_StimlnBehStats.mat', animdirect, prefix);
%     %savefile = [prefix '_LTP_Days' num2str(min(days)) ':' num2str(max(days)) 'Tet' num2str(tet)];
%     save(savefile);
%
% end
%
%
