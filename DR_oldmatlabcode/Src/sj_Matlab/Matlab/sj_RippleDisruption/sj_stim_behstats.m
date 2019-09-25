
function [ep] = sj_stim_behstats (animdirect,prefix,days,allepochs,figopt1,figopt2, figopt3, saveg1, saveg2, savedata)

%% Save stimstats file
%% Get and plot Position and Velocity of Stimulations in run epochs for
%% given animal - Borrow from sj_stim_behstats_onlyvelocity
% eg
% sj_stim_behstats('/data25/sjadhav/RippleInterruption/RE1_direct','RE1',1:10,[2 4]);
% sj_stim_behstats('/data25/sjadhav/RippleInterruption/RCb_direct','RCb',1:8,[2 4],0,0,0);

% figopt1,2: Make plots for posn(1) and vel(2) mapping of stimulations
% saveg1,2: Save graphs
% figopt3, saveg2: Plots across days
% savedata: Save mat file with data

%%
if nargin<5,
    figopt1 = 0;
end
if nargin<6,
    figopt2 = 0;
end
if nargin<7,
    figopt3 = 0;
end
if nargin<8,
    saveg1 = 0;
end
if nargin<9,
    saveg2 = 0;
end
if nargin<10,
    savedata = 0;
end

%% Fixed parameters
Fs=29.97; %video sampling rate
tsamp=1/Fs; % in sec
lockout=0.25; % in sec
dur = 15*60; % Standard epoch duration in secs


% Use Calebs speed filter if necessary
defaultfilter = 'velocitydayprocess_filter.mat';
eval(['load ', defaultfilter]);
L = length(velocityfilter.kernel);

% Variable parameters
thrsvel=5;  %<5cm per sec is still on track (3/2 for sleep session in sj_stimresp2_withvel.m)
pastsecs = 1; % For looking at history of Beh, eg. velocity: used in stimresp2
thrsdistwell = 20; %in cm, for defining well boundary
thrsdistint = 15; %in cm, for defining intersection boundary
armthrsdist = 20;

%% Initialize
directoryname = animdirect;
if (animdirect(end) == '/')
    animdirect = animdirect(1:end-1);
end
cd(animdirect);

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
wellclr = {'r','m','g'};

% ---------------------------------------

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
        figure(day); hold on; % for position at stim
        if prefix == 'sjc'
            redimscreen_2versubplots
            orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
        else
            redimscreen_2x2subplots
            orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
        end
    end
    
    if (figopt2==1)
        figure(21); hold on; % for vel
        %redimscreen_portrait
        if prefix=='sjc'
            redimscreen_2versubplots
            orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
        else
            redimscreen_2horsubplots
            orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
        end
    end
    
    % Epoch 2,4... - Run 1,2,...
    % Loop over epochs
    for ep = 1:length(allepochs)
        
        % Get Stimln information
        % ------------------------
        epoch = allepochs(ep);
        stim = DIO{day}{epoch}{16};
        if isempty(stim),
            stim = DIO{day}{epoch}{15};
        end
        stim_starttime = stim.pulsetimes(:,1); %in tens of ms
        stim_endtime = stim.pulsetimes(:,2);
        stim_length = stim.pulselength;
        stim_isi = stim.timesincelast(2:end);
        % Get rid of errors
        stim_isi=diff(stim_starttime);
        rem=find(stim_isi<=lockout*10000);
        stim_starttime(rem+1)=[];
        stim_isi=diff(stim_starttime)./10000; % in sec
        pt = stim.pulsetimes ./ 10000;    % in sec
    
        % Get Position and Velocity
        % ----------------------------
        currtime = pos{day}{epoch}.data(:,1);
        totaltime = length(currtime)/Fs; % in sec. This will be inaccurate-limited by Fs. ~100ms 
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
        
        % Stimulation indices in pos and vel vector
        % ------------------------------------------
        % Get indices for stimulation time in the pos vector
        posind = lookup(pt(:,1), currtime);
        pos_stim = currpos(posind,:);
        % Skip first few and last few if needed
        posind=posind(1:end);
        pos_stim = currpos(posind,:);
        % look at instantaneous velocity at stimulation times
        velindk = posind; npre=1;npost=1;
        for i=1:length(velindk)
            % if ((velindk(i)-npre>0) & (velindk(i)+npost<length(smooth_currvel)))
            % vel_stim(i) = min(smooth_currvel(velindk(i)-npre:velindk(i)+npost));
            % else
            % vel_stim(i) = min(smooth_currvel(velindk(i):velindk(i)));
            % end
            vel_stim(i) = smooth_currvel(velindk(i));
            %if vel_stim(i)>20, vel_stim(i)=vel_stim(i)-20; end
        end
        
        % Vel-Moving vs Still
        % --------------------
        % Divide into moving stimulations and still stimulations
        moving_stimidx = find(vel_stim>thrsvel); moving_stimvel = vel_stim(find(vel_stim>thrsvel));
        still_stimidx = find(vel_stim<=thrsvel); still_stimvel = vel_stim(find(vel_stim<=thrsvel));
        % Divide all velocity into moving and still
        moving_vel =  smooth_currvel(find( smooth_currvel>thrsvel));
        still_vel =  smooth_currvel(find( smooth_currvel<=thrsvel));
        time_moving = length(moving_vel)/Fs;
        time_still = length(still_vel)/Fs;
   
        % Stim rate - this is no longer used
        % ----------
        ratecnt=ratecnt+1; % each session (d,ep) is a measure
        stimrate(ratecnt)=length(pos_stim)./(totaltime);
        
        % Position at stimulation calculations
        % -------------------------------------
        % Get Well and Intersection positions
        wellpos = linpos{day}{epoch}.wellSegmentInfo.wellCoord;
        intpos(1,:) = linpos{day}{epoch}.segmentInfo.segmentCoords(1,3:4);
        intpos(2,:) = linpos{day}{epoch}.segmentInfo.segmentCoords(2,3:4);
        intpos(3,:) = linpos{day}{epoch}.segmentInfo.segmentCoords(4,3:4);
        % Find number of stimulations at well, at intersection and center-arm run
        sstimarm=[]; pstimarm=[];
        for i=1:3
            currwellpos = wellpos(i,:);
            currintpos = intpos(i,:);
            nstimwell(d,ep,i) = length(find(dist(pos_stim,repmat(currwellpos,length(pos_stim),1)) <= thrsdistwell));
            nstimint(d,ep,i) = length(find(dist(pos_stim,repmat(currintpos,length(pos_stim),1)) <= thrsdistint));
            
            switch i
                case 2+fix
                    stimarm = pos_stim(find(pos_stim(:,1)<currwellpos(1)+armthrsdist),:); % Find curr arm by using current well X-position
                    sstimarm{i}=stimarm;
                    
                case 1
                    stimarm = pos_stim(find((pos_stim(:,1)<currwellpos(1)+armthrsdist) & (pos_stim(:,1)>currwellpos(1)-armthrsdist)),:);
                    sstimarm{i}=stimarm;
                    
                case 3-fix
                    stimarm = pos_stim(find(pos_stim(:,1)>currwellpos(1)-armthrsdist),:);
                    sstimarm{i}=stimarm;
            end
            
            nstimarm(d,ep,i) = length(find( (dist(stimarm,repmat(currintpos,length(stimarm),1))>thrsdistint) & (dist(stimarm,repmat(currwellpos,length(stimarm),1))>thrsdistwell) ));
            pstimarm{i} = stimarm(find( (dist(stimarm,repmat(currintpos,length(stimarm),1))>thrsdistint) & (dist(stimarm,repmat(currwellpos,length(stimarm),1))>thrsdistwell) ),:);
            
        end
        edgewell=[2,3]; centerwell=1;
        nstimedgewell(d,ep,:) = nstimwell(d,ep,edgewell);
        nstimcenterwell(d,ep,1) = nstimwell(d,ep,centerwell);
        nstimcenterarm(d,ep,1) = nstimarm(d,ep,centerwell);
        
        % Figure for position at stimulation for current Day and Epoch if asked for
        % ---------------------------------------------
        if figopt1==1
            figure(day); hold on;
            %redimscreen_2x2subplots
            % Plot stimulation spreads
            subplot(2,length(allepochs),ep); hold on;
            plot(pos_stim(:,1),pos_stim(:,2),'.')
            % plot Well and Intersection positions, and circles around it
            for i=1:3
                plot(wellpos(i,1),wellpos(i,2),[wellclr{i} 's'],'MarkerSize',12,'MarkerFaceColor',[wellclr{i}]);
                plot(intpos(i,1),intpos(i,2),[wellclr{i} 's'],'MarkerSize',12,'MarkerFaceColor',[wellclr{i}]);
                plot(pstimarm{i}(:,1),pstimarm{i}(:,2),[wellclr{i} '.']);
            end
            
            %plot(wellpos(:,1),wellpos(:,2),'rs','MarkerSize',12,'MarkerFaceColor','r');
            %plot(intpos(:,1),intpos(:,2),'cs','MarkerSize',12,'MarkerFaceColor','c');
            
            % Plot circle
            for i=1:3
                N=256;
                circle(wellpos(i,:),thrsdistwell,256,'r-');
                circle(intpos(i,:),thrsdistint,256,'c-');
            end
            %axis([-5 150 -5 150]);
            xlabel('X-position (mm)');
            ylabel('Y-position (mm)');
            title(['Posn at Stim: Day' num2str(day) ' - Ep' num2str(ep) '; Nstim: ' num2str(length(pos_stim)) '; Rate: ' num2str(roundn(stimrate(ratecnt))) ' Hz'],'FontSize',20,'Fontweight','normal');
            
            % Plot stimulation statistics
            subplot(2,length(allepochs),length(allepochs)+ep); hold on;
            one = [nstimwell(d,ep,1), nstimarm(d,ep,1), nstimint(d,ep,1)];
            two = [nstimwell(d,ep,2), nstimarm(d,ep,2), nstimint(d,ep,2) ];
            thr = [nstimwell(d,ep,3), nstimarm(d,ep,3), nstimint(d,ep,3)];
            Y = [two; one; thr];
            bar(Y,'group');
            title(['No. of stim']);
            if ep==1,
                legend('Well','Arm','Interstn','Location','NorthEastOutside');
            end
            xlabel(['Section']);
            set(gca,'xtick',[1 2 3],'xticklabel',{'2';'1';'3'});
        end
    
        % --------------------
        % Save for day-epoch
        % --------------------
        % Stim-time and rate. Can also get stim-isi from this
        dayep_stimtime{d}{ep} = pt; % stimtimes in current day and epoch
        dayep_totalstim(d,ep) = length(pt);
        dayep_totaltime(d,ep) = totaltime; % in sec
        dayep_stimrate(d,ep) = length(pt)./(totaltime); % in Hz
        % Vel and position at stimulation
        dayep_velstim{d}{ep} = vel_stim;
        dayep_posstim{d}{ep} = pos_stim;
        % Vel and pos in entire epoch
        dayep_vel{d}{ep} = smooth_currvel;
        dayep_pos{d}{ep} = currpos;
        % Well and intersection positions for current epoch - if you need to calculate pos_stim
        dayep_wellpos{d}{ep} = wellpos;
        dayep_intpos{d}{ep} = intpos;
        % Old - Stuff you have calculated above about velocity
        dayep_velstim_moving{d}{ep} = moving_stimvel;
        dayep_velstim_still{d}{ep} = still_stimvel;
        dayep_time_moving(d,ep) = time_moving; % in sec
        dayep_time_still(d,ep) = time_still; % in sec
        dayep_vel_moving{d}{ep} = moving_vel;
        dayep_vel_still{d}{ep} = still_vel;
        
        % Reset values
        vel_stim=[]; % Most Imp. Have to do this.
        pos_stim=[]; currtime=[]; wellpos=[]; intpos=[]; sstimarm=[]; pstimarm=[];
        
        
    end  %% end epoch
    
    
    %  Plot velocity distribution at stimulation if asked for current day
    % -------------------------------------------------------------------
    if (figopt2==1)
        smooth_currvelday = [dayep_vel{d}{1}; dayep_vel{d}{2}];
        vel_stimday = [dayep_velstim{d}{1} dayep_velstim{d}{2}];
        
        figure(21); hold on;
        subplot(1,2,1); hold on;
        xhist = [0:1:50];
        vel_hist = hist(smooth_currvelday,xhist);
        norm_vel_hist = vel_hist ./ max(vel_hist);
        store_velhist(d,:) = norm_vel_hist;
        
        %             plot([0:1:50],hist(smooth_currvelday,[0:1:50]),[clr{d} '.-'],'MarkerSize',8,'LineWidth',2);
        %             plot(thrsvel*ones(size([0:10:max(hist(smooth_currvelday,[0:1:50]))])),[0:10:max(hist(smooth_currvelday,[0:1:50]))],'k--','Linewidth',1);
        %             plot(2*thrsvel*ones(size([0:10:max(hist(smooth_currvelday,[0:1:50]))])),[0:10:max(hist(smooth_currvelday,[0:1:50]))],'k--','Linewidth',1);
        plot(xhist,norm_vel_hist,[clr{d} '.-'],'MarkerSize',8,'LineWidth',2);
        plot(thrsvel*ones(size([0:0.1:1])),[0:0.1:1],'k--','Linewidth',1);
        plot(2*thrsvel*ones(size([0:0.1:1])),[0:0.1:1],'k--','Linewidth',1);
        if d==length(days)
            xlabel('Speed (cm/sec)');
            ylabel(['No']);
            title([prefix '- Speed Distr: Day ' num2str(min(days)) ' to ' num2str(max(days)) ],'FontSize',18,'Fontweight','normal');
        end
        
        subplot(1,2,2); hold on;
        %xhist = [0:1:floor(max(vel_stimday))];
        xhist = [0:1:50];
        yhist = hist(vel_stimday,xhist);
        norm_yhist = yhist./max(yhist);
        plot(xhist,norm_yhist,[clr{d} '.-'],'MarkerSize',16,'LineWidth',2);
        plot(thrsvel*ones(size([0:0.1:1])),[0:0.1:1],'k--','Linewidth',1);
        plot(2*thrsvel*ones(size([0:0.1:1])),[0:0.1:1],'k--','Linewidth',1);
        store_yhist(d,:) = norm_yhist;
        if d==length(days)
            xlabel('Speed (cm/sec)');
            ylabel(['No of Stimulations']);
            title([prefix '- Speed at Stim: Day ' num2str(min(days)) ' to ' num2str(max(days))],'FontSize',20,'Fontweight','normal');
        end
        % At end, plot Means of Distributions
        if d==length(days) && ep==2
            
            figure(24); hold on; redimscreen_figforppt1;            
            day_stimrate = mean((reshape(stimrate,2,length(days))),1);
            mean_stimrate = mean(day_stimrate);
            mean_store_yhist = mean(store_yhist,1);
            frac_below5 = sum(mean_store_yhist(1:7)) ./ sum(mean_store_yhist);
            frac_below5 = roundn(frac_below5,-3);
            per_below5 = frac_below5*100;
            
            frac_below10 = sum(mean_store_yhist(1:17)) ./ sum(mean_store_yhist);
            frac_below10 = roundn(frac_below10,-3);
            per_below10 = frac_below10*100;
            
            plot(xhist,mean(store_velhist,1),'b','MarkerSize',8,'LineWidth',3);
            plot(xhist,mean(store_yhist,1),'r','MarkerSize',8,'LineWidth',3);
            plot(thrsvel*ones(size([0:0.1:1])),[0:0.1:1],'k--','Linewidth',1);
            plot(2*thrsvel*ones(size([0:0.1:1])),[0:0.1:1],'k--','Linewidth',1);
            xlabel('Speed (cm/sec)');
            ylabel(['No']);
            text(20, 0.9,['(Avg Stim Rate = ' num2str(roundn(mean_stimrate)) ' Hz)'],...
                'FontSize',16,'FontWeight','normal');
            text(20, 0.7,['NStim below ' num2str(thrsvel) ' cm/sec = ' num2str(per_below5) '%'],'FontSize',18,'Fontweight','normal');
            text(20, 0.55,['NStim below ' num2str(2*thrsvel) ' cm/sec = ' num2str(per_below10) '%'],'FontSize',18,'Fontweight','normal');
            title([prefix '- Speed Distributions: Day ' num2str(min(days)) ' to ' num2str(max(days))],'FontSize',20,'Fontweight','normal');
        end
    end

    % Saving figures for current day
    % ------------------------------
    if saveg1==1,
        figure(day); hold on;
        figfile = [figdir,prefix,'_PosnAtStimStats_Day',num2str(day)];
        print('-dpdf', figfile);
        print('-djpeg', figfile);
        saveas(gcf,figfile,'fig');
    end
    
    % ---------------------------------
    % Save for Day - Avg across epochs
    % --------------------------------
    % Stimrate
    day_totalstim(d) = mean(dayep_totalstim(d,:),2);
    day_totaltime(d) = mean(dayep_totaltime(d,:),2); % in sec
    day_stimrate(d) = mean( dayep_stimrate(d,:),2); % in Hz
    % Vel and pos at stimln
    day_velstim{d} = [dayep_velstim{d}{1}';dayep_velstim{d}{2}']; % nX1 vector
    day_posstim{d} = [dayep_posstim{d}{1};dayep_posstim{d}{2}]; % nX2 matrix of x-y positions
    % Vel and pos in entire epoch
    day_vel{d} = [dayep_vel{d}{1}; dayep_vel{d}{2}]; % nX1 vector
    day_pos{d} = [dayep_pos{d}{1}; dayep_pos{d}{2}]; % nX2 matrix of x-y positions
    % Old - Stuff you have calculated above about velocity
    day_velstim_moving{d} = [dayep_velstim_moving{d}{1} dayep_velstim_moving{d}{2}];
    day_velstim_still{d} = [dayep_velstim_still{d}{1} dayep_velstim_still{d}{2}];
    day_time_moving(d) = mean(dayep_time_moving(d,:),2); % in sec
    day_time_still(d) = mean(dayep_time_still(d,:),2); % in sec
    day_vel_moving{d} = [dayep_vel_moving{d}{1}; dayep_vel_moving{d}{2}];
    day_vel_still{d} = [dayep_vel_still{d}{1}; dayep_vel_still{d}{2}];
    
    % Reminder - what is saved for each epoch
    % % Stim-time and rate. Can also get stim-isi from this
    % dayep_stimtime{d}{ep} = pt; % stimtimes in current day and epoch
    % dayep_totalstim(d,ep) = length(pt);
    % dayep_totaltime(d,ep) = totaltime; % in sec
    % dayep_stimrate(d,ep) = length(pt)./(dur*60); % in Hz
    % % Vel and position at stimulation
    % dayep_velstim{d}{ep} = vel_stim;
    % dayep_posstim{d}{ep} = pos_stim;
    % % Vel and pos in entire epoch
    % dayep_vel{d}{ep} = smooth_currvel;
    % dayep_pos{d}{ep} = currpos;
    % % Well and intersection positions for current epoch - if you need to calculate pos_stim
    % dayep_wellpos{d}{ep} = wellpos;
    % dayep_intpos{d}{ep} = intpos;
    %% Old - Stuff you have calculated above about velocity
    % dayep_velstim_moving{d}{ep} = moving_stimvel;
    % dayep_velstim_still{d}{ep} = still_stimvel;
    % dayep_time_moving(d,ep) = time_moving; % in sec
    % dayep_time_still(d,ep) = time_still; % in sec
    % dayep_vel_moving{d}{ep} = moving_vel;
    % dayep_vel_still{d}{ep} = still_vel;
    
end % end day



%%%%%%%% SAVING DATA %%%%%%%%%%%%%%%%%%%%
if savedata ==1
    savefile = sprintf('%s/%s_stimstats.mat', animdirect, prefix);
    %savefile = sprintf('%s/ProcessedData/%s_StimlnBehStats.mat', animdirect, prefix);
    save(savefile,'day_totalstim','day_totaltime','day_stimrate','day_velstim','day_velstim_moving','day_velstim_still',...
        'day_vel','day_vel_moving','day_vel_still','day_time_moving','day_time_still',...
        'dayep_stimtime','dayep_totalstim','dayep_totaltime','dayep_stimrate','dayep_velstim','dayep_posstim',...
        'dayep_vel','dayep_pos','dayep_wellpos','dayep_intpos');
    save(savefile,'day_*','dayep_*');
end

% ----------------
% Save Vel figure
% ---------------
if saveg2==1,
    figure(21); hold on;
    figfile = [figdir,prefix,'_VelAtStimStats_AcrossDays'];
    print('-dpdf', figfile);
    print('-djpeg', figfile);
    saveas(gcf,figfile,'fig');
end

% -------------------------------
% Plot across day comparisons
% ---------------------------
if figopt3==1
    figure(101); hold on;
    redimscreen_2versubplots
    orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
    
    wellintvsarm = []; well1=[]; well2 = []; well3=[];
    for d = 1:length(days)
        
        day = days(d);
        
        % Well and intersection vs. arm across days
        currday  = [sum(sum(nstimwell(d,:,:),3)) + sum(sum(nstimint(d,:,:),3)), sum(sum(nstimarm(d,:,:),3))];
        wellintvsarm = [wellintvsarm; currday];
        
        % Wells across days
        well1 = [well1, sum(nstimwell(d,:,1))];
        well2 = [well2, sum(nstimwell(d,:,2))];
        well3 = [well3, sum(nstimwell(d,:,3))];
        
    end
    
    % Well and intersection vs. arm
    subplot(2,1,1); hold on;
    bar(wellintvsarm,'group');
    title([prefix ': Arm stimlns vs. Well+Interstn stimlns']);
    legend('All Wells+Interstns','All Arms','Location','NorthEastOutside');
    xlabel(['Days']); ylabel(['No of stimulations']);
    %set(gca,'xtick',[1 2 3],'xticklabel',{'2';'1';'3'});
    %set(gca,'xtick',[days],'xticklabel',{num2str(days)'});
    
    % Wells across days
    subplot(2,1,2); hold on;
    Y = [well2; well1; well3];
    bar(Y,'group');
    title([prefix ': No. of well stimulns across days']);
    %legend('All Wells+Interstns','All Arms','Location','NorthEastOutside');
    xlabel(['Wells']); ylabel(['No of stimulations']);
    set(gca,'xtick',[1 2 3],'xticklabel',{'2';'1';'3'});
end

if saveg2==1,
    figure(21); hold on;
    orient(gcf,'landscape'); hold on; set(gcf, 'PaperPositionMode', 'auto');
    saveas(gcf,[prefix '_StimlnBehStats_AcrossDays'],'fig');
    saveas(gcf,[prefix '_StimlnBehStats_AcrossDays'],'jpg');
end


