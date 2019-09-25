function [params] = sj_HPexpt_eventtrigspecgram_gen(prefix, day, epochs, tet, trigtype, riptet, dospeed, figopt)
% Shantanu - Nov 2012. Generic version of event trigerred spectrogram - For ripples
% Get and plot ripple trig spec - Zscore using output generated by sj_HPexpt_baselinespecgram

if nargin<1,
    keyboard
    error('Please enter Expt Prefix and Day No!');
end
if nargin<2,
    keyboard
    error('Please enter Day No!');
end
if nargin<3,
    epochs=1; %% Epochs
end
if nargin<4,
    tet=1; %
end
if trigtype<5,
    trigtype='rip'; %
end
if nargin<6,
    riptet=4; % Have to start using getripples to look at ripples across multiple tetrodes!
end
if nargin<7,
    dospeed=0; %
end
if nargin<8,
    figopt=0; %
end
saveg1=0;


do_random = 0; % To align to random triggers



% Define Chronux params
% -------------------------------------------
win = [0.5 0.5]; % Window around triggering events / ripples
               % Change as appropriate for trigger type below 
movingwin = [100 10]/1000; cwin = movingwin;
params.Fs = 1500;
params.err = [2 0.05];
params.fpass = [0 400];
%params.trialave = 0; % This is 0 by default
params.tapers = [3 5]; % Should I put this in or let it use default tapers

% Speed parameters
lowsp_thrs = 2; % cm/sec
highsp_thrs = 7; % cm/sec

% Sleep
if ismember(epochs,[1 3 5 7]),
    lowsp_thrs = 0.5; %cm/sec
    highsp_thrs = 2; % cm/sec
end


% SET DATA
% -------------------------------------------

switch prefix
    case 'HPa'
        rawdir = '/data25/sjadhav/HPExpt/HPa/';
        directoryname = '/data25/sjadhav/HPExpt/HPa_direct/';
end
dir2=directoryname;

if (day<10)
    daystring = ['0',num2str(day)];
else
    daystring = num2str(day);
end

if (tet<10)
    tetstring = ['0',num2str(tet)];
else
    tetstring = num2str(tet);
end


S_all = []; Sgnd_all = [];
Smean_all=[]; Sgndmean_all=[];

for ep=1:length(epochs)
    
    epoch = epochs(ep);
    
    switch trigtype
        case 'rip'
            
            
            % Get Ripple Times - HAVE TO SWITCH TO USING getripples FOR LOOKING ACROSS MULTIPLE TETRODES
            % -------------------------------------------------------------------------------------------
            % SHOULD SWITCH THIS TO ALLTET. If alltet does not exist, then use given tet/tets.
            ripfile = sprintf('%s/%sripples%02d.mat', directoryname, prefix, day);
            load(ripfile);
            rip_starttime=[]; rip_sizes=[];
            for i=1:length(riptet)
                currriptet=riptet(i);
                rip_starttime = [rip_starttime; ripples{day}{epoch}{currriptet}.starttime];   % in sec
                rip_sizes = [rip_sizes; ripples{day}{epoch}{currriptet}.maxthresh];   % in units of std dev
            end
            
            rem = find(rip_sizes<4);
            rip_starttime(rem) = [];
            rip_sizes(rem) = [];
            
            [rip_starttime,sortidx] = sort(rip_starttime);
            rip_sizes = rip_sizes(sortidx);
            % Define triggering events as the start of each ripple
            triggers = rip_starttime;
            
            % Implement speed criterion
            if ~isempty(dospeed)
                posfile = sprintf('%s/%spos%02d.mat', directoryname, prefix, day);
                load(posfile);
                absvel = abs(pos{day}{epoch}.data(:,5)); % Can also use field 9
                postime = pos{day}{epoch}.data(:,1); % in secs
                
                pidx = lookup(triggers,postime);
                speed_atrip = absvel(pidx);
                lowsp_idx = find(speed_atrip <= lowsp_thrs);
                highsp_idx = find(speed_atrip >= highsp_thrs);
                
                if strcmp(dospeed,'low')==1
                    triggers = triggers(lowsp_idx);
                end
                
                if strcmp(dospeed,'high')==1
                    triggers = triggers(highsp_idx);
                end
            end
     
    end % end switch
    
    
    
    % Get the Baseline Spectrogram values
    % -----------------------------------
    cd([directoryname,'/EEGSpec/']);
    eegspecfile = [dir2,'/EEGSpec/',prefix,'eegspec',daystring,'-Tet',tetstring];
    load(eegspecfile);
    % Data
    %spec = eegspec{day}{epoch}{tet}.specgram;
    meanspec = eegspec{day}{epoch}{tet}.meanspec;
    stdspec = eegspec{day}{epoch}{tet}.stdspec;
    % Mean and std for whole day
    %lastepoch = length(eegspec{day});
    meandayspec = eegspec{day}{1}{tet}.meandayspec; % Stored in first epoch
    stddayspec = eegspec{day}{1}{tet}.stddayspec;
    %clear spec
    
    % Now get EEG for the given tet and process
    % --------------------------------------
    
    cd([directoryname,'/EEG/']);
    curreegfile = [dir2,'/EEG/',prefix,'eeg',daystring,'-',num2str(epoch),'-',tetstring];
    load(curreegfile);
    lfp = eeg{day}{epoch}{tet}.data;
    starttime = eeg{day}{epoch}{tet}.starttime; % This is in secs
    endtime = starttime + (length(lfp)-1) * (1 / params.Fs);
    clear eeg
    
    % Update triggers
    % -----------------
    % Subtract startime of epoch to get this with start at 0.
    % You only input eeg vector to mtspecgram, which starts at index 1
    
    triggers = triggers-starttime; endtime = endtime - starttime;
    if do_random ==1, triggers_rand = triggers_rand-starttime; end
    
    %Remove triggering events that are too close to the beginning or end
    while triggers(1)<win(1)
        triggers(1) = [];
    end
    rem = find(triggers + win(2) > endtime);
    triggers(rem) = [];
    %while triggers(end)> endtime-win(2)
    %    triggers(end) = [];
    %end
    %triggers(end)=[];
    
    % Calculate event triggered spectrogram
    % -----------------------------------
    disp(['Doing event-triggered specgram. Ntriggers = ',num2str(length(triggers))]);
    [S,Stime,Sfreq] = mtspecgramtrigc(lfp,triggers,[win(1) win(2)],[cwin(1) cwin(2)],params);
    Stime = Stime - win(1); % This will make time start -win(1) instead of 0
    % Alternatively
    %--------------
    % Can cut eeg in windows around triggers and build up a matrix, and then do specgram for
    % each cut piece separately. See event_spectrogram.m from kkay for an example
    clear lfp
    
    % Z-score the event-triggered  spectrogram using data from continuous spectrogram
    % ------------------------------------------------------------------------------
    S = bsxfun(@minus,S,meanspec); % Can use meandayspec and stddayspec instead
    S = bsxfun(@rdivide,S,stdspec);
    Smean = mean(S,3); % mean across events
    %S_all = [S_all;S];
    Smean_all(:,:,ep) = Smean;
    if do_random==1
        S_rand = bsxfun(@minus,S_rand,meanspec); S_rand = bsxfun(@rdivide,S_rand,stdspec);
        Smean_rand = mean(S_rand,3);
        %S_all = [S_all;S];
        Smean_all_rand(:,:,ep) = Smean_rand;
    end
    
    % Alternative way
    % for i = 1:size(S,1)
    %     for j = 1:size(S,3)
    %     	S(i,:,j) = (S(i,:,j) - meanspec)./stdspec;
    %     end
    % end
    
end % end epochs

% If multiple epochs, take mean across epochs
if length(epochs) > 1
    Smean = mean(Smean_all,3); % Mean across epochs
end

%% ------------------------------------------------
% PLOT
% -------------------------------------------------

if figopt ==1
    
    % ------------------------------
    % Figure and Font Sizes
    forppr = 1;
    % If yes, everything set to redimscreen_figforppr1
    % If not, everything set to redimscreen_figforppt1
    
    figdir = '/data25/sjadhav/';
    %figdir = '/data25/sjadhav/RippleInterruption/Figures/01AugSep11_RippleDisFigs/DisruptnCalibrationAndEgs/CalibrationEgs/';
    datadir = '/data25/sjadhav/';
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
    % ---------------------------------------
    
    
    
    %% Specgram EEG
    % -------------
    
    figure; hold on;
    set(gcf,'Position',[55 660 560 420]);
    imagesc(Stime,Sfreq,Smean'); colorbar;
    title(['Tet ',num2str(tet),' Specgram aligned to ',trigtype],'FontSize',18,'Fontweight','normal');
    ylabel('Freq','FontSize',20,'Fontweight','normal');
    xlabel('Time(s)','FontSize',20,'Fontweight','normal');
    set(gca,'XLim',[min(Stime) max(Stime)]);
    set(gca,'YLim',[min(Sfreq) max(Sfreq)]);
    
    % Plot Line at 0 ms - Start of ripple
    ypts = Sfreq;
    xpts = 0*ones(size(ypts));
    plot(xpts , ypts, 'k--','Linewidth',3);
    
    % Plot lines at 100ms
    xpts = (100)*ones(size(ypts));
    plot(xpts , ypts, 'k--','Linewidth',2);
    
    if saveg1==1,
        figfile = [figdir,prefix,'_Day',num2str(day),'_Tet',num2str(tetu),'_HpCellMatrix'];
        print('-dpdf', figfile);
        print('-djpeg', figfile);
        saveas(gcf,figfile,'fig');
    end
    
end % end figopt

i=1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


