function [basefr0] = sj_HPexpt_rewardalign_lintrack(prefix, day, epoch, tet, cell, binsize, saveg1)
% Shantanu - Nov 2012. From   sj_HPexpt_rewardalign
% Align spiking to reward - esp. for PFC
% sj_HPexpt_rewardalign_lintrack('HPb',1,[2],12,2,100,0);

% For Wtrack. it was:
% sj_HPexpt_rewardalign('HPb',1,[4 6],9,3,100,0);

% sj_HPexpt_rewardalign('HPa',2,[2 4],15,1,100,0);
% sj_HPexpt_rewardalign('HPa',8,[2],17,1,100,0);

%---------------------------------------------------------


if nargin<1,
    keyboard
    error('Please enter Expt Prefix and Day No!');
end
if nargin<2,
    keyboard
    error('Please enter Day No!');
end
if nargin<3,
    epoch=2; %% Epoch
end
if nargin<4,
    tet=1; %
end
if nargin<5,
    cell=1; %
end
if nargin<6,
    binsize=10;
end
if nargin<7
    saveg1=0;
end

doraster_mu=0;
plotline1 = 100; plotline2=[]; % Where to plot lines after reward
Screen = get(0,'Screensize');

% --------------- Parameters ---------------
pret=4000; postt=6000; %% Times to plot - For 2 sec
plotline1 = 100;
if isempty(binsize)
    binsize = 100;  %% ms
end
binsize_plot=binsize; %% If you put binsize_plot=1000, then units are Nspikes/binsize, not inst. firing rate in Hz
timeaxis = -pret:binsize:postt;
bins_center = find(abs(timeaxis)<=100);
bins_resp = find(timeaxis>=0 & timeaxis<=100);
smwin=100; %Smoothing Window
% Smooth over binsize*3 ms (30ms for 10ms; 15ms for 5 ms)
%nstd = round(binsize*3/binsize); g1 = gaussian(nstd, 5*nstd+1);
nstd = round(binsize*2/binsize); g1 = gaussian(nstd, 1*nstd+1);

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


%% -----------------------------------------
% SET DATA
% -------------------------------------------

switch prefix
    case 'HPa'
        rawdir = '/data25/sjadhav/HPExpt/HPa/';
        directoryname = '/data25/sjadhav/HPExpt/HPa_direct/';
    case 'HPb'
        rawdir = '/data25/sjadhav/HPExpt/HPb/';
        directoryname = '/data25/sjadhav/HPExpt/HPb_direct/';
end


clr = {'b','g','c','m','y','k','r'};

%% --------------------------------------------------
%  %Align Single Cell Firing Rate to Reward
% -------------------------------------------------

% Load times file - Not Used here
% -----------------
currdir = pwd;
cd(rawdir);
dayfolders = dir;
daystr = sprintf('%02d', day);
for i = 3:length(dayfolders)
    if dayfolders(i).isdir
        if strcmp(daystr,dayfolders(i).name(1:2))
            disp(upper(dayfolders(i).name))
            cd(dayfolders(i).name);
            load times;
        end
    end
end
cd(currdir);


% Load reward file and Spikes File
% ---------------------------------
rewfile = sprintf('%s/%srewardinfo%02d.mat', directoryname, prefix, day);
load(rewfile);

spikefile = sprintf('%s/%sspikes%02d.mat', directoryname, prefix, day);
if exist(spikefile)==2
    load(spikefile);
end


% For epoch
cntallrew=0; cntctrrew=0; cntsiderew=0; 
%cntside2rew=0; cntside2norew=0; cntside0rew=0; cntside0norew=0; cntallnorew=0; cntsidenorew=0;
for e=1:length(epoch)
    ep=epoch(e);
    % Now Getting Range directly from times file - Not Used
    userange = ranges(ep+1,:); % 1st row is allepochs. So you need +1
    usename = names{ep+1}(end-15:end);

    currrewinfo = rewardinfo{day}{ep};
    %Convert rewardtime to ms
    currrewinfo(:,2) = currrewinfo(:,2)./10; % ms % This is from inputs
    currrewinfo(:,5) = currrewinfo(:,5)./10; % ms % This is with reward times updated from outputs 
    useidx=5; % Which one to use. Only Matters for W-track
    
    % 1)
    allwelltime = currrewinfo(:,useidx);
    % 2) For lintrack, this is same as 1)
    allrewidx = find(currrewinfo(:,3)==1); allrewtime = currrewinfo(allrewidx,useidx); 
    % 3)
    ctridx = find(currrewinfo(:,1)==1); ctrwelltime = currrewinfo(ctridx,useidx);
    ctrrewtime = ctrwelltime; % Always rewarded
    % 4)
    sideidx = find(currrewinfo(:,1)~=1); sidewelltime = currrewinfo(sideidx,useidx);
    sidewell_logic = currrewinfo(sideidx,3); % This will be always 1 for now
    siderewtime = sidewelltime(find(sidewell_logic==1));
    
    
    % Do Spikes
    % -----------
    % Get spikes in ms - Multiply by 1000
    if ~isempty(spikes{day}{ep}{tet})
        if ~isempty(spikes{day}{ep}{tet}{cell})
            
            % Display tag
            currtag = spikes{day}{ep}{tet}{cell}.tag;
            disp(['Cell Tag is :' currtag]);
            
            if ~isempty(spikes{day}{ep}{tet}{cell}.data)
                cmd=sprintf('spikeu = spikes{day}{ep}{tet}{%d}.data(:,1)*1000;',cell); eval(cmd);
            else
                disp(['Error - No .data for tet-cell. Epoch: ',num2str(ep)]); keyboard;
            end
        else
            disp(['Error - No {cell} for this tet. Epoch: ',num2str(ep)]); keyboard;
        end
    else
        disp(['Error - No {tet} for this day-epoch. Epoch: ',num2str(ep)]); keyboard;
    end
    
    
    
    % ALINGING TO REW Time
    % ----------------------
    
    % All Rewarded
    
    for i=1:length(allrewtime)
        cntallrew=cntallrew+1;
        currtime = (allrewtime(i));
        currspks = spikeu(find( (spikeu>=(currtime-pret)) & (spikeu<=(currtime+postt)) ));
        %currspks = currspks-(currtime-pret); % Set to 0 at -pret
        %histspks = histc(currspks,[0:binsize:pret+postt]); % Make histogram
        currspks = currspks-(currtime); % Set to 0 at reward
        histspks = histc(currspks,[-pret:binsize:postt]); % Make histogram
        histspks = smoothvect(histspks, g1); % Smmoth histogram
        % Save spktimes
        allrew_spks{cntallrew}=currspks;
        % Save histogram
        allrew_spkshist(cntallrew,:)=histspks;
    end
    
    % Ctr Well - always rewarded
    for i=1:length(ctrrewtime)
        cntctrrew=cntctrrew+1;
        currtime = (ctrrewtime(i));
        currspks = spikeu(find( (spikeu>=(currtime-pret)) & (spikeu<=(currtime+postt)) ));
        %currspks = currspks-(currtime-pret); % Set to 0 at -pret
        %histspks = histc(currspks,[0:binsize:pret+postt]); % Make histogram
        currspks = currspks-(currtime); % Set to 0 at reward
        histspks = histc(currspks,[-pret:binsize:postt]); % Make histogram
        histspks = smoothvect(histspks, g1); % Smmoth histogram
        % Save spktimes
        ctrrew_spks{cntctrrew}=currspks;
        % Save histogram
        ctrrew_spkshist(cntctrrew,:)=histspks;
    end
    
    % Side Well - rewarded
    for i=1:length(siderewtime)
        cntsiderew=cntsiderew+1;
        currtime = (siderewtime(i));
        currspks = spikeu(find( (spikeu>=(currtime-pret)) & (spikeu<=(currtime+postt)) ));
        %currspks = currspks-(currtime-pret); % Set to 0 at -pret
        %histspks = histc(currspks,[0:binsize:pret+postt]); % Make histogram
        currspks = currspks-(currtime); % Set to 0 at reward
        histspks = histc(currspks,[-pret:binsize:postt]); % Make histogram
        histspks = smoothvect(histspks, g1); % Smmoth histogram
        % Save spktimes
        siderew_spks{cntsiderew}=currspks;
        % Save histogram
        siderew_spkshist(cntsiderew,:)=histspks;
    end
    
end % End Epoch


%% ------------------------------------------------
% PLOT
% -------------------------------------------------

% All Rewarded vs Well1 Rewarded vs Well0 Rewarded
% ---------------------------------------------------

figure; hold on;
%redimscreen_figforppt1;
%set(gcf,'Position',[0 Screen(4)*0.55 Screen(3)*0.2 Screen(4)*0.4])
%redimscreen_2versubplots
set(gcf,'Position',[16 44 900 1045]);

% Raster
% ---------
subplot(3,1,[1 2]); hold on; % Raster - rewarded and then unrewarded
% for c=1:length(allrew_spks)
%     tmps = allrew_spks{c};
%     plotraster(tmps,c*ones(size(tmps)),0.8,[],'LineWidth',2,'Color','b');
% end
% push = length(allrew_spks);
push=0;
for c=1:length(ctrrew_spks)
    tmps = ctrrew_spks{c};
    plotraster(tmps,(push+c)*ones(size(tmps)),0.8,[],'LineWidth',2,'Color','c');
end
push = push+length(ctrrew_spks);
for c=1:length(siderew_spks)
    tmps = siderew_spks{c};
    plotraster(tmps,(push+c)*ones(size(tmps)),0.8,[],'LineWidth',2,'Color','g');
end
    
set(gca,'XLim',[-pret+300 postt-300]);
set(gca,'YLim',[0 push+length(siderew_spks)+1]);
set(gca,'XTick',[]); set(gca,'YTick',[]);
ypts = 0:1:push+length(siderew_spks)+1;
xpts = 0*ones(size(ypts));
plot(xpts , ypts, 'k--','Linewidth',2);


subplot(3,1,3); hold on; % Histogram - rewarded
%fr = mean(allrew_spkshist)*(1000/binsize);
xaxis = -pret:binsize:postt;
%plot(xaxis,fr,'b','Linewidth',3); % Convert to firing rate in Hz
fr = mean(ctrrew_spkshist)*(1000/binsize);
plot(xaxis,fr,'c','Linewidth',3); % Convert to firing rate in Hz
fr2 = mean(siderew_spkshist)*(1000/binsize);
plot(xaxis,fr2,'g','Linewidth',3); 

xlabel('Time(ms)','FontSize',20,'Fontweight','normal');
ylabel('Firing Rate (Hz)','FontSize',16,'Fontweight','normal');
%set(gca,'xtick',[(-pret/binsize):(1000/binsize):(postt/binsize)],'xticklabel',{num2str([-pret:1000:postt]')},...
    %'FontSize',16,'Fontweight','normal');
legend('Well1Rew','Well0Rew','Location','NorthWest');
%set(gca,'XLim',[1 ((pret+postt)/binsize)+1]);
set(gca,'XLim',[-pret+300 postt-300]);
set(gca,'YLim',[min([fr, fr2]) max([fr, fr2])+0.5]);
ypts = min([fr, fr2]):0.1:max([fr, fr2])+0.5;
xpts = 0*ones(size(ypts));
% Plot Line at 0 ms - Onset of reward
plot(xpts , ypts, 'k--','Linewidth',2);

% Plot lines at plotline1 ms 
%xpts = ((pret+plotline1+binsize)/binsize)*ones(size(ypts));
%plot(xpts , ypts, 'k--','Linewidth',1);

if saveg1==1,
    figfile = [figdir,prefix,'_Day',num2str(day),'_Tet',num2str(tetu),'_LinRew'];
    print('-dpdf', figfile);
    print('-djpeg', figfile);
    saveas(gcf,figfile,'fig');
end

% ----------------------------------------------------------------------------





i=1;
























%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



