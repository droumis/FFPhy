function [basefr0] = sj_HPexpt_rewardalign(prefix, day, epoch, tet, cell, binsize, saveg1)
% Shantanu - Nov 2012. From  sj_HPexpt_burstaligntorip
% Align spiking to reward - esp. for PFC
% sj_HPexpt_rewardalign('HPb',1,[4 6],9,3,100,0);

% sj_HPexpt_rewardalign('HPa',2,[2 4],15,1,100,0);
% sj_HPexpt_rewardalign('HPa',8,[2],17,1,100,0);

%---------------------------------------------------------
%[ Shantanu - May 2012. Also used to align to optical stimulation for animal JW7.
% From sj_multiunit_align, which is derived from sj_rippower_multiunit_calib.
% Change to getting times directly from epoch file. From given tetrode, plot mulitunit align, and spike align
% for either given cell or all cells on tetrode.
% Also, adding ability to use MU cluster rather than the multi file.] - Old

% eg. sj_multiunitandspike_align('JW7', 5, 1, 8, 1, [], 200, 0);


% Only aligning multiunit activity to DIO times.
% Enter times manually. Implement loading times file later

% From sj_rippower_stimampl3, which split into sj_rippower_multiunit_run and sj_rippower_multiunit_calib
% Only for calibration now, unlike sj_rippower_stimampl3  - needs manual input of time ranges



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




% % Load dio file
% %------------------
% DIOfile = sprintf('%s/%sDIO%02d.mat', directoryname, prefix, day);
% load(DIOfile);
% stim = DIO{day}{epoch}{16};
% stim_starttime = stim.pulsetimes(:,1)./10; %ms
% stim_endtime = stim.pulsetimes(:,2)./10; %ms
% stim_length = stim.pulselength;
% stim_isi = stim.timesincelast(2:end)./10; %ms

% Load reward file and Spikes File
% ---------------------------------
rewfile = sprintf('%s/%srewardinfo%02d.mat', directoryname, prefix, day);
load(rewfile);

spikefile = sprintf('%s/%sspikes%02d.mat', directoryname, prefix, day);
if exist(spikefile)==2
    load(spikefile);
end


% For epoch
cntallrew=0; cntallnorew=0; cntctrrew=0; cntsiderew=0; cntsidenorew=0;
cntside2rew=0; cntside2norew=0; cntside0rew=0; cntside0norew=0;
for e=1:length(epoch)
    ep=epoch(e);
    % Now Getting Range directly from times file - Not Used
    userange = ranges(ep+1,:); % 1st row is allepochs. So you need +1
    usename = names{ep+1}(end-15:end);

    currrewinfo = rewardinfo{day}{ep};
    %Convert rewardtime to ms
    currrewinfo(:,2) = currrewinfo(:,2)./10; % ms % This is from inputs
    currrewinfo(:,5) = currrewinfo(:,5)./10; % ms % This is with reward times updated from outputs 
    useidx=5; % Which one to use
    
    % 1)
    allwelltime = currrewinfo(:,useidx);
    % 2)
    allrewidx = find(currrewinfo(:,3)==1); allrewtime = currrewinfo(allrewidx,useidx);
    % 3)
    allnorewidx = find(currrewinfo(:,3)==0); allnorewtime = currrewinfo(allnorewidx,useidx);
    % 4)
    ctridx = find(currrewinfo(:,1)==1); ctrwelltime = currrewinfo(ctridx,useidx);
    ctrrewtime = ctrwelltime; % Always rewarded
    % 5)
    sideidx = find(currrewinfo(:,1)~=1); sidewelltime = currrewinfo(sideidx,useidx);
    sidewell_logic = currrewinfo(sideidx,3);
    siderewtime = sidewelltime(find(sidewell_logic==1));
    sidenorewtime = sidewelltime(find(sidewell_logic==0));
    % 6)
    side2idx = find(currrewinfo(:,1)==2); side2welltime = currrewinfo(side2idx,useidx);
    side2well_logic = currrewinfo(side2idx,3);
    side2rewtime = side2welltime(find(side2well_logic==1));
    side2norewtime = side2welltime(find(side2well_logic==0));
    % 7)
    side0idx = find(currrewinfo(:,1)==0); side0welltime = currrewinfo(side0idx,useidx);
    side0well_logic = currrewinfo(side0idx,3);
    side0rewtime = side0welltime(find(side0well_logic==1));
    side0norewtime = side0welltime(find(side0well_logic==0));
    
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
    
    % All UnRewarded
    for i=1:length(allnorewtime)
        cntallnorew=cntallnorew+1;
        currtime = (allnorewtime(i));
        currspks = spikeu(find( (spikeu>=(currtime-pret)) & (spikeu<=(currtime+postt)) ));
        %currspks = currspks-(currtime-pret); % Set to 0 at -pret
        %histspks = histc(currspks,[0:binsize:pret+postt]); % Make histogram
        currspks = currspks-(currtime); % Set to 0 at reward
        histspks = histc(currspks,[-pret:binsize:postt]); % Make histogram
        histspks = smoothvect(histspks, g1); % Smmoth histogram
        % Save spktimes
        allnorew_spks{cntallnorew}=currspks;
        % Save histogram
        allnorew_spkshist(cntallnorew,:)=histspks;
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
    
    % Side Well - NOT rewarded
    for i=1:length(sidenorewtime)
        cntsidenorew=cntsidenorew+1;
        currtime = (sidenorewtime(i));
        currspks = spikeu(find( (spikeu>=(currtime-pret)) & (spikeu<=(currtime+postt)) ));
        %currspks = currspks-(currtime-pret); % Set to 0 at -pret
        %histspks = histc(currspks,[0:binsize:pret+postt]); % Make histogram
        currspks = currspks-(currtime); % Set to 0 at reward
        histspks = histc(currspks,[-pret:binsize:postt]); % Make histogram
        histspks = smoothvect(histspks, g1); % Smmoth histogram
        % Save spktimes
        sidenorew_spks{cntsidenorew}=currspks;
        % Save histogram
        sidenorew_spkshist(cntsidenorew,:)=histspks;
    end
    
     % Side2 Well - rewarded
    for i=1:length(side2rewtime)
        cntside2rew=cntside2rew+1;
        currtime = (side2rewtime(i));
        currspks = spikeu(find( (spikeu>=(currtime-pret)) & (spikeu<=(currtime+postt)) ));
        %currspks = currspks-(currtime-pret); % Set to 0 at -pret
        %histspks = histc(currspks,[0:binsize:pret+postt]); % Make histogram
        currspks = currspks-(currtime); % Set to 0 at reward
        histspks = histc(currspks,[-pret:binsize:postt]); % Make histogram
        histspks = smoothvect(histspks, g1); % Smmoth histogram
        % Save spktimes
        side2rew_spks{cntside2rew}=currspks;
        % Save histogram
        side2rew_spkshist(cntside2rew,:)=histspks;
    end
    
    % Side2 Well - NOT rewarded
    for i=1:length(side2norewtime)
        cntside2norew=cntside2norew+1;
        currtime = (side2norewtime(i));
        currspks = spikeu(find( (spikeu>=(currtime-pret)) & (spikeu<=(currtime+postt)) ));
        %currspks = currspks-(currtime-pret); % Set to 0 at -pret
        %histspks = histc(currspks,[0:binsize:pret+postt]); % Make histogram
        currspks = currspks-(currtime); % Set to 0 at reward
        histspks = histc(currspks,[-pret:binsize:postt]); % Make histogram
        histspks = smoothvect(histspks, g1); % Smmoth histogram
        % Save spktimes
        side2norew_spks{cntside2norew}=currspks;
        % Save histogram
        side2norew_spkshist(cntside2norew,:)=histspks;
    end
    
      % Side0 Well - rewarded
    for i=1:length(side0rewtime)
        cntside0rew=cntside0rew+1;
        currtime = (side0rewtime(i));
        currspks = spikeu(find( (spikeu>=(currtime-pret)) & (spikeu<=(currtime+postt)) ));
        %currspks = currspks-(currtime-pret); % Set to 0 at -pret
        %histspks = histc(currspks,[0:binsize:pret+postt]); % Make histogram
        currspks = currspks-(currtime); % Set to 0 at reward
        histspks = histc(currspks,[-pret:binsize:postt]); % Make histogram
        histspks = smoothvect(histspks, g1); % Smmoth histogram
        % Save spktimes
        side0rew_spks{cntside0rew}=currspks;
        % Save histogram
        side0rew_spkshist(cntside0rew,:)=histspks;
    end
    
    % Side0 Well - NOT rewarded
    for i=1:length(side0norewtime)
        cntside0norew=cntside0norew+1;
        currtime = (side0norewtime(i));
        currspks = spikeu(find( (spikeu>=(currtime-pret)) & (spikeu<=(currtime+postt)) ));
        %currspks = currspks-(currtime-pret); % Set to 0 at -pret
        %histspks = histc(currspks,[0:binsize:pret+postt]); % Make histogram
        currspks = currspks-(currtime); % Set to 0 at reward
        histspks = histc(currspks,[-pret:binsize:postt]); % Make histogram
        histspks = smoothvect(histspks, g1); % Smmoth histogram
        % Save spktimes
        side0norew_spks{cntside0norew}=currspks;
        % Save histogram
        side0norew_spkshist(cntside0norew,:)=histspks;
    end
    
    
end % End Epoch


%% ------------------------------------------------
% PLOT
% -------------------------------------------------

% All Rewarded vs Not Rewarded
% ----------------------------

figure; hold on;
%redimscreen_figforppt1;
%set(gcf,'Position',[0 Screen(4)*0.55 Screen(3)*0.2 Screen(4)*0.4])
%redimscreen_2versubplots
set(gcf,'Position',[16 44 900 1045]);

% Raster
% ---------
subplot(3,1,[1 2]); hold on; % Raster - rewarded and then unrewarded
for c=1:length(allrew_spks)
    tmps = allrew_spks{c};
    plotraster(tmps,c*ones(size(tmps)),0.8,[],'LineWidth',2,'Color','b');
end
push = length(allrew_spks);
for c=1:length(allnorew_spks)
    tmps = allnorew_spks{c};
    plotraster(tmps,(push+c)*ones(size(tmps)),0.8,[],'LineWidth',2,'Color','r');
end
    
set(gca,'XLim',[-pret+300 postt-300]);
set(gca,'YLim',[0 push+length(allnorew_spks)+1]);
set(gca,'XTick',[]); set(gca,'YTick',[]);
ypts = 0:1:push+length(allnorew_spks)+1;
xpts = 0*ones(size(ypts));
plot(xpts , ypts, 'k--','Linewidth',2);


subplot(3,1,3); hold on; % Histogram - rewarded
fr = mean(allrew_spkshist)*(1000/binsize); % It is mean - so normalized to Ntrials
xaxis = -pret:binsize:postt;
plot(xaxis,fr,'b','Linewidth',3); % Convert to firing rate in Hz
fr2 = mean(allnorew_spkshist)*(1000/binsize);
plot(xaxis,fr2,'r','Linewidth',3); % Convert to firing rate in Hz
xlabel('Time(ms)','FontSize',20,'Fontweight','normal');
ylabel('Firing Rate (Hz)','FontSize',16,'Fontweight','normal');
%set(gca,'xtick',[(-pret/binsize):(1000/binsize):(postt/binsize)],'xticklabel',{num2str([-pret:1000:postt]')},...
    %'FontSize',16,'Fontweight','normal');
legend('Rew','NoRew','Location','NorthWest');
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
    figfile = [figdir,prefix,'_Day',num2str(day),'_Tet',num2str(tetu),'_Rew'];
    print('-dpdf', figfile);
    print('-djpeg', figfile);
    saveas(gcf,figfile,'fig');
end

% ----------------------------------------------------------------------------

% Ctr Rewarded, Side Rewarded vs Side Not Rewarded
% -------------------------------------------------

figure; hold on;
set(gcf,'Position',[984 42 900 1045]);

subplot(3,1,[1 2]); hold on; % 
for c=1:length(ctrrew_spks) % Ctr Rew Spks
    tmps = ctrrew_spks{c};
    plotraster(tmps,c*ones(size(tmps)),0.8,[],'LineWidth',2,'Color','b');
end
push = length(ctrrew_spks);
for c=1:length(siderew_spks) % Side Rew spks
    tmps = siderew_spks{c};
    plotraster(tmps,(push+c)*ones(size(tmps)),0.8,[],'LineWidth',2,'Color','c');
end
push = push + length(siderew_spks);
for c=1:length(sidenorew_spks) % Side NoRew spks
    tmps = sidenorew_spks{c};
    plotraster(tmps,(push+c)*ones(size(tmps)),0.8,[],'LineWidth',2,'Color','r');
end
set(gca,'XLim',[-pret+300 postt-300]);
set(gca,'YLim',[0 push+length(sidenorew_spks)+1]);
set(gca,'XTick',[]); set(gca,'YTick',[]);
ypts = 0:1:push+length(allnorew_spks)+1;
xpts = 0*ones(size(ypts));
plot(xpts , ypts, 'k--','Linewidth',2);

subplot(3,1,3); hold on; % Histogram - rewarded
fr = mean(ctrrew_spkshist)*(1000/binsize);
xaxis = -pret:binsize:postt;
plot(xaxis,fr,'b','Linewidth',3); % Convert to firing rate in Hz
fr2 = mean(siderew_spkshist)*(1000/binsize);
plot(xaxis,fr2,'c','Linewidth',3); % Convert to firing rate in Hz
fr3 = mean(sidenorew_spkshist)*(1000/binsize);
plot(xaxis,fr3,'r','Linewidth',3); % Convert to firing rate in Hz
xlabel('Time(ms)','FontSize',20,'Fontweight','normal');
ylabel('Firing Rate (Hz)','FontSize',16,'Fontweight','normal');
legend('CtrRew','SideRew','SideNoRew','Location','NorthWest');
set(gca,'XLim',[-pret+300 postt-300]);
set(gca,'YLim',[min([fr, fr2, fr3])-0.5 max([fr, fr2, fr3])+0.5]);
ypts = min([fr, fr2, fr3])-0.5:0.1:max([fr, fr2, fr3])+0.5;
xpts = 0*ones(size(ypts));
plot(xpts , ypts, 'k--','Linewidth',2);


if saveg1==1,
    figfile = [figdir,prefix,'_Day',num2str(day),'_Tet',num2str(tetu),'_Rew'];
    print('-dpdf', figfile);
    print('-djpeg', figfile);
    saveas(gcf,figfile,'fig');
end

% ----------------------------------------------------------------------------

% % Ctr Rewarded, Side2 Rewarded vs Side0 Rewarded
% % ----------------------------

figure; hold on;
set(gcf,'Position',[984 42 900 1045]);

subplot(3,1,[1 2]); hold on; % 
for c=1:length(ctrrew_spks) % Ctr Rew Spks
    tmps = ctrrew_spks{c};
    plotraster(tmps,c*ones(size(tmps)),0.8,[],'LineWidth',2,'Color','b');
end
push = length(ctrrew_spks);
for c=1:length(side2rew_spks) % Side2 Rew spks
    tmps = side2rew_spks{c};
    plotraster(tmps,(push+c)*ones(size(tmps)),0.8,[],'LineWidth',2,'Color','c');
end
push = push + length(side2rew_spks);
for c=1:length(side0rew_spks) % Side0 Rew spks
    tmps = side0rew_spks{c};
    plotraster(tmps,(push+c)*ones(size(tmps)),0.8,[],'LineWidth',2,'Color','g');
end


set(gca,'XLim',[-pret+300 postt-300]);
set(gca,'YLim',[0 push+length(side0rew_spks)+1]);
set(gca,'XTick',[]); set(gca,'YTick',[]);
ypts = 0:1:push+length(allnorew_spks)+1;
xpts = 0*ones(size(ypts));
plot(xpts , ypts, 'k--','Linewidth',2);

subplot(3,1,3); hold on; % Histogram - rewarded
fr = mean(ctrrew_spkshist)*(1000/binsize);
xaxis = -pret:binsize:postt;
plot(xaxis,fr,'b','Linewidth',3);
fr2 = mean(side2rew_spkshist)*(1000/binsize);
plot(xaxis,fr2,'c','Linewidth',3);
fr3 = mean(side0rew_spkshist)*(1000/binsize);
plot(xaxis,fr3,'g','Linewidth',3); 

xlabel('Time(ms)','FontSize',20,'Fontweight','normal');
ylabel('Firing Rate (Hz)','FontSize',16,'Fontweight','normal');
legend('CtrRew','Side2Rew','Side0Rew','Location','NorthWest');
set(gca,'XLim',[-pret+300 postt-300]);
set(gca,'YLim',[min([fr, fr2, fr3])-0.5 max([fr, fr2, fr3])+0.5]);
ypts = min([fr, fr2, fr3])-0.5:0.1:max([fr, fr2, fr3])+0.5;
xpts = 0*ones(size(ypts));
plot(xpts , ypts, 'k--','Linewidth',2);


if saveg1==1,
    figfile = [figdir,prefix,'_Day',num2str(day),'_Tet',num2str(tetu),'_Rew'];
    print('-dpdf', figfile);
    print('-djpeg', figfile);
    saveas(gcf,figfile,'fig');
end





% ----------------------------------------------------------------------------

% % Ctr Rewarded, Side2/3 Rewarded vs Side2/3 Not Rewarded
% % ----------------------------
% 
% figure; hold on;
% set(gcf,'Position',[984 42 900 1045]);
% 
% subplot(3,1,[1 2]); hold on; % 
% for c=1:length(ctrrew_spks) % Ctr Rew Spks
%     tmps = ctrrew_spks{c};
%     plotraster(tmps,c*ones(size(tmps)),0.8,[],'LineWidth',2,'Color','b');
% end
% push = length(ctrrew_spks);
% for c=1:length(side2rew_spks) % Side2 Rew spks
%     tmps = side2rew_spks{c};
%     plotraster(tmps,(push+c)*ones(size(tmps)),0.8,[],'LineWidth',2,'Color','c');
% end
% push = push + length(side2rew_spks);
% for c=1:length(side2norew_spks) % Side2 NoRew spks
%     tmps = side2norew_spks{c};
%     plotraster(tmps,(push+c)*ones(size(tmps)),0.8,[],'LineWidth',2,'Color','r');
% end
% push = push + length(side2norew_spks);
% for c=1:length(side0rew_spks) % Side0 Rew spks
%     tmps = side0rew_spks{c};
%     plotraster(tmps,(push+c)*ones(size(tmps)),0.8,[],'LineWidth',2,'Color','g');
% end
% push = push + length(side0rew_spks);
% for c=1:length(side0norew_spks) % Side0 NoRew spks
%     tmps = side0norew_spks{c};
%     plotraster(tmps,(push+c)*ones(size(tmps)),0.8,[],'LineWidth',2,'Color','m');
% end
% 
% 
% set(gca,'XLim',[-pret+300 postt-300]);
% set(gca,'YLim',[0 push+length(side0norew_spks)+1]);
% set(gca,'XTick',[]); set(gca,'YTick',[]);
% ypts = 0:1:push+length(allnorew_spks)+1;
% xpts = 0*ones(size(ypts));
% plot(xpts , ypts, 'k--','Linewidth',2);
% 
% subplot(3,1,3); hold on; % Histogram - rewarded
% fr = mean(ctrrew_spkshist)*(1000/binsize);
% xaxis = -pret:binsize:postt;
% plot(xaxis,fr,'b','Linewidth',3);
% fr2 = mean(side2rew_spkshist)*(1000/binsize);
% plot(xaxis,fr2,'c','Linewidth',3);
% fr3 = mean(side2norew_spkshist)*(1000/binsize);
% plot(xaxis,fr3,'r','Linewidth',3);
% fr4 = mean(side0rew_spkshist)*(1000/binsize);
% plot(xaxis,fr4,'g','Linewidth',3); 
% fr5 = mean(side0norew_spkshist)*(1000/binsize);
% plot(xaxis,fr5,'m','Linewidth',3); 
% 
% xlabel('Time(ms)','FontSize',20,'Fontweight','normal');
% ylabel('Firing Rate (Hz)','FontSize',16,'Fontweight','normal');
% legend('CtrRew','Side2Rew','Side2NoRew','Side0Rew','Side0NoRew','Location','NorthWest');
% set(gca,'XLim',[-pret+300 postt-300]);
% set(gca,'YLim',[min([fr, fr2, fr3, fr4, fr5])-0.5 max([fr, fr2, fr3, fr4, fr5])+0.5]);
% ypts = min([fr, fr2, fr3, fr4, fr5])-0.5:0.1:max([fr, fr2, fr3, fr4, fr5])+0.5;
% xpts = 0*ones(size(ypts));
% plot(xpts , ypts, 'k--','Linewidth',2);
% 
% 
% if saveg1==1,
%     figfile = [figdir,prefix,'_Day',num2str(day),'_Tet',num2str(tetu),'_Rew'];
%     print('-dpdf', figfile);
%     print('-djpeg', figfile);
%     saveas(gcf,figfile,'fig');
% end




i=1;
























%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



