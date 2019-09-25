function [pret] = sj_ripegs_run3(prefix, day, epoch,nostim, saveg1)
% Align to Allripples - mainly for control. Get ripple time from allripples. 
% sj_ripegs_run2 aligns to ripples detected on individual tetrode,
% sj_ripegs_run which aligns to stimulation times 
% like sj_rippower_multiunit_run
% plot example LFPs 
%
% Autoload times file to get time ranges


if nargin<1,
    keyboard
    error('Please enter Expt Prefix and Day No!');
end
if nargin<2,
    keyboard
    error('Please enter Day No!');
end
if nargin<3,
    epoch=2; %% Epoch - 2 or 4 for runs
end
if nargin<4
    nostim=0; 
end
if nargin<5
    saveg1=0; % Save summary figure
end

% --------------- Parameters ---------------
%rem_mustimart=0; % Remove stim artifact in MU and rippower firing while plotting - for exp and con, not nor

iscon=1; % Whether control or exptal stimulation - for rippower alignment
sd=7; %% SD threshold for plotting 

if iscon==1
    pret=50;
    postt=400; %% Times to plot
else
    pret=100;
    postt=100;
end

% ------------------------------
% Figure and Font Sizes
forppr = 1;
% If yes, everything set to redimscreen_figforppr1
% If not, everything set to redimscreen_figforppt1

figdir = '/data25/sjadhav/RippleInterruption/Figures/01AugSep11_RippleDisFigs/DisruptnCalibrationAndEgs/DisruptionEgs/';
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

datadir = '/data25/sjadhav/RippleInterruption/ProcessedData/';
% ---------------------------------------


%% -----------------------------------------
% SET DATA
% -------------------------------------------

switch prefix
    case 'REb'
        directoryname = '/data25/sjadhav/RippleInterruption/REb_direct';
        dire = '/data25/sjadhav/RippleInterruption/REb';
    case 'RE1'
        directoryname = '/data25/sjadhav/RippleInterruption/RE1_direct';
        dire = '/data25/sjadhav/RippleInterruption/RE1';
        riptetlist=[1,5,7];
    case 'REc'
        directoryname = '/data25/sjadhav/RippleInterruption/REc_direct';
        dire = '/data25/sjadhav/RippleInterruption/REc';
        riptetlist = [4,6,8,9,10];
    case 'REd'
        directoryname = '/data25/sjadhav/RippleInterruption/REd_direct';
        dire = '/data25/sjadhav/RippleInterruption/REd';
        riptetlist = [3,4,5,6,10,11];
        riptetlist = [4,5,6,10];
    case 'REe'
        directoryname = '/data25/sjadhav/RippleInterruption/REe_direct';
        dire = '/data25/sjadhav/RippleInterruption/REe';
        %riptetlist = [3,4,6,11,12,13];
        riptetlist = [3,4,6,11,12,13];
        riptetlist = [3,4,6,11];
    case 'REf'
        directoryname = '/data25/sjadhav/RippleInterruption/REf_direct';
        dire = '/data25/sjadhav/RippleInterruption/REf';
        riptetlist = [1,5,9,10,11,12];
        riptetlist = [1,5,9,11];
    case 'RCa'
        directoryname = '/data25/sjadhav/RippleInterruption/RCa_direct';
        dire = '/data25/sjadhav/RippleInterruption/RCa';
        riptetlist = [2,3,4,6,9];
    case 'RCb'
        directoryname = '/data25/sjadhav/RippleInterruption/RCb_direct';
        dire = '/data25/sjadhav/RippleInterruption/RCb';
        %riptetlist = [3,4,9,10,11,12];
        riptetlist = [3,4,9,10];
    case 'RCc'
        directoryname = '/data25/sjadhav/RippleInterruption/RCc_direct';
        dire = '/data25/sjadhav/RippleInterruption/RCc';
        riptetlist = [3,4,5,6,11,13];
        riptetlist = [3,4,5,6];
    case 'RCd'
        directoryname = '/data25/sjadhav/RippleInterruption/RCd_direct';
        dire = '/data25/sjadhav/RippleInterruption/RCd';
        riptetlist = [1,2,3,4,5,6];
end

currdir = pwd;
if (directoryname(end) == '/')
    directoryname = directoryname(1:end-1);
end
if (dire(end) == '/')
    dire = dire(1:end-1);
end

if (day < 10)
    daystring = ['0',num2str(day)];
else
    daystring = num2str(day);
end

%% --------------------------------------------------
% Align MU Firing Rate to Ripples
% -------------------------------------------------

% Load dio file
%------------------
% DIOfile = sprintf('%s/%sDIO%02d.mat', directoryname, prefix, day);
% load(DIOfile);
% stim = DIO{day}{epoch}{16};
% % if isempty(stim)
% %             stim = DIO{day}{epoch}{15};
% %         end
% stim_starttime = stim.pulsetimes(:,1)./10; %ms
% stim_endtime = stim.pulsetimes(:,2)./10; %ms
% stim_length = stim.pulselength;
% stim_isi = stim.timesincelast(2:end)./10; %ms
% pt = stim.pulsetimes ./ 10000; % in s

% Load Allripple file
ripfile = sprintf('%s/%sripplesall%02d.mat', directoryname, prefix, day);
load(ripfile);

% Get ripstarttime or midtime for current tetrodes in given day and epoch
rip_starttime = ripplesall{day}{epoch}{2}.starttime;   % in sec - Ntet=2 condition
ripsize = ripplesall{day}{epoch}{2}.maxthresh;
pt = rip_starttime;
%eind = lookup(pt(:,1), ti);
% if iscon==1
%     eind = eind - round(150*e.samprate/1000); % Align to ripple det in Con
% end

tets=riptetlist;

for te=1:length(tets)
    
    currtet=tets(te);
    cnteeg=0; cntrip=0;
    
    % Load EEG and ripple LFP file
    %-------------------------
    if nostim==0 % stimulation artifact not removed  
        EEGfile = sprintf('%s/EEG/%seeg%02d-%01d-%02d.mat', directoryname, prefix, day,epoch,currtet);
        load(EEGfile);
        e = eeg{day}{epoch}{currtet};
        if te==1
            ti = geteegtimes(e);
            eind = lookup(pt, ti);
            e.samprate=round(e.samprate);
        end
        ripfile = sprintf('%s/EEG/%sripple%02d-%01d-%02d.mat', directoryname, prefix, day,epoch,currtet);
        load(ripfile);
        ripamp = ripple{day}{epoch}{currtet}.data(:,1);
        ripenv = ripple{day}{epoch}{currtet}.data(:,3);
    else
        EEGnostimfile = sprintf('%s/EEG/%seegnostim%02d-%01d-%02d.mat', directoryname, prefix, day,epoch,tets(te));
        load(EEGnostimfile);
        e = eeg{day}{epoch}{currtet};
        if te==1
            ti = geteegtimes(e);
            eind = lookup(pt(:,1), ti);
            e.samprate=round(e.samprate);
        end     
        ripnostimfile = sprintf('%s/EEG/%sripplenostim%02d-%01d-%02d.mat', directoryname, prefix, day,epoch,currtet);
        load(ripnostimfile);
        ripamp = ripple{day}{epoch}{currtet}.data(:,1);
        ripenv = ripple{day}{epoch}{currtet}.data(:,3);  
    end % end nostim   
    
    
             
    % Align EEG and Ripple Band to ripple time
    %------------------------------------------
   
    nelements = length(1000-round((pret/1000)*e.samprate):1000+round((postt/1000)*e.samprate));
    for i=1:length(pt) % Need to Skip initial and final indices?
        i;
        cnteeg=cnteeg+1;
        currriptime = pt(i); currripsize = ripsize(i);
        currind = eind(i);
        if ( (currind-round((pret/1000)*e.samprate) <=0) || (currind+round((postt/1000)*e.samprate)>length(e.data)) )
            e_stim{te}(cnteeg,:)=0*(1:nelements);
            ripamp_stim{te}(cnteeg,:)=0*(1:nelements);
            ripenv_stim{te}(cnteeg,:)=0*(1:nelements);
        else
            e_stim{te}(cnteeg,:)=e.data(currind-round((pret/1000)*e.samprate):currind+round((postt/1000)*e.samprate));
            ripamp_stim{te}(cnteeg,:)=double(ripamp(currind-round((pret/1000)*e.samprate):currind+round((postt/1000)*e.samprate)));
            ripenv_stim{te}(cnteeg,:)=double(ripenv(currind-round((pret/1000)*e.samprate):currind+round((postt/1000)*e.samprate)));
            ripsize_stim(cnteeg,te) = currripsize; 
        end
    end
    
end % end tets

% Get time ranges from times file - No need
%----------------------------------
cd(dire);
dirlist = dir('*');
for i=2:length(dirlist)
    currname = dirlist(i).name;
    if strcmp(currname(1:2),daystring)
        daydir = currname;
        break;
    end
end
timesfile = sprintf('%s/%s/times.mat',dire,daydir);
load(timesfile);
nranges=1;
range1=ranges(epoch+1,:);


% ----------------------------------------------------------------
%% Plot example single-trial LFPs aligned to ripples
% ----------------------------------------------------------------

nplots=length(tets);
taxis = [1:size(ripenv_stim{1},2)]*1000/e.samprate;
taxis = taxis-pret;

if iscon==1, 
    lineat=150; 
else
    lineat=100;
end
    
for i=1:size(e_stim{1},1)
    
    % Only do if atleast one ripple is bigger than threshold sd
    if max(ripsize_stim(i,:)) >= sd
        figure; hold on;
        redimscreen_halfvert;
        
        for n=1:nplots
            eegtet = e_stim{n};
            riptet = ripamp_stim{n};
            ripenvtet = ripenv_stim{n};
            curreeg = eegtet(i,:);
            currrip = riptet(i,:);
            currsize = roundn(ripsize_stim(i,n),-2);
            % EEg
            subplot(nplots,2,2*n-1); hold on;
            plot(taxis,curreeg,'k-','LineWidth',2);
            axis tight
            [ylimits] = get(gca,'YLim'); y1=ylimits(1); y2=ylimits(2);
            ypts = ylimits(1):0.1:ylimits(2);
            xpts = 0*ones(size(ypts)); % Plot Line at 0ms
            plot(xpts , ypts, 'r--','Linewidth',2);
            xpts = lineat*ones(size(ypts)); % Plot lines at 100ms
            plot(xpts , ypts, 'r:','Linewidth',2);
            set(gca,'XLim',[-pret postt]);
            if n==1
                title(['Stim No: ' num2str(i)]);
            end
            % Rippleband
            subplot(nplots,2,2*n); hold on;
            plot(taxis,currrip,'k-','LineWidth',2);
            axis tight
            [ylimits] = get(gca,'YLim'); y1=ylimits(1); y2=ylimits(2);
            ypts = ylimits(1):0.1:ylimits(2);
            xpts = 0*ones(size(ypts)); % Plot Line at 0ms
            plot(xpts , ypts, 'r--','Linewidth',2);
            xpts = lineat*ones(size(ypts)); % Plot lines at 100ms
            plot(xpts , ypts, 'r:','Linewidth',2);
            set(gca,'XLim',[-pret postt]);
            title(['Rip Size: ' num2str(currsize)]);
            
            %         if n==nplots,
            %             xlabel('Time(ms)','FontSize',xfont,'Fontweight','normal');
            %         else
            %             axis off;
            %         end
            
        end
        keyboard; % Pause after each plot
        
    end
end

if saveg1==1
    figfile = [figdir,'n',prefix,'_Day',num2str(day),'Ep',num2str(epoch),'Stim',num2str(i),'_RipEg1'];
    print('-dpdf', figfile);
    print('-djpeg', figfile);
    saveas(gcf,figfile,'fig');
end
  
% Manually plot
tet = 1; stim=75;
figure; hold on; 
eegtet=e_stim{tet};
x=eegtet(stim,:); % eg 75 tet4
plot(taxis,x,'k-','LineWidth',2);
set(gca,'XLim',[-50 50]);
figfile = [figdir,prefix,'_Day',num2str(day),'Ep',num2str(epoch),'Tet',num2str(tet),'Stim',num2str(stim),'_lfpeg1'];
print('-dpdf', figfile);
print('-djpeg', figfile);
saveas(gcf,figfile,'fig');

figure; hold on; 
riptet=ripamp_stim{tet};
y=riptet(stim,:); % eg 75 tet4
plot(taxis,y,'k-','LineWidth',2);
set(gca,'XLim',[-50 50]);

figfile = [figdir,prefix,'_Day',num2str(day),'Ep',num2str(epoch),'Tet',num2str(tet),'Stim',num2str(stim),'_ripeg1'];
print('-dpdf', figfile);
print('-djpeg', figfile);
saveas(gcf,figfile,'fig');

% *************************************************************************
cd(datadir);
keyboard;

