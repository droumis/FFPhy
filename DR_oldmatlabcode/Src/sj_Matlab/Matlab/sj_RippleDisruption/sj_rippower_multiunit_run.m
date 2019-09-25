function [basefr0] = sj_rippower_multiunit_run(prefixes, days, epochs, tets, nostim, figopt_single, saveg_single, figopt1, saveg1)
% Loop over and get z-normalized MU rate and Rip Power
% If only 1 session is entered, then do plot for that
% Plot Ripple Power and multinunit fr aligned to stimulation
% Can also align to ripples if you set dorip to 1
%
% Shantanu - derived from sj_rippower_stimampl3.
% This one only for run, no calibration Curve
% Autoload times file to get time ranges
% What tetrodes to use - can get from below or notes, or use riptetlist

% sj_rippower_multiunit_run('REe',9,4,[4 6 11],0,1,0);
% sj_rippower_multiunit_run('REe',9,4,11,0,1);
% sj_rippower_multiunit_run('REe',9,4,[2 4 6 11, 3 12 13],0,1);
% sj_rippower_multiunit_run('REd',9,4,[3 4 5 6 10 11],0,1);
% sj_rippower_multiunit_run('REf',9,4,[5 9 12, 1 10 11],0,1);

% Old eg
% sj_rippower_stimampl3('REd', 1,2,1,10,[3 4 5 6 10 11],0);
% sj_rippower_stimampl3('RCb', 4,4,5,1,10,[3 4 9 10 11 12],0);
% sj_rippower_stimampl3('RCb', 8,2,0,10,[],0);
% sj_rippower_stimampl3('RCb', 8,2,1,10,[3 4 9 10 11 12],0);

% Used for calibration
% sj_rippower_stimampl3('REd', 18,1,1,10,[3 4 5 6 10 11],0);
% sj_rippower_stimampl3('REd', 18,1,0,10,[],1);
% sj_rippower_stimampl3('REc', 20,1,0,9,[],0);
% sj_rippower_stimampl3('REc', 18,1,0,10,[],0);
% sj_rippower_stimampl3('REc', 19,1,0,10,[],1);
% sj_rippower_stimampl3('REc', 15,1,0,10,[],0);
% sj_rippower_stimampl3('REb', 27,1,1,4,[3 4 8 10 11 12],[],1);
% sj_rippower_stimampl3('RCa', 21,1,,0,4,[],1);
% sj_rippower_stimampl3('RCa', 21,1,,0,4,[],0);
% sj_rippower_stimampl3('RCb', 15,1,1,12,[3 4 9 10 11 12],0);
% sj_rippower_stimampl3('RCc', 17,1,1,13,[2 3 5 6 11 13],0);
% sj_rippower_stimampl3('RCc', 16,1,1,13,[2 3 5 6 11 13],0);


if nargin<1,
    keyboard
    error('Please enter Expt Prefix and Day No!');
end
if nargin<2,
    keyboard
    error('Please enter Day No!');
end
if nargin<3,
    epoch=[2 4]; %% Epoch - 2 or 4 for runs
end
if nargin<4,
    tets=[3 9 11]; %
end
if nargin<5
    nostim=0; % Set nostim=1 if you want to use artifact removed files for EEG
end
if nargin<6
    figopt_single=0; % Plot figure if single prefix,day,epoch is entered
end
if nargin<7
    saveg_single=0; % Save figure if single prefix,day,epoch is entered
end
if nargin<8
    figopt1=0; % Plot summary figure
end
if nargin<9
    saveg1=0; % Save summary figure
end

% --------------- Parameters ---------------
rem_mustimart=1; % Remove stim artifact in MU and rippower firing while plotting - for exp and con, not nor
iscon=0; % Whether control or exptal stimulation - for rippower alignment

plot_ex=0; % For plotting example single-trial MU rate plots
dorip = 0; % SET TO 1 IF YOU ALSO WANT TO ALIGN TO RIPPLES - Not implemented for across sessions yet
sd=3; %% SD for ripples
pret=220; postt=230; %% Times to plot
binsize = 5;  %% ms, for MU Fir Rate
binsize_plot=binsize; %% If you put binsize_plot=1000, then units are Nspikes/binsize, not inst. firing rate in Hz

% ------------------------------
% Figure and Font Sizes

forppr = 1;
% If yes, everything set to redimscreen_figforppr1
% If not, everything set to redimscreen_figforppt1

figdir = '/data25/sjadhav/RippleInterruption/Figures/01AugSep11_RippleDisFigs/DisruptnCalibrationAndEgs/';
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

if (size(prefixes,1)>1 || length(days)>1 || length(epochs)>1)
    figopt_single=0; saveg_single=0;
end
if (size(prefixes,1)==1 && length(days)==1 && length(epochs)==1)
    figopt1=0; saveg1=0;
end

cntsess=0; daysempty=0; tetsempty=0;
for an = 1:size(prefixes,1)
    
    prefix = prefixes(an,:);
    
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
            riptetlist = [1,5,7];
            usedays=1:8;
        case 'REc'
            directoryname = '/data25/sjadhav/RippleInterruption/REc_direct';
            dire = '/data25/sjadhav/RippleInterruption/REc';
            riptetlist = [4,6,8,9,10];
            usedaytets = [4,6,8,9,10];
            usedays=1:8; % MUrate
        case 'REd'
            directoryname = '/data25/sjadhav/RippleInterruption/REd_direct';
            dire = '/data25/sjadhav/RippleInterruption/REd';
            riptetlist = [3,4,5,6,10,11];
            usedaytets = [3,4,5,6,10,11];
            usedays=[1:8]; %1,4,7,8 or 1,8
        case 'REe'
            directoryname = '/data25/sjadhav/RippleInterruption/REe_direct';
            dire = '/data25/sjadhav/RippleInterruption/REe';
            riptetlist = [3,4,6,11,12,13];
            usedaytets = [3,4,6,11,12,13];
            usedays=[1,3]; % MUrate
            usedays=[1:8]; % rippower
        case 'REf'
            directoryname = '/data25/sjadhav/RippleInterruption/REf_direct';
            dire = '/data25/sjadhav/RippleInterruption/REf';
            riptetlist = [1,5,9,10,11,12];
            usedaytets = [1,5,9,10,11,12];
            usedays=1:8;
        case 'REg'
            directoryname = '/data25/sjadhav/RippleInterruption/REg_direct/StimAmpl';
            dire = '/data25/sjadhav/RippleInterruption/REg/StimAmpl';
            riptetlist = [2,5];
            usedaytets = [2,5];
            usedays=1:8;
        case 'RCa'
            directoryname = '/data25/sjadhav/RippleInterruption/RCa_direct';
            dire = '/data25/sjadhav/RippleInterruption/RCa';
            riptetlist = [2,3,4,6,9];
        case 'RCb'
            directoryname = '/data25/sjadhav/RippleInterruption/RCb_direct';
            dire = '/data25/sjadhav/RippleInterruption/RCb';
            riptetlist = [3,4,9,10,11,12];
            usedaytets = [3,4,9,10,11,12];
            usedays=[1:8];
        case 'RCc'
            directoryname = '/data25/sjadhav/RippleInterruption/RCc_direct';
            dire = '/data25/sjadhav/RippleInterruption/RCc';
            riptetlist = [3,4,5,6,11,13];
            usedaytets = [3,4,5,6,11,13];
            usedays = [1:8];
        case 'RCd'
            directoryname = '/data25/sjadhav/RippleInterruption/RCd_direct';
            dire = '/data25/sjadhav/RippleInterruption/RCd';
            riptetlist = [1,2,3,4,5,6];
            usedaytets = [1,2,3,4,5,6];
            usedays=[1:8]; % For MU, only 1 for rip power
    end
    
    currdir = pwd;
    if (directoryname(end) == '/')
        directoryname = directoryname(1:end-1);
    end
    if (dire(end) == '/')
        dire = dire(1:end-1);
    end
    
    if isempty(days)
        daysempty=1;
        days=usedays;
    end
    
    for d = 1:length(days)
        day = days(d)
        
        if (day < 10)
            daystring = ['0',num2str(day)];
        else
            daystring = num2str(day);
        end
        
        for ep=1:length(epochs)
            epoch = epochs(ep);
            if isempty(tets)
                tetsempty=1;
                tets = usedaytets;
            end
            
            cntsess = cntsess+1;
            %% For LFP
            eeg_pre = []; eegnostim_pre = []; rip_pre=[]; ripnostim_pre=[];
            eeg_post = []; eegnostim_post = []; rip_post=[]; ripnostim_post=[];
            eeg_run = []; eegnostim_run = []; rip_run=[]; ripnostim_run=[];
            dio_pre = [];
            s_pre = []; s_run = []; s_post = [];
            t_pre = []; t_post = []; t_run=[];
            
            ripamp=[]; ripenv=[];
            %% --------------------------------------------------
            % Align MU Firing Rate to Stimulations and Ripples
            % -------------------------------------------------
            
            % Set Counters
            cnt=0; % for MU spikes,
            cntrip=0; cnteeg=0;
            
            % Load dio file
            %------------------
            DIOfile = sprintf('%s/%sDIO%02d.mat', directoryname, prefix, day);
            load(DIOfile);
            stim = DIO{day}{epoch}{16};
            % if isempty(stim)
            %             stim = DIO{day}{epoch}{15};
            %         end
            stim_starttime = stim.pulsetimes(:,1)./10; %ms
            stim_endtime = stim.pulsetimes(:,2)./10; %ms
            stim_length = stim.pulselength;
            stim_isi = stim.timesincelast(2:end)./10; %ms
            
            % Load EEG and ripple LFP file
            %-------------------------
            if nostim==0 % stimulation artifact not removed
                
                EEGfile = sprintf('%s/EEG/%seeg%02d-%01d-%02d.mat', directoryname, prefix, day,epoch,tets(1));
                load(EEGfile);
                e = eeg{day}{epoch}{tets(1)};
                t = geteegtimes(e);
                pt = stim.pulsetimes ./ 10000; % in s
                eind = lookup(pt(:,1), t);
                e.samprate=round(e.samprate);
                
                for t=1:length(tets),
                    currtet=tets(t);
                    ripfile = sprintf('%s/EEG/%sripple%02d-%01d-%02d.mat', directoryname, prefix, day,epoch,currtet);
                    load(ripfile);
                    ripamp(t,:) = ripple{day}{epoch}{currtet}.data(:,1);
                    ripenv(t,:) = ripple{day}{epoch}{currtet}.data(:,3);
                end
            else
                
                EEGnostimfile = sprintf('%s/EEG/%seegnostim%02d-%01d-%02d.mat', directoryname, prefix, day,epoch,tets(1));
                load(EEGnostimfile);
                e = eeg{day}{epoch}{tets(1)};
                t = geteegtimes(e);
                pt = stim.pulsetimes ./ 10000;  % in s
                eind = lookup(pt(:,1), t);
                e.samprate=round(e.samprate);
                
                for t=1:length(tets),
                    currtet=tets(t);
                    ripnostimfile = sprintf('%s/EEG/%sripplenostim%02d-%01d-%02d.mat', directoryname, prefix, day,epoch,currtet);
                    load(ripnostimfile);
                    ripamp(t,:) = ripple{day}{epoch}{currtet}.data(:,1);
                    ripenv(t,:) = ripple{day}{epoch}{currtet}.data(:,3);
                end
            end % end nostim
            
            % Load MU file
            %------------------
            multifile = sprintf('%s/%smulti%02d.mat', directoryname, prefix, day);
            load(multifile);
            for t=1:length(tets),
                currtet=tets(t);
                cmd=sprintf('multi%d = multi{day}{epoch}{%d}/10;',currtet,currtet); eval(cmd);
            end
            
            %% Set which multiunit firing rate AND which Ripple Power to use
            %---------------------------------------------
            
            ripampu=sum(ripamp,1); % Sum across tetrodes if multiple exist
            ripenvu=sum(ripenv,1);
            multiu=[];
            for t=1:length(tets),
                currtet=tets(t);
                cmd=sprintf('curr_multi = multi%d;',currtet); eval(cmd);
                multiu=[multiu;curr_multi];
                % eg. multiu=[multi3; multi4; multi5; multi10; multi11; multi12];
            end
            tetu=tets;

            % Align spikes to stimulation
            %---------------------------
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
            
            % Align EEG and Ripple Power to stimulation
            %------------------------------------------
            for i =1:length(stim_starttime) % Need to Skip initial and final indices?
                i;
                cnteeg=cnteeg+1;
                currstim = stim_starttime(i); currind = eind(i);
                
                nelements = length(1000-round((pret/1000)*e.samprate):1000+round((postt/1000)*e.samprate));
                if ( (currind-round((pret/1000)*e.samprate) <=0) || (currind+round((postt/1000)*e.samprate)>length(e.data)) )
                    e_stim(cnteeg,:)=0*(1:nelements);
                    ripamp_stim(cnteeg,:)=0*(1:nelements);
                    ripenv_stim(cnteeg,:)=0*(1:nelements);
                else
                    e_stim(cnteeg,:)=e.data(currind-round((pret/1000)*e.samprate):currind+round((postt/1000)*e.samprate));
                    ripamp_stim(cnteeg,:)=double(ripampu(currind-round((pret/1000)*e.samprate):currind+round((postt/1000)*e.samprate)));
                    ripenv_stim(cnteeg,:)=double(ripenvu(currind-round((pret/1000)*e.samprate):currind+round((postt/1000)*e.samprate)));
                end
                
            end
            
            % Align spike to ripples if asked to
            %---------------------------------
            if dorip==1
                
                % Load extracted ripple file
                ripfile = sprintf('%s/%sripples%02d.mat', directoryname, prefix, day);
                load(ripfile);
                % Individual for each ripple tetrode
                for rt=1:length(riptetlist)
                    currriptet = riptetlist(rt);
                    rip_starttime{currriptet} = 1000*ripples{day}{epoch}{currriptet}.starttime;   % in msec
                end
                
                % Get ripple starttime across all tetrodes
                allripfile = sprintf('%s/%sripplesalltet%02d.mat', directoryname, prefix, day);
                load(allripfile);
                ripalltet_starttime = 1000*ripples{day}{epoch}{1}.starttime; % 1 tet condition
                %ripalltet_starttime = 1000*ripples{day}{epoch}{2}.starttime; % 2 tet condition
                
                % Either use the ripalltetstarttime, or pick a tetrode
                % userip_starttime = ripalltet_starttime;
                userip_starttime = rip_starttime{2};
                for i =2:length(userip_starttime)-1
                    i;
                    cntrip=cntrip+1;
                    currrip = userip_starttime(i);
                    currspks =  multiu(find( (multiu>=(currrip-pret)) & (multiu<=(currrip+postt)) ));
                    currspks = currspks-(currrip-pret);
                    histspks = histc(currspks,[0:binsize:pret+postt]);
                    rip_spks{cntrip}=currspks;
                    rip_spkshist(cntrip,:) = histspks;
                end
            end
            
            % Get time ranges from times file
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
            
            % There is only 1 range. So no need to divide by stim ampl like you do for calibration
            % -------------------------------------------------------------------------
            % Each stimulus start time is converted to 0.1ms resolution to get indexright
            
            stim_starttime_comp = stim_starttime*10; % Convert from ms to 10000 pts in sec
            % LFP
            eeghistall = e_stim;
            ripamp_histall = ripamp_stim;
            ripenv_histall = ripenv_stim;
            % Spikes
            stimhistall = stim_spkshist;
            if rem_mustimart==1
                stimhistall(:,(pret/binsize)+1)=0; % Stim artifact=0
                stimhistall(:,(pret/binsize)+2)=0;
                stimhistall(:,(pret/binsize)+3)=0;
            end
            
            % Save for combining if there are multiple sessions
            % --------------------------------------------------
            % MU Rate
            % --------
            MUrate_mean = mean(stimhistall,1)*(1000/binsize_plot); % in Hz
            MUrate_err = sem(stimhistall,1)*(1000/binsize_plot);
            
            % Get baselines and z-score data
            if iscon==0,
                basedata = MUrate_mean(1:100/binsize); % 1st 100ms for Exp and Nor
                %basedata = MUrate_mean(1:(pret/2)/binsize);
            else
                basedata = MUrate_mean((pret+130)/binsize:(pret+230)/binsize); % 130-230ms after stimulation for Con
                %basedata = MUrate_mean((3*pret/2)/binsize:(2*pret)/binsize);
            end
            basefr0 = mean(basedata);
            stdfr0 = std(basedata);
            MUrate_z(cntsess,:) = (MUrate_mean - basefr0)./stdfr0;
            MUrate_err_z(cntsess,:) = abs((MUrate_err - basefr0)./stdfr0); % sem in z-scores for current session
            if rem_mustimart==1
                MUrate_z(cntsess,(pret/binsize)+1:(pret/binsize)+3)=0;
            end
            
            % Ripple Power
            % ------------
            rippower_mean = mean(ripenv_histall,1);
            rippower_err = sem(ripenv_histall,1);
            % Get baselines and z-score data
            if iscon==0,
                basedata = rippower_mean(1:round(0.1*e.samprate)); % 1st 100ms for Exp and Nor
                %basedata = rippower_mean(1:round(0.1*e.samprate));
            else
                basedata = rippower_mean(round((pret/1000+0.13)*e.samprate):round((pret/1000+0.23)*e.samprate)); % 130-230ms after stimulation for Con
                %basedata = rippower_mean(round(0.33*e.samprate):round(0.43*e.samprate));
            end
            basepow0 = mean(basedata);
            stdpow0 = std(basedata);
            %stdpow0 = std(basedata);
            Rippower_z(cntsess,:) = (rippower_mean - basepow0)./stdpow0;
            Rippower_err_z(cntsess,:) = abs((rippower_err - basepow0)./stdpow0); % sem in z-scores for current session
        end % end epochs  
    end % end days
    
    % At end of each animal
    if daysempty==1
        days=[];
    end
    if tetsempty==1
        tets=[];
    end
    
end % end prefixes




% ***************** FIGURES FOR SUMMARY **********************
% *************************************************************************

if figopt1==1
    % -------------------------------------------------
    %% Plot Ripple Power aligned to stimulation
    % ----------------------------------------------
    figure; hold on;
    redimscreen_figforppt1;
    %redimscreen_halfvert(0);
    % Plot
    % ----
    taxis = [1:size(ripenv_stim,2)]*1000/e.samprate;
    taxis = taxis-pret;
    stimptidx = round((pret/1000)*e.samprate); % 200ms=300pts, 250ms=375pts;50ms=75pts,40ms=60pts, 30ms=45pts; 5ms=8pts; 35ms=53 pts;
    yplot = mean(Rippower_z);
    yerr = std(Rippower_z);
   
    % -------------------------------------
    % Eliminating stimulation artifact in figure if it is regular file (not nostim)
    % ------------------------------------
    % For Exp
    if ((rem_mustimart == 1) && (nostim==0))
        if iscon==0
            plot(taxis(1:stimptidx-rempts), yplot(1:stimptidx-rempts),'k','Linewidth',4);
            plot(taxis(stimptidx+rempts:end-10), yplot(stimptidx+rempts:end-10),'k','Linewidth',4);
            jbfill(taxis(1:stimptidx-rempts),yplot(1:stimptidx-rempts)+yerr(1:stimptidx-rempts),yplot(1:stimptidx-rempts)-yerr(1:stimptidx-rempts),'k','k',0.3, 0.3);
            jbfill(taxis(stimptidx+rempts:end-10),yplot(stimptidx+rempts:end-10)+yerr(stimptidx+rempts:end-10),yplot(stimptidx+rempts:end-10)-yerr(stimptidx+rempts:end-10),'k','k',0.3, 0.3);
            % Lines
            plot(taxis(1:stimptidx-rempts),yplot(1:stimptidx-rempts)+yerr(1:stimptidx-rempts),'k--','Linewidth',2);
            plot(taxis(1:stimptidx-rempts),yplot(1:stimptidx-rempts)-yerr(1:stimptidx-rempts),'k--','Linewidth',2);
            plot(taxis(stimptidx+rempts:end-10),yplot(stimptidx+rempts:end-10)+yerr(stimptidx+rempts:end-10),'k--','Linewidth',2);
            plot(taxis(stimptidx+rempts:end-10),yplot(stimptidx+rempts:end-10)-yerr(stimptidx+rempts:end-10),'k--','Linewidth',2);
        else       
            plot(taxis(1:stimptidx-rempts), yplot(1:stimptidx-rempts),'k','Linewidth',4);
            plot(taxis(stimptidx+rempts:end-10), yplot(stimptidx+rempts:end-10),'k','Linewidth',4);
            jbfill(taxis(1:stimptidx-rempts),yplot(1:stimptidx-rempts)+yerr(1:stimptidx-rempts),yplot(1:stimptidx-rempts)-yerr(1:stimptidx-rempts),'k','k',1, 1);
            jbfill(taxis(stimptidx+rempts:end-10),yplot(stimptidx+rempts:end-10)+yerr(stimptidx+rempts:end-10),yplot(stimptidx+rempts:end-10)-yerr(stimptidx+rempts:end-10),'k','k',1, 1);
            % Lines
            plot(taxis(1:stimptidx-rempts),yplot(1:stimptidx-rempts)+yerr(1:stimptidx-rempts),'k--','Linewidth',2);
            plot(taxis(1:stimptidx-rempts),yplot(1:stimptidx-rempts)-yerr(1:stimptidx-rempts),'k--','Linewidth',2);
            plot(taxis(stimptidx+rempts:end-10),yplot(stimptidx+rempts:end-10)+yerr(stimptidx+rempts:end-10),'k--','Linewidth',2);
            plot(taxis(stimptidx+rempts:end-10),yplot(stimptidx+rempts:end-10)-yerr(stimptidx+rempts:end-10),'k--','Linewidth',2);
        end
    else
        plot(taxis, yplot,'k','Linewidth',4);
        jbfill(taxis,yplot+yerr,yplot-yerr,'k','k',1, 1);
        % Lines
        jbfill(taxis,yplot+yerr,'k--','Linewidth',2);
        jbfill(taxis,yplot-yerr,'k--','Linewidth',2);
    end
    plot(taxis, 0*ones(size(taxis)),'r--','Linewidth',2 ); % Line at 0 in z-score, basepow0 in raw
    % Axes
    %--------
    [ylimits] = get(gca,'YLim'); y1=ylimits(1); y2=ylimits(2);
    ypts = ylimits(1):0.1:ylimits(2);
    xpts = 0*ones(size(ypts)); % Plot Line at 0ms
    plot(xpts , ypts, 'r--','Linewidth',2);
    xpts = 100*ones(size(ypts)); % Plot lines at 100ms and 200ms
    plot(xpts , ypts, 'r:','Linewidth',2);
    xpts = 200*ones(size(ypts));
    plot(xpts , ypts, 'r:','Linewidth',2);
    title(['Summ-Ripple Power aligned to stimln'],'FontSize',tfont,'Fontweight','normal');
    ylabel(['Ripple Power z-score'],'FontSize',yfont,'Fontweight','normal');
    xlabel('Time(ms)','FontSize',xfont,'Fontweight','normal');
    set(gca,'XLim',[-pret-binsize postt+binsize]);
    %set(gca,'YLim',[-5 y2]);
    % Save Graph if asked to
    if saveg1==1,
        figfile = [figdir,'ExpSumm_RipPowerStim'];
        print('-dpdf', figfile);
        print('-djpeg', figfile);
        saveas(gcf,figfile,'fig');
    end
    
    
    % -------------------------------------------------
    %% Plot MU Fir Rate Hist Aligned to Stimulation
    % ----------------------------------------------
    figure; hold on;
    redimscreen_figforppt1;
    %redimscreen_halfvert(1);
    % Plot
    % ----
    taxis = [-pret+binsize:binsize:postt+binsize];
    yplot=mean(MUrate_z);
    yerr=std(MUrate_z);
    stimptidx = (pret/binsize); 
   
    % Skip n bins after stimulation for stim artifact
    plot(taxis(1:stimptidx), yplot(1:stimptidx),'k','Linewidth',4);
    jbfill(taxis(1:stimptidx),yplot(1:stimptidx)+yerr(1:stimptidx),yplot(1:stimptidx)-yerr(1:stimptidx),'k','k',1, 1);
    plot(taxis(stimptidx+ceil(15/binsize)+1:end-2), yplot(stimptidx+ceil(15/binsize)+1:end-2),'k','Linewidth',4);
    jbfill(taxis(stimptidx+ceil(15/binsize)+1:end-2),yplot(stimptidx+ceil(15/binsize)+1:end-2)+yerr(stimptidx+ceil(15/binsize)+1:end-2),yplot(stimptidx+ceil(15/binsize)+1:end-2)-yerr(stimptidx+ceil(15/binsize)+1:end-2),'k','k',1, 1);
    % Dotted lines also along stdev
    plot(taxis(1:stimptidx),yplot(1:stimptidx)+yerr(1:stimptidx),'k--','Linewidth',2);
    plot(taxis(1:stimptidx),yplot(1:stimptidx)-yerr(1:stimptidx),'k--','Linewidth',2);
    plot(taxis(stimptidx+ceil(15/binsize)+1:end-2),yplot(stimptidx+ceil(15/binsize)+1:end-2)+yerr(stimptidx+ceil(15/binsize)+1:end-2),'k--','Linewidth',2);
    plot(taxis(stimptidx+ceil(15/binsize)+1:end-2),yplot(stimptidx+ceil(15/binsize)+1:end-2)-yerr(stimptidx+ceil(15/binsize)+1:end-2),'k--','Linewidth',2);
    
    plot(taxis, 0*ones(size(taxis)),'r--','Linewidth',2 );  % Line at 0 in z-score, basefr0 in raw
    % Axes
    %--------
    [ylimits] = get(gca,'YLim');
    y1=ylimits(1); y2=ylimits(2);
    ypts = ylimits(1):0.1:ylimits(2);
    xpts = 0*ones(size(ypts)); % Plot Line at 0ms
    xpts = 0*ones(size(ypts)); % Plot Line at 0ms
    plot(xpts , ypts, 'r--','Linewidth',2);
    xpts = 100*ones(size(ypts));
    plot(xpts , ypts, 'r:','Linewidth',2); % Plot lines at 100ms and 200ms
    xpts = 200*ones(size(ypts));
    plot(xpts , ypts, 'r:','Linewidth',2);
    title(['Summ-MU Rate aligned to stimln'],'FontSize',tfont,'Fontweight','normal');
    ylabel(['MU Fir Rate z-score in ',num2str(binsize),'ms bins'],'FontSize',yfont,'Fontweight','normal');
    xlabel('Time(ms)','FontSize',xfont,'Fontweight','normal');
    %text(-pret+binsize,0.8*y2,['Fir Rate ' num2str(basefr0) 'Hz'], 'FontSize',10,'Fontweight','normal')
    set(gca,'XLim',[-pret-binsize postt+binsize]);
    %set(gca,'YLim',[-5 y2]);
    
    % Save Graph if asked to
    if saveg1==1,
        figfile = [figdir,'ExpSumm_MURateStim_Day'];
        print('-dpdf', figfile);
        print('-djpeg', figfile);
        saveas(gcf,figfile,'fig');
    end
    
    
end % end figopt1


% ***************** FIGURES FOR SINGLE PREFIX,DAY,EP **********************
% *************************************************************************

if figopt_single==1
    
    % -------------------------------------------------
    %% Plot Ripple Power aligned to stimulation
    % ----------------------------------------------
    
    figure; hold on;
    redimscreen_figforppt1;
    %redimscreen_halfvert(0);
    
    % Plot
    % ----
    taxis = [1:size(ripenv_stim,2)]*1000/e.samprate;
    taxis = taxis-pret;
    yplot=Rippower_z;
    %taxis = taxis - (pret/1000)*e.samprate;
    %set(gca,'XLim',[-pret-binsize postt+binsize]);
    
    stimptidx = round((pret/1000)*e.samprate); % 200ms=300pts, 250ms=375pts;50ms=75pts,40ms=60pts, 30ms=45pts; 5ms=8pts; 35ms=53 pts;
    rempts = 60; % remove artifact pts
    % --------
    % Eliminating stimulation artifact in figure if it is regular file (not nostim)
    % ---------------------------------
    % For Exp
    if ((rem_mustimart == 1) && (nostim==0))
        if iscon==0
            plot(taxis(1:stimptidx-30), yplot(1:stimptidx-30),'k','Linewidth',4);
            plot(taxis(stimptidx+rempts:end), yplot(stimptidx+rempts:end),'k','Linewidth',4);
            %jbfill(taxis(1:stimptidx-8),yplot(1:stimptidx-8)+yerr(1:stimptidx-8),yplot(1:stimptidx-8)-yerr(1:stimptidx-8),'r','r',0.3, 0.3);
            %jbfill(taxis(stimptidx+53:end-10),yplot(stimptidx+53:end-10)+yerr(stimptidx+53:end-10),yplot(stimptidx+53:end-10)-yerr(stimptidx+53:end-10),'r','r',0.3, 0.3);
        else
            plot(taxis(1:stimptidx-rempts), yplot(1:stimptidx-rempts),'k','Linewidth',4);
            plot(taxis(stimptidx+rempts:end-5), yplot(stimptidx+rempts:end-5),'k','Linewidth',4);
            %jbfill(taxis(1:stimptidx-53),yplot(1:stimptidx-53)+yerr(1:stimptidx-53),yplot(1:stimptidx-53)-yerr(1:stimptidx-53),'r','r',0.3, 0.3);
            %jbfill(taxis(stimptidx+53:end-10),yplot(stimptidx+53:end-10)+yerr(stimptidx+53:end-10),yplot(stimptidx+53:end-10)-yerr(stimptidx+53:end-10),'r','r',0.3, 0.3);
        end
    else
        plot(taxis, yplot,'k','Linewidth',4);
        %jbfill(taxis,yplot+sem(ripenv_hist,1),yplot-sem(ripenv_hist,1),'r','r',0.3, 0.3);
    end
    
    plot(taxis,0*ones(size(taxis)),'r--','Linewidth',2 ); % Line at 0 in z-score, basepow0 in raw/unnor
    plot(taxis, (3)*ones(size(taxis)),'c--','Linewidth',2 ); % Line at n in z-score, basepow0+n*stdpow0 in raw/unnor
    plot(taxis, (4)*ones(size(taxis)),'c--','Linewidth',2 );
    plot(taxis, (5)*ones(size(taxis)),'c--','Linewidth',2 );
    % Also can plot the baseline got from file
    %plot(taxis, rippower_baseline*ones(size(taxis)),'k','Linewidth',2 );
    %plot(taxis, rippower_thresh*ones(size(taxis)),'c','Linewidth',2 );
    
    % Axes
    %--------
    [ylimits] = get(gca,'YLim'); y1=ylimits(1); y2=ylimits(2);
    ypts = ylimits(1):0.1:ylimits(2);
    xpts = 0*ones(size(ypts)); % Plot Line at 0ms
    plot(xpts , ypts, 'r--','Linewidth',2);
    xpts = 100*ones(size(ypts)); % Plot lines at 100ms and 200ms
    plot(xpts , ypts, 'r:','Linewidth',2);
    xpts = 200*ones(size(ypts));
    plot(xpts , ypts, 'r:','Linewidth',2);
    title([prefix,'-Ripple Power aligned to stimln - Tets ' num2str(tetu)],'FontSize',tfont,'Fontweight','normal');
    ylabel(['Ripple Power z'],'FontSize',yfont,'Fontweight','normal');
    xlabel('Time(ms)','FontSize',xfont,'Fontweight','normal');
    %text(-pret+0.05*e.samprate,0.97*y2,['Power ' num2str(basefr0)],'FontSize',10,'Fontweight','normal');
    set(gca,'XLim',[-pret-binsize postt+binsize]);
    %set(gca,'YLim',[-5 y2]);
    
    % Save Graph if asked to
    if saveg_single==1,
        if length(tetu)==1
            figfile = [figdir,prefix,'_RipPowerStimNor_Day',num2str(day),'Ep',num2str(epoch),'Tets',num2str(tetu)];
        else
            figfile = [figdir,prefix,'_RipPowerStimNor_Day',num2str(day),'Ep',num2str(epoch),'TetsAll1'];
        end
        print('-dpdf', figfile);
        print('-djpeg', figfile);
        saveas(gcf,figfile,'fig');
    end
    
    % -------------------------------------------------
    %% Plot MU Fir Rate Hist Aligned to Stimulation
    % ----------------------------------------------
    
    figure; hold on;
    redimscreen_figforppt1;
    %redimscreen_halfvert(1);
    
    % Plot
    % ----
    taxis = taxis - pret/binsize;
    yplot=MUrate_mean;
    bar(taxis, yplot,'FaceColor','k','EdgeColor','k');
    plot(taxis, basefr0*ones(size(taxis)),'r--','Linewidth',2 );  % Line at 0 in z-score, basefr0 in raw/unnor
    plot(taxis, (basefr0+3*stdfr0)*ones(size(taxis)),'c--','Linewidth',2 ); % Line at n in z-score, (basefr0+n*stdfr0) in raw/unnor
    plot(taxis, (basefr0+4*stdfr0)*ones(size(taxis)),'c--','Linewidth',2 );
    plot(taxis, (basefr0+5*stdfr0)*ones(size(taxis)),'c--','Linewidth',2 );
    % Axes
    %--------
    [ylimits] = get(gca,'YLim');
    y1=ylimits(1); y2=ylimits(2);
    ypts = ylimits(1):0.1:ylimits(2);
    xpts = 0*ones(size(ypts)); % Plot Line at 0ms
    xpts = 0*ones(size(ypts)); % Plot Line at 0ms
    plot(xpts , ypts, 'r--','Linewidth',2);
    xpts = 100*ones(size(ypts));
    plot(xpts , ypts, 'r:','Linewidth',2); % Plot lines at 100ms and 200ms
    xpts = 200*ones(size(ypts));
    plot(xpts , ypts, 'r:','Linewidth',2);
    title([prefix,'-MU Rate aligned to stimln - Tets ' num2str(tetu)],'FontSize',tfont,'Fontweight','normal');
    ylabel(['MU Fir Rate in ',num2str(binsize),'ms bins'],'FontSize',yfont,'Fontweight','normal');
    xlabel('Time(ms)','FontSize',xfont,'Fontweight','normal');
    %text(-pret+binsize,0.8*y2,['Fir Rate ' num2str(basefr0) 'Hz'], 'FontSize',10,'Fontweight','normal')
    set(gca,'XLim',[-pret-binsize postt+binsize]);
    %set(gca,'YLim',[-5 y2]);
    
    % Save Graph if asked to
    if saveg_single==1,
        if length(tetu)==1
            figfile = [figdir,prefix,'_MURateStimNor_Day',num2str(day),'Ep',num2str(epoch),'Tets',num2str(tetu)];
        else
            figfile = [figdir,prefix,'_MURateStimUnNor_Day',num2str(day),'Ep',num2str(epoch),'TetsAll1'];
        end
        print('-dpdf', figfile);
        print('-djpeg', figfile);
        saveas(gcf,figfile,'fig');
    end
    
    
    
    % ----------------------------------------------------------------
    %% Plot example single-trial firing rate and ripple power aligned to stimulation
    % ----------------------------------------------------------------
    
    
    if plot_ex==1
        figure; hold on;
        orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
        redimscreen_figforppt1;
        [ns]=[1,twoplot_idx]; % 0uA and 2nd amplitude
        for i=1:length(ns),
            n=ns(i);
            cmd=sprintf('stimhistall = stimhistall%d;',n); eval(cmd);
            cmd=sprintf('amp = amp%d;',n); eval(cmd);
            
            % Get 5 random indices from trials in stimhistall
            vec=1:size(stimhistall,1);
            r=randperm(length(vec));
            tr=vec(r(1:5));
            for p=1:5
                subplot(5,length(ns),2*(p-1)+i); hold on;
                stimhist=stimhistall(tr(p),:);
                yplot = (1000/binsize_plot)*stimhist; %Multiunit fr in "binsize"ms bins
                taxis = [-pret:binsize:postt];
                bar(taxis, yplot,'r');
                set(gca,'XLim',[-pret-binsize postt+binsize]);
                if n==1
                    [ylimits] = get(gca,'YLim');
                    y1=ylimits(1); y2=ylimits(2);
                    basefr0 = mean(yplot(1:(pret/2)/binsize));
                    basefr0 = round(10*basefr0)./10; % Because roundn does not exist in Lab matlab
                    std0 = std(yplot(1:(pret/2)/binsize));
                    std0 = roundn(std0,-1);
                    plot(taxis, (basefr0+3*std0)*ones(size(taxis)),'c-.','Linewidth',2 );
                    text(-pret+binsize,0.8*y2,['Fir Rate ' num2str(basefr0) 'Hz'], 'FontSize',14,'Fontweight','normal');
                end
                % Plot Line at 0ms
                ypts = 0:1:y2;
                xpts = 0*ones(size(ypts));
                % Plot lines at 100ms and 200ms
                plot(xpts , ypts, 'k-.','Linewidth',2);
                % Plot lines at 100ms and 200ms
                xpts = 100*ones(size(ypts));
                plot(xpts , ypts, 'k:','Linewidth',3);
                xpts = 200*ones(size(ypts));
                plot(xpts , ypts, 'k:','Linewidth',3);
                % Set Yaxis based on 0uA
                set(gca,'YLim',[y1 y2]);   % Same axis - All Graphs, based on 0uA
                % Plot Baseline FR and STD
                plot(taxis, basefr0*ones(size(taxis)),'k--','Linewidth',2 );
                if n~=1,
                    curr_basefr = mean(yplot(1:(pret/2)/binsize));
                    curr_basefr = round(10*curr_basefr)./10; % Because roundn does not exist in Lab matlab
                    cstd=std(yplot(1:(pret/2)/binsize));
                    cstd=roundn(cstd,-1);
                    plot(taxis, curr_basefr*ones(size(taxis)),'b--','Linewidth',2 );
                    plot(taxis, (curr_basefr+3*cstd)*ones(size(taxis)),'g-.','Linewidth',2 );
                    text(-pret+binsize,0.8*y2,['Fir Rate ' num2str(curr_basefr) 'Hz'],...
                        'FontSize',14,'Fontweight','normal','Color','b');
                end
                if p==1
                    title([num2str(amp) 'uA; ' num2str(binsize) 'ms bins; Tet ' num2str(tetu)],'FontSize',24,'Fontweight','normal');
                    %ylabel([ num2str(amp) 'uA'],'FontSize',14,'Fontweight','bold');
                end
                if p==5,
                    xlabel('Time(ms)','FontSize',24,'Fontweight','normal');
                end
            end
        end
    end
    
    % -------------------------------------
    %% Plot Multiunit rate around ripples
    % -------------------------------------
    
    if dorip==1
        %plot(taxis, 2*ripnostim_pre(i,:),['r-'],'Linewidth',2,'Markersize',6);
        ypts = 0:1.1*max(yplot);
        xpts = 200*ones(size(ypts));
        plot(xpts , ypts, 'k--','Linewidth',2);
        title(['Multiunit Firing around ripples-Day' num2str(day)],...
            'FontSize',24,'Fontweight','normal');
        %axis([0 800 -800 600]);
        ylabel('Instantaeous Multiunit Firing Rate','FontSize',24,'Fontweight','normal');
        xlabel('Time(ms)','FontSize',20,'Fontweight','normal');
        %text( 4, 2450,['DetRate(4,6,7) =' num2str(round(DetRateper*100)/100)
        %'%'],'FontSize', 24, 'FontWeight','bold');
        if saveg_single==1,
            if length(tetu)==1
                figfile = [figdir,prefix,'_MURateRip_Day',num2str(day),'Ep',num2str(epoch),'Tets',num2str(tetu)];
            else
                figfile = [figdir,prefix,'_MURateRip_Day',num2str(day),'Ep',num2str(epoch),'TetsAll1'];
            end
            print('-dpdf', figfile);
            print('-djpeg', figfile);
            saveas(gcf,figfile,'fig');
        end
    end
    
end % figopt_single

% *************************************************************************
cd(datadir);
keyboard;

