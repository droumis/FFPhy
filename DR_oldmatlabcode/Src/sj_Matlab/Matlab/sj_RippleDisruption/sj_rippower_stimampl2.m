function [basefr0] = sj_rippower_stimampl2(prefix, day, epoch, binsize, manytet, maintet, tets, nostim, saveg1)
%% Plot Ripple Power Calibration Curve - Change from version 1:
% Data can be calibration data, and now also data taken during runs. So
% include an epoch field in function input. Also, Get ripple power by
% combining from all tetrodes that were used in ripple detection.

% Be careful what you input! If Stimulator is on, get rid of artifact (ripple_nostim files)
% Options from Multiunit Data also in here (sj_multiunitdr_rip_stimampl1.m)
% eg
% sj_rippower_stimampl2('REd', 1,2,5,1,10,[3 4 5 6 10 11],0);
% sj_rippower_stimampl2('REd', 18,1,5,1,10,[3 4 5 6 10 11],0);
% sj_rippower_stimampl2('REd', 18,1,5,0,10,[],1);
% sj_rippower_stimampl2('REc', 20,1,5,0,9,[],0);
% sj_rippower_stimampl2('REc', 18,1,5,0,10,[],0);
% sj_rippower_stimampl2('REc', 19,1,5,0,10,[],1);
% sj_rippower_stimampl2('REc', 15,1,5,0,10,[],0);
% sj_rippower_stimampl2('REb', 27,1,10,1,4,[3 4 8 10 11 12],[],1);
% sj_rippower_stimampl2('RCa', 21,1,5,0,4,[],1);
% sj_rippower_stimampl2('RCa', 21,1,5,0,4,[],0);

if nargin<1,
    keyboard
    error('Please enter Expt Prefix and Day No!');
end
if nargin<2,
    keyboard
    error('Please enter Day No!');
end
if nargin<3,
    epoch=1; %% Epoch - Usually 1 for calibration
end
if nargin<4,
    binsize=5;
end
if nargin<5,
    manytet=0; % Set to 1 if you want combine MU firing across tetrodes, 0 otherwise
end
if nargin<6,
    maintet=4; % If you want to look at MU firing on single tetrode
end
if nargin<7,
    tets=[4 5]; % If you want to look at MU firing on multiple tetrodes combined
end
if nargin<8
    nostim=0; % Set nostim=1 if you want to use artifact removed files for EEG
end
if nargin<9
    saveg1=0;
end

% for the plot comparing 0uA with 2nd stimulation amplitude, choose which index
if epoch~=1
    twoplot_idx=1; % If Expt Day, always a single stimulation
else
    twoplot_idx=4;
end

plot_ex=0; % For plotting example single-trial MU rate plots

set(0,'defaultaxesfontsize',14);set(0,'defaultaxesfontweight','bold'); set(0,'defaultaxeslinewidth',2);



%% -----------------------------------------
% SET DATA
% -------------------------------------------

switch prefix
    case 'REb'
        directoryname = '/data25/sjadhav/RippleInterruption/REb_direct/StimAmpl';
    case 'RCa'
        directoryname = '/data25/sjadhav/RippleInterruption/RCa_direct/StimAmpl';
    case 'RE1'
        directoryname = '/data25/sjadhav/RippleInterruption/RE1_direct/StimAmpl';
    case 'REc'
        if day>10
            directoryname = '/data25/sjadhav/RippleInterruption/REc_direct/StimAmpl';
        else
            directoryname = '/data25/sjadhav/RippleInterruption/REc_direct';
        end
    case 'REd'
        if day>12
            directoryname = '/data25/sjadhav/RippleInterruption/REd_direct/StimAmpl';
        else
            directoryname = '/data25/sjadhav/RippleInterruption/REd_direct';
        end
        
    case 'REe'
        if day>12
            directoryname = '/data25/sjadhav/RippleInterruption/REe_direct/StimAmpl';
        else
            directoryname = '/data25/sjadhav/RippleInterruption/REe_direct';
        end
        
    case 'REf'
        if day>12
            directoryname = '/data25/sjadhav/RippleInterruption/REf_direct/StimAmpl';
        else
            directoryname = '/data25/sjadhav/RippleInterruption/REf_direct';
        end
    
     case 'REg'
        if day>12
            directoryname = '/data25/sjadhav/RippleInterruption/REg_direct/StimAmpl';
        else
            directoryname = '/data25/sjadhav/RippleInterruption/REg_direct';
        end
        
     case 'RCb'
        if day>12
            directoryname = '/data25/sjadhav/RippleInterruption/RCb_direct/StimAmpl';
        else
            directoryname = '/data25/sjadhav/RippleInterruption/RCb_direct';
        end
        
end

dorip = 0; % SET TO 1 IF YOU ALSO WANT TO ALIGN TO RIPPLES

sd=4; %% SD for ripples

pret=200; postt=250; %% Times to plot
%binsize = 5;  %% ms, for MU Fir Rate
binsize_plot=binsize; %% If you put binsize_plot=1000, then units are Nspikes/binsize, not inst. firing rate in Hz

%% For LFP
eeg_pre = []; eegnostim_pre = []; rip_pre=[]; ripnostim_pre=[];
eeg_post = []; eegnostim_post = []; rip_post=[]; ripnostim_post=[];
eeg_run = []; eegnostim_run = []; rip_run=[]; ripnostim_run=[];
dio_pre = [];
s_pre = []; s_run = []; s_post = [];
t_pre = []; t_post = []; t_run=[];
clr = {'b','g','c','m','y','k','r'};

%% --------------------------------------------------
% Align MU Firing Rate to Stimulations and Ripples
% -------------------------------------------------

% Set Counters
cnt=0; cntrip=0; cnteeg=0;

% Load extracted ripple file
if dorip ==1
    ripfile = sprintf('%s/%sripples%02dstd%02d.mat', directoryname, prefix, day, sd);
    load(ripfile);
    rip_starttime = 1000* ripples{day}{epoch}{tet}.starttime;   % in msec
end

% Load dio file
DIOfile = sprintf('%s/%sDIO%02d.mat', directoryname, prefix, day);
load(DIOfile);
stim = DIO{day}{epoch}{16};
stim_starttime = stim.pulsetimes(:,1)./10; %ms
stim_endtime = stim.pulsetimes(:,2)./10; %ms
stim_length = stim.pulselength;
stim_isi = stim.timesincelast(2:end)./10; %ms

% Load EEG and ripple file
if nostim==0
    
    EEGfile = sprintf('%s/EEG/%seeg%02d-%01d-%02d.mat', directoryname, prefix, day,epoch,maintet);
    load(EEGfile);
    e = eeg{day}{epoch}{maintet};
    t = geteegtimes(e);
    pt = stim.pulsetimes ./ 10000; % in s
    eind = lookup(pt(:,1), t);
    
    e.samprate=round(e.samprate);
    
    if manytet==0
        ripfile = sprintf('%s/EEG/%sripple%02d-%01d-%02d.mat', directoryname, prefix, day,epoch,maintet);
        load(ripfile);
        ripamp = ripple{day}{epoch}{maintet}.data(:,1);
        ripenv = ripple{day}{epoch}{maintet}.data(:,3);
    else
        for t=1:length(tets),
            currtet=tets(t);
            ripfile = sprintf('%s/EEG/%sripple%02d-%01d-%02d.mat', directoryname, prefix, day,epoch,currtet);
            load(ripfile);
            ripamp(t,:) = ripple{day}{epoch}{currtet}.data(:,1);
            ripenv(t,:) = ripple{day}{epoch}{currtet}.data(:,3);
        end
    end
    
else
    
    EEGnostimfile = sprintf('%s/EEG/%seegnostim%02d-%01d-%02d.mat', directoryname, prefix, day,epoch,maintet);
    load(EEGnostimfile);
    e = eeg{day}{epoch}{maintet};
    t = geteegtimes(e);
    pt = stim.pulsetimes ./ 10000;  % in s
    eind = lookup(pt(:,1), t);
    
    e.samprate=round(e.samprate);
    
    if manytet==0
        ripnostimfile = sprintf('%s/EEG/%sripplenostim%02d-%01d-%02d.mat', directoryname, prefix, day,epoch,maintet);
        load(ripnostimfile);
        ripamp = ripple{day}{epoch}{maintet}.data(:,1);
        ripenv = ripple{day}{epoch}{maintet}.data(:,3);
    else
        for t=1:length(tets),
            currtet=tets(t);
            ripnostimfile = sprintf('%s/EEG/%sripplenostim%02d-%01d-%02d.mat', directoryname, prefix, day,epoch,currtet);
            load(ripnostimfile);
            ripamp(t,:) = ripple{day}{epoch}{currtet}.data(:,1);
            ripenv(t,:) = ripple{day}{epoch}{currtet}.data(:,3);
        end
    end
end


% Load MU file
multifile = sprintf('%s/%smulti%02d.mat', directoryname, prefix, day);
load(multifile);

% Get Multiunit timestamps (divide by 10 for ms)
cmd=sprintf('multi%d = multi{day}{epoch}{%d}/10;',maintet,maintet); eval(cmd);
for t=1:length(tets),
    currtet=tets(t);
    cmd=sprintf('multi%d = multi{day}{epoch}{%d}/10;',currtet,currtet); eval(cmd);
end
%multi3 = multi{day}{epoch}{3}/10 ; multi4 = multi{day}{epoch}{4}/10 ;

%% Set which multiunit firing rate AND which Ripple Power to use
if manytet==0
    ripampu=ripamp;
    ripenvu=ripenv;
    
    cmd=sprintf('multi_this = multi%d;',maintet); eval(cmd); % eg %multiu=[multi4];
    multiu=[multi_this];
    tetu=maintet;
else
    
    ripampu=sum(ripamp,1); % Sum across tetrodes
    ripenvu=sum(ripenv,1);
    
    multiu=[];
    for t=1:length(tets),
        currtet=tets(t);
        cmd=sprintf('curr_multi = multi%d;',currtet); eval(cmd);
        multiu=[multiu;curr_multi];
        % eg. multiu=[multi3; multi4; multi5; multi10; multi11; multi12];
    end
    tetu=tets;
end

% Align spikes to stimulation
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
for i =1:length(stim_starttime) % Need to Skip initial and final indices?
    i;
    cnteeg=cnteeg+1;
    currstim = stim_starttime(i); currind = eind(i);
    
    nelements = length(1000-round(0.2*e.samprate):1000+round(0.25*e.samprate));
    if ( (currind-round(0.2*e.samprate) <=0) || (currind+round(0.25*e.samprate)>length(e.data)) )
        e_stim(cnteeg,:)=0*(1:nelements);
        ripamp_stim(cnteeg,:)=0*(1:nelements);
        ripenv_stim(cnteeg,:)=0*(1:nelements);
    else
        e_stim(cnteeg,:)=e.data(currind-round(0.2*e.samprate):currind+round(0.25*e.samprate));
        ripamp_stim(cnteeg,:)=double(ripampu(currind-round(0.2*e.samprate):currind+round(0.25*e.samprate)));
        ripenv_stim(cnteeg,:)=double(ripenvu(currind-round(0.2*e.samprate):currind+round(0.25*e.samprate)));
    end
    
end

% Align spike to ripples if asked to
if dorip==1
    for i =2:length(rip_starttime)-1
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


%% Set ranges for each stimulation amplitude Or Other Parameter

switch prefix      % BY ANIMAL
    
     case 'REg'
        
        switch day
            
            case 16
                
                nranges=8;
                [range1] = timetrans({'00:03:58', '00:06:59'},10000,2); amp1=0; %0uA
                [range2] = timetrans({'00:07:05', '00:10:06'},10000,2); amp2=10; 
                [range3] = timetrans({'00:10:25', '00:13:26'},10000,2); amp3=20; 
                [range4] = timetrans({'00:13:40', '00:16:41'},10000,2); amp4=30; 
                [range5] = timetrans({'00:16:55', '00:19:56'},10000,2); amp5=40; 
                [range6] = timetrans({'00:20:05', '00:23:06'},10000,2); amp6=50; 
                [range7] = timetrans({'00:23:10', '00:26:11'},10000,2); amp7=60; 
                [range8] = timetrans({'00:26:30', '00:29:31'},10000,2); amp8=70; 
                
                
            case 17
                
                nranges=7;
                [range1] = timetrans({'00:10:49', '00:16:20'},10000,2); amp1=0; %0uA
                [range2] = timetrans({'00:16:29', '00:19:30'},10000,2); amp2=5; %5uA
                [range3] = timetrans({'00:19:39', '00:22:40'},10000,2); amp3=10; %10uA
                [range4] = timetrans({'00:23:15', '00:26:16'},10000,2); amp4=15; %15uA
                [range5] = timetrans({'00:26:30', '00:29:31'},10000,2); amp5=20; 
                [range6] = timetrans({'00:29:39', '00:32:40'},10000,2); amp6=25; 
                [range7] = timetrans({'00:32:50', '00:35:51'},10000,2); amp7=30; 
  
                
        end     
        
    
    case 'RCb'
        
        switch day
            
            
            case 1
                % Day1 - Control Stim On; Rip Detectn on Tets 3,4,9,10,11,12
                % 150ms delay
                
                nranges=1;
                switch epoch
                    case 2
                        [range1] = timetrans({'00:24:16', '00:39:17'},10000,2); amp1=100;
                    case 4
                        [range1] = timetrans({'00:58:57', '01:13:58'},10000,2); amp1=100;
                        
                end % end switch epoch
                
            case 15
                
                nranges=7;
                [range1] = timetrans({'05:02:20', '05:05:20'},10000,2); amp1=0;
                %[range1] = timetrans({'05:25:40', '05:28:54'},10000,2);%amp1=0;
                [range2] = timetrans({'05:05:25', '05:09:25'},10000,2); amp2=20;
                [range3] = timetrans({'05:09:30', '05:12:30'},10000,2); amp3=40;
                [range4] = timetrans({'05:12:40', '05:15:40'},10000,2); amp4=60;
                [range5] = timetrans({'05:15:50', '05:18:50'},10000,2); amp5=80;
                [range6] = timetrans({'05:19:05', '05:22:05'},10000,2); amp6=100;
                [range7] = timetrans({'05:22:35', '05:25:35'},10000,2); amp7=120;
                
            case 16
                
                % Con Dis 100 uA
                
                nranges=5;
                [range1] = timetrans({'00:18:26', '00:22:22'},10000,2); amp1=0;
                %[range1] = timetrans({'05:25:40', '05:28:54'},10000,2);%amp1=0;
                [range2] = timetrans({'00:06:02', '00:09:02'},10000,2); amp2=150;
                [range3] = timetrans({'00:09:05', '00:12:05'},10000,2); amp3=200;
                [range4] = timetrans({'00:12:20', '00:15:20'},10000,2); amp4=200; % +10ms jitter
                [range5] = timetrans({'00:15:25', '00:18:25'},10000,2); amp5=150;% +10ms jitter
                
            case 17
                
                % Con Dis 100 uA
                
                nranges=4;
                [range1] = timetrans({'00:02:27', '00:05:27'},10000,2); amp1=150; %100uA
                [range2] = timetrans({'00:05:30', '00:08:30'},10000,2); amp2=150; %0uA
                [range3] = timetrans({'00:08:35', '00:11:35'},10000,2); amp3=200; %100uA
                [range4] = timetrans({'00:11:36', '00:14:36'},10000,2); amp4=200; %0uA
                
            case 18
                
                % Con Dis 100 uA
                
                nranges=4;
                [range1] = timetrans({'01:27:54', '01:30:54'},10000,2); amp1=150; %100uA
                [range2] = timetrans({'01:30:55', '01:34:00'},10000,2); amp2=150; %0uA
                [range3] = timetrans({'01:34:17', '01:37:17'},10000,2); amp3=200; %100uA
                [range4] = timetrans({'01:37:20', '01:40:24'},10000,2); amp4=200; %0uA
                
        end
    
    
        
   
        
    case 'REf'
        
        switch day
            
            case 3
                % Day3 - Stim On; Rip Detectn on Tets
                
                nranges=1;
                switch epoch
                    case 2
                        [range1] = timetrans({'00:20:05', '00:35:36'},10000,2); amp1=200;
                    case 4
                        [range1] = timetrans({'00:54:11', '01:09:12'},10000,2); amp1=200;
                        
                end % end switch epoch
                
            case 4
                % Day4 - Stim On; Rip Detectn on Tets
                
                nranges=1;
                switch epoch
                    case 2
                        [range1] = timetrans({'00:31:08', '00:46:09'},10000,2); amp1=200;
                    case 4
                        [range1] = timetrans({'01:05:55', '01:20:56'},10000,2); amp1=200;
                        
                end % end switch epoch
                
            case 5
                % Day5 - Stim On; Rip Detectn on Tets
                
                nranges=1;
                switch epoch
                    case 2
                        [range1] = timetrans({'00:32:07', '00:47:08'},10000,2); amp1=200;
                    case 4
                        [range1] = timetrans({'01:05:20', '01:20:21'},10000,2); amp1=200;
                        
                end % end switch epoch
                
             case 8
                % Day8 - Stim On; Rip Detectn on Tets
                
                nranges=1;
                switch epoch
                    case 2
                        [range1] = timetrans({'00:27:00', '00:42:01'},10000,2); amp1=300;
                    case 4
                        [range1] = timetrans({'01:00:17', '01:17:40'},10000,2); amp1=300;
                        
                end % end switch epoch
                
             case 9
                % Day9 - Stim On; Rip Detectn on Tets
                
                nranges=1;
                switch epoch
                    case 2
                        [range1] = timetrans({'00:36:10', '00:51:31'},10000,2); amp1=300;
                    case 4
                        [range1] = timetrans({'01:09:40', '01:24:41'},10000,2); amp1=300;
                        
                end % end switch epoch
                
                
            case 15
                % 15ripdistest1 - Stim On; Rip Detectn on Tet 1
                % 200uA;  Nthr=1
                nranges=1;
                [range1] = timetrans({'00:22:05', '00:25:12'},10000,2); amp1=200;
                
            case 16
                % 16ripdistest2- Stim On; Rip Detectn on Tet 1
                % 200uA;  Nthr=1
                % Epoch = 1; 0uA 0:02:03-0:04:41
                % Epoch = 2; 6 amplitudes
                switch epoch
                    case 1
                        nranges=1;
                        [range1] = timetrans({'00:02:03', '00:04:41'},10000,2); amp1=0;
                    case 2
                        nranges=6;
                        [range1] = timetrans({'00:05:15', '00:08:00'},10000,2); amp1=50;
                        [range2] = timetrans({'00:08:24', '00:11:00'},10000,2); amp2=100;
                        [range3] = timetrans({'00:11:05', '00:13:40'},10000,2); amp3=150;
                        [range4] = timetrans({'00:13:45', '00:16:20'},10000,2); amp4=200;
                        [range5] = timetrans({'00:16:28', '00:19:05'},10000,2); amp5=250;
                        [range6] = timetrans({'00:19:15', '00:22:02'},10000,2); amp6=300;
                        
                        
                end % end switch epoch
        end
        
        
    case 'REe'
        
        switch day
            
            case 1
                % Day1 - Stim On; Rip Detectn on Tets 3,4,6,11,12,13
                
                nranges=1;
                switch epoch
                    case 2
                        [range1] = timetrans({'00:23:55', '00:38:56'},10000,2); amp1=200;
                    case 4
                        [range1] = timetrans({'00:57:15', '01:12:16'},10000,2); amp1=200;
                        
                end % end switch epoch
                
            case 2
                % Day2 - Stim On; Rip Detectn on Tets 3,4,6,11,12,13
                
                nranges=1;
                switch epoch
                    case 2
                        [range1] = timetrans({'00:39:45', '00:54:46'},10000,2); amp1=350;
                    case 4
                        [range1] = timetrans({'01:16:15', '01:31:16'},10000,2); amp1=350;
                        
                end % end switch epoch
                
            case 3
                % Day3 - Stim On; Rip Detectn on Tets 3,4,6,11,12,13
                
                nranges=1;
                switch epoch
                    case 2
                        [range1] = timetrans({'01:00:20', '01:15:21'},10000,2); amp1=350;
                    case 4
                        [range1] = timetrans({'01:35:50', '01:50:51'},10000,2); amp1=350;
                        
                end % end switch epoch
                
            case 4
                % Day4 - Stim On; Rip Detectn on Tets 3,4,6,11,12,13
                
                nranges=1;
                switch epoch
                    case 2
                        [range1] = timetrans({'00:31:00', '00:46:01'},10000,2); amp1=350;
                    case 4
                        [range1] = timetrans({'01:06:58', '01:21:59'},10000,2); amp1=350;
                        
                end % end switch epoch
                
            case 5
                % Day5 - Stim On; Rip Detectn on Tets 3,4,6,11,12,13
                
                nranges=1;
                switch epoch
                    case 2
                        [range1] = timetrans({'00:29:27', '00:44:35'},10000,2); amp1=350;
                    case 4
                        [range1] = timetrans({'01:05:12', '01:20:13'},10000,2); amp1=350;
                        
                end % end switch epoch
                
             case 6
                % Day6 - Stim On; Rip Detectn on Tets 3,4,6,11,12,13
                
                nranges=1;
                switch epoch
                    case 2
                        [range1] = timetrans({'00:34:25', '00:49:26'},10000,2); amp1=350;
                    case 4
                        [range1] = timetrans({'01:08:35', '01:23:36'},10000,2); amp1=350;
                        
                end % end switch epoch
                
            case 7
                % Day7 - Stim On; Rip Detectn on Tets 3,4,6,11,12,13
                
                nranges=1;
                switch epoch
                    case 2
                        [range1] = timetrans({'00:23:25', '00:38:33'},10000,2); amp1=350;
                    case 4
                        [range1] = timetrans({'00:56:30', '01:11:31'},10000,2); amp1=350;
                        
                end % end switch epoch
                
            case 8
                % Day8 - Stim On; Rip Detectn on Tets 3,4,6,11,12,13
                
                nranges=1;
                switch epoch
                    case 2
                        [range1] = timetrans({'00:30:15', '00:45:16'},10000,2); amp1=350;
                    case 4
                        [range1] = timetrans({'01:04:07', '01:19:08'},10000,2); amp1=350;
                        
                end % end switch epoch
                
            case 9
                % Day9 - Stim On; Rip Detectn on Tets 3,4,6,11,12,13
                
                nranges=1;
                switch epoch
                    case 2
                        [range1] = timetrans({'00:28:57', '00:43:58'},10000,2); amp1=350;
                    case 4
                        [range1] = timetrans({'01:03:59', '01:19:00'},10000,2); amp1=350;
                        
                end % end switch epoch
                
                
            case 15
                
                nranges=1;
                switch epoch
                    case 1
                        [range1] = timetrans({'04:31:32', '04:36:37'},10000,2); amp1=200;
                    case 2
                        [range1] = timetrans({'04:38:05', '04:42:10'},10000,2); amp1=200;
                    case 3
                        [range1] = timetrans({'04:43:35', '04:48:41'},10000,2); amp1=200;
                    case 4
                        [range1] = timetrans({'04:50:11', '04:55:16'},10000,2); amp1=200;
                        
                        
                end % end switch epoch
                
            case 16
                % 16ripdistest2 - Stim On; Rip Detectn on Tet 2
                % 200uA; SD=5/4; Nthr=1
                nranges=1;
                [range1] = timetrans({'00:22:25', '00:29:16'},10000,2); amp1=200;
                
            case 17
                % 17ripdistest3 - Stim On; Rip Detectn on Tet
                % 3,4,6,11,12,13
                % 200uA; SD=5; Nthr=2
                nranges=1;
                [range1] = timetrans({'01:23:48', '01:27:56'},10000,2); amp1=200;
                
            case 18
                % 18ripdistest4 - Stim On; Rip Detectn on Tet
                % 3,4,6,11,12,13
                % 200uA; SD=5; Nthr=2
                nranges=1;
                [range1] = timetrans({'00:06:46', '00:12:08'},10000,2); amp1=200;
                
                
        end
        
        
        
    case 'REd'
        
        switch day        % BY DAY WITHIN ANIMAL
            
            case 1
                % Day1 - Stim On; Rip Detectn on Tets 3,4,5,6,10,11
                
                nranges=1;
                switch epoch
                    case 2
                        [range1] = timetrans({'00:27:48', '00:42:49'},10000,2); amp1=80;
                    case 4
                        [range1] = timetrans({'01:04:52', '01:19:53'},10000,2); amp1=80;
                        
                end % end switch epoch
                
            case 2
                % Day2 - Stim On; Rip Detectn on Tets 3,4,5,6,10,11
                % 150uA; SD=5-6; Nthr=2
                
                nranges=1;
                switch epoch
                    case 2
                        [range1] = timetrans({'00:56:23', '01:11:24'},10000,2); amp1=150;
                    case 4
                        [range1] = timetrans({'01:35:10', '01:51:30'},10000,2); amp1=150;
                        
                end % end switch epoch
                
            case 3
                % Day3 - Stim On; Rip Detectn on Tets 3,4,5,6,10,11
                % 150uA; SD=5-6; Nthr=2
                
                nranges=1;
                switch epoch
                    case 2
                        [range1] = timetrans({'00:36:04', '00:51:05'},10000,2); amp1=150;
                    case 4
                        [range1] = timetrans({'01:10:55', '01:25:56'},10000,2); amp1=150;
                        
                end % end switch epoch
                
            case 4
                % Day4 - Stim On; Rip Detectn on Tets 3,4,5,6,10,11
                % 150uA; SD=5-6; Nthr=2
                
                nranges=1;
                switch epoch
                    case 2
                        [range1] = timetrans({'00:53:03', '01:08:04'},10000,2); amp1=150;
                    case 4
                        [range1] = timetrans({'01:32:23', '01:47:24'},10000,2); amp1=150;
                        
                end % end switch epoch
                
            case 5
                % Day5 - Stim On; Rip Detectn on Tets 3,4,5,6,10,11
                % 150uA; SD=5-6; Nthr=2
                
                nranges=1;
                switch epoch
                    case 2
                        [range1] = timetrans({'00:44:01', '00:59:02'},10000,2); amp1=150;
                    case 4
                        [range1] = timetrans({'01:18:08', '01:33:09'},10000,2); amp1=150;
                        
                end % end switch epoch
                
            case 6
                % Day6 - Stim On; Rip Detectn on Tets 3,4,5,6,10,11
                % 150uA; SD=5-6; Nthr=2
                
                nranges=1;
                switch epoch
                    case 2
                        [range1] = timetrans({'00:39:20', '00:54:21'},10000,2); amp1=150;
                    case 4
                        [range1] = timetrans({'01:14:40', '01:29:41'},10000,2); amp1=150;
                        
                end % end switch epoch
                
            case 7
                % Day7 - Stim On; Rip Detectn on Tets 3,4,5,6,10,11
                % 150uA; SD=5-6; Nthr=2
                
                nranges=1;
                switch epoch
                    case 2
                        [range1] = timetrans({'00:45:18', '01:00:19'},10000,2); amp1=150;
                    case 4
                        [range1] = timetrans({'01:21:32', '01:36:33'},10000,2); amp1=150;
                        
                end % end switch epoch
                
            case 8
                % Day8 - Stim On; Rip Detectn on Tets 3,4,5,6,10,11
                % 150uA; SD=5-6; Nthr=2
                
                nranges=1;
                switch epoch
                    case 2
                        [range1] = timetrans({'00:53:10', '01:08:11'},10000,2); amp1=150;
                    case 4
                        [range1] = timetrans({'01:30:00', '01:45:01'},10000,2); amp1=150;
                        
                end % end switch epoch
                
            case 9
                % Day9 - Stim On; Rip Detectn on Tets 3,4,5,6,10,11
                % 150uA; SD=5-6; Nthr=2
                
                nranges=1;
                switch epoch
                    case 2
                        [range1] = timetrans({'00:40:09', '00:55:23'},10000,2); amp1=150;
                    case 4
                        [range1] = timetrans({'01:17:48', '01:32:49'},10000,2); amp1=150;
                        
                end % end switch epoch
                
                
            case 15
                % 15ripdistest2 - Stim On; Rip Detectn on Tet 10
                % 100uA; SD=5/4; Nthr=1
                nranges=1;
                [range1] = timetrans({'00:07:48', '00:13:58'},10000,2); amp1=100;
                
            case 16
                % 16ripdistest3 - Stim On; Rip Detectn on Tet 10
                % 250uA/300uA; SD=5/4; Nthr=1
                nranges=2;
                [range1] = timetrans({'01:02:06', '01:06:06'},10000,2); amp1=250;
                [range2] = timetrans({'01:06:15', '01:10:15'},10000,2); amp2=300;
                
            case 17
                % 17ripdistest4 - Stim On; Rip Detectn on Tet 3 in cell
                % layer
                % 100/150/200/0uA; SD=5/4; Nthr=1
                nranges=4;
                [range1] = timetrans({'00:18:25', '00:27:29'},10000,2); amp1=0;
                [range2] = timetrans({'00:06:29', '00:10:59'},10000,2); amp2=100;
                [range3] = timetrans({'00:11:05', '00:14:20'},10000,2); amp3=150;
                [range4] = timetrans({'00:14:25', '00:18:20'},10000,2); amp4=200;
                
            case 18
                % 17ripdistest4 - Stim On; Rip Detectn on Tet 11 in cell
                % layer -Spikes on Tet10,6,11
                % 20/40/60/80uA/100/0; SD=4.5; Nthr=1
                nranges=6;
                [range1] = timetrans({'00:36:42', '00:40:11'},10000,2); amp1=0;
                [range2] = timetrans({'00:20:29', '00:23:29'},10000,2); amp2=20;
                [range3] = timetrans({'00:23:35', '00:26:40'},10000,2); amp3=40;
                [range4] = timetrans({'00:27:00', '00:30:12'},10000,2); amp4=60;
                [range5] = timetrans({'00:30:24', '00:33:24'},10000,2); amp5=80;
                [range6] = timetrans({'00:33:40', '00:36:40'},10000,2); amp6=100;
                
        end
        
        
    case 'REc'
        
        switch day        % BY DAY WITHIN ANIMAL
            
            case 15
                % 15ripdistest1 - Stim Off; Rip Detectn on 5 Tets
                % 7,8,9,10,4; SD=4; Nthr=1 and Nthr =2;
                nranges=2;
                [range1] = timetrans({'00:21:52', '00:24:52'},10000,2); amp1=1;
                [range2] = timetrans({'00:25:00', '00:28:00'},10000,2); amp2=2;
                
            case 16
                % 16ripdistest3 - Stim Off; Rip Detecn on Tet 10; Nthr=1;
                nranges=1;
                [range1] = timetrans({'00:05:22', '00:08:22'},10000,2); amp1=1; %0
                
            case 18
                % 18lintest1 - Stim Off; Rip Detecn on 5 Tets 7,8,9,10,4;
                % Nthr=2; Lots of cable noise. 5 epochs, use only 1st one
                nranges=1;
                [range1] = timetrans({'00:47:11', '01:02:17'},10000,2); amp1=1; %0
                
            case 19
                % 19ripdistest4 - 900uA Stim ; Rip Detecn on 5 Tets 6,7,8,9,10;
                % Nthr=3; Cable noise, but rat mostly sleeping.
                nranges=1;
                [range1] = timetrans({'00:24:49', '00:27:49'},10000,2); amp1=1; %0
                
            case 20
                % 20ripdistest5 Stim ; Rip Detecn on 4 Tets 6,8,9,10;
                % Nthr=2; Cable noise, but rat mostly sleeping.
                nranges=4;
                [range1] = timetrans({'00:48:15', '00:49:26'},10000,2); amp1=0; %0
                [range2] = timetrans({'00:45:37', '00:48:37'},10000,2); amp2=300;
                [range3] = timetrans({'00:42:20', '00:45:20'},10000,2); amp3=500;
                [range4] = timetrans({'00:39:05', '00:42:05'},10000,2); amp4=700;
                
            case 21
                % 21ripdistest6 Stim ; Rip Detecn on 5 Tets 4,6,8,9,10;
                % Nthr=2; Rat mostly sleeping.
                nranges=4;
                [range1] = timetrans({'01:37:05', '01:40:05'},10000,2); amp1=0; %0
                [range2] = timetrans({'01:34:00', '01:37:00'},10000,2); amp2=1500;
                [range3] = timetrans({'01:30:48', '01:33:48'},10000,2); amp3=3000;
                [range4] = timetrans({'01:27:32', '01:30:32'},10000,2); amp4=6000;
                
            case 01
                % 01_112310 Stim ; Rip Detecn on 4 Tets 6,8,9,10;
                % Nthr=2; Cable noise
                nranges=1;
                switch epoch
                    case 2
                        [range1] = timetrans({'00:24:15', '00:39:16'},10000,2); amp1=900;
                    case 4
                        [range1] = timetrans({'01:02:20', '01:17:20'},10000,2); amp1=900;
                end % end switch epoch
                
            case 02
                % 02_112410 Stim ; Rip Detecn on 4 Tets 6,8,9,10;
                % Nthr=2;
                nranges=1;
                switch epoch
                    case 2
                        [range1] = timetrans({'00:22:30', '00:37:31'},10000,2); amp1=1500; % epoch 2
                    case 4
                        [range1] = timetrans({'01:15:00', '01:30:01'},10000,2); amp1=2000; % epoch 4
                end
                
            case 03
                % 03_112510 Stim ; Rip Detecn on 5 Tets 4,6,8,9,10;
                % Nthr=2;
                nranges=1;
                switch epoch
                    case 2
                        [range1] = timetrans({'00:23:25', '00:38:26'},10000,2); amp1=3000; % epoch 2
                    case 4
                        [range1] = timetrans({'01:01:05', '01:16:06'},10000,2); amp1=3000; % epoch 4
                end
                
            case 04
                % 04_112610 Stim ; Rip Detecn on 5 Tets 4,6,8,9,10;
                % Nthr=2;
                nranges=1;
                switch epoch
                    case 2
                        [range1] = timetrans({'00:40:58', '00:55:59'},10000,2); amp1=4000; % epoch 2
                    case 4
                        [range1] = timetrans({'01:17:25', '01:32:26'},10000,2); amp1=4000; % epoch 4
                end
                
            case 05
                % 05_112710 Stim ; Rip Detecn on 5 Tets 4,6,8,9,10;
                % Nthr=2;
                nranges=1;
                switch epoch
                    case 2
                        [range1] = timetrans({'00:21:53', '00:36:54'},10000,2); amp1=6000; % epoch 2
                    case 4
                        [range1] = timetrans({'01:01:12', '01:16:13'},10000,2); amp1=6000; % epoch 4
                end
                
            case 06
                % 06_112810 Stim ; Rip Detecn on 5 Tets 4,6,8,9,10;
                % Nthr=2;
                nranges=1;
                switch epoch
                    case 2
                        [range1] = timetrans({'00:34:40', '00:49:41'},10000,2); amp1=5000; % epoch 2
                    case 4
                        [range1] = timetrans({'01:11:07', '01:26:08'},10000,2); amp1=5000; % epoch 4
                end
                
            case 07
                % 07_112910 Stim ; Rip Detecn on 5 Tets 4,6,8,9,10;
                % Nthr=2;
                nranges=1;
                switch epoch
                    case 2
                        [range1] = timetrans({'00:20:31', '00:35:32'},10000,2); amp1=5000; % epoch 2
                    case 4
                        [range1] = timetrans({'00:54:06', '01:09:07'},10000,2); amp1=5000; % epoch 4
                end
                
        end % end switch day
        
        
    case 'REb'
        
        switch day        % BY DAY WITHIN ANIMAL
            
            case 16
                % 16ripdistest1
                nranges=7; % 0uA to 120 uA
                [range1] = timetrans({'00:20:45', '00:22:50'},10000,2); amp1=0; %0
                [range2] = timetrans({'00:39:00', '00:42:00'},10000,2); amp2=20; %20
                [range3] = timetrans({'00:35:50', '00:38:50'},10000,2); amp3=40;%40
                [range4] = timetrans({'00:32:40', '00:35:40'},10000,2); amp4=60;%60
                [range5] = timetrans({'00:29:29', '00:32:29'},10000,2); amp5=80;%80
                [range6] = timetrans({'00:22:50', '00:25:50'},10000,2); amp6=100;%100
                [range7] = timetrans({'00:25:58', '00:28:58'},10000,2); amp7=120;%120
                
            case 17
                % 17riptest
                nranges=4;
                [range1] = timetrans({'00:01:37', '00:02:40'},10000,2); amp1=0;
                [range2] = timetrans({'00:02:48', '00:04:02'},10000,2); amp2=60;
                [range3] = timetrans({'00:04:14', '00:05:14'},10000,2); amp3=80;
                [range4] = timetrans({'00:05:23', '00:06:33'},10000,2); amp4=100;
                
            case 19
                %19ripdistest2
                nranges=5; % 0uA, 80, 100, 130 150
                [range1] = timetrans({'00:50:51', '00:52:51'},10000,2); amp1=0;
                [range2] = timetrans({'00:52:51', '00:54:51'},10000,2); amp2=80;
                [range3] = timetrans({'00:44:24', '00:46:24'},10000,2); amp3=100;
                [range4] = timetrans({'00:46:40', '00:48:40'},10000,2); amp4=130;
                [range5] = timetrans({'00:48:50', '00:50:50'},10000,2); amp5=150;
                
            case 20
                % 20ripdistest3
                nranges=5; % 0uA, 150, 300, 500 700
                [range1] = timetrans({'00:11:39', '00:13:45'},10000,2); amp1=0; %0
                [range2] = timetrans({'00:09:35', '00:11:35'},10000,2); amp2=150;%150
                [range3] = timetrans({'00:07:16', '00:09:16'},10000,2); amp3=300;%300
                [range4] = timetrans({'00:02:35', '00:04:35'},10000,2); amp4=500;%500
                [range5] = timetrans({'00:04:55', '00:06:56'},10000,2); amp5=700;%700
                
            case 21
                % 21ripdistest4
                nranges=4; % 0uA, 150, 300, 500
                [range1] = timetrans({'00:08:04', '00:10:20'},10000,2); amp1=0;%0
                [range2] = timetrans({'00:01:26', '00:03:30'},10000,2); amp2=150;%150
                [range3] = timetrans({'00:03:42', '00:05:46'},10000,2); amp3=300;%300
                [range4] = timetrans({'00:05:58', '00:08:02'},10000,2); amp4=500;%500
                
        end % end switch day
        
    case 'RCa'
        
        switch day
            
            case 20
                % 20Riptesttime
                nranges=4;
                [range1] = timetrans({'00:01:50', '00:03:50'},10000,2);
                [range2] = timetrans({'00:03:52', '00:05:52'},10000,2);
                [range3] = timetrans({'00:05:55', '00:07:55'},10000,2);
                [range4] = timetrans({'00:07:57', '00:10:00'},10000,2);
                
            case 21
                % 21Ripdis
                nranges=1;
                [range1] = timetrans({'00:02:05', '00:07:10'},10000,2); amp1=250; % 300us
                
            case 23
                %23RCaRipamptest_051910
                nranges=5;
                [range1] = timetrans({'00:29:25', '00:31:25'},10000,2);
                [range2] = timetrans({'00:31:26', '00:33:31'},10000,2);
                [range3] = timetrans({'00:33:42', '00:35:42'},10000,2);
                [range4] = timetrans({'00:35:50', '00:38:05'},10000,2);
                [range5] = timetrans({'00:38:13', '00:40:13'},10000,2);
                
            case 26
                %26RCaRipamptestn_051910
                nranges=7;
                [range1] = timetrans({'00:07:42', '00:09:42'},10000,2);
                [range2] = timetrans({'00:09:44', '00:11:44'},10000,2);
                [range3] = timetrans({'00:11:54', '00:13:57'},10000,2);
                [range4] = timetrans({'00:15:03', '00:17:03'},10000,2);
                [range5] = timetrans({'00:17:11', '00:19:11'},10000,2);
                [range6] = timetrans({'00:19:18', '00:21:18'},10000,2);
                [range7] = timetrans({'00:21:25', '00:23:30'},10000,2);
                
            case 27
                %27RCaRipamptestn_052110
                nranges=7;
                [range1] = timetrans({'00:03:25', '00:05:25'},10000,2);
                [range2] = timetrans({'00:05:27', '00:07:27'},10000,2);
                [range3] = timetrans({'00:07:28', '00:09:28'},10000,2);
                [range4] = timetrans({'00:09:30', '00:11:35'},10000,2);
                [range6] = timetrans({'00:11:45', '00:13:45'},10000,2);
                [range5] = timetrans({'00:13:50', '00:15:50'},10000,2);
                [range7] = timetrans({'00:15:58', '00:18:09'},10000,2);
                
        end
        
end % end switch prefix

%% Divide by Stimulation amplitude - Arrange in matrix form and bar graph form

% Each stimulus start time is converted to 0.1ms resolution to get index right
stim_starttime_comp = stim_starttime*10; % Convert from ms to 10000 pts in sec

Stimhist_matr=[];
for n=1:nranges
    
    cmd=sprintf('startu=range%d(1);',n); eval(cmd);
    cmd=sprintf('endu=range%d(2);',n); eval(cmd);
    % Index range by comparison of stim times to range start and end
    currrange = find(stim_starttime_comp>=startu & stim_starttime_comp<=endu);
    
    temp = stim_spkshist( currrange,:);
    temp_matrix= stim_spkshist( currrange,:); % Separate for matrix to optimize color plot
    
    tempeeg = e_stim(currrange,:);
    tempripamp = ripamp_stim(currrange,:);
    tempripenv = ripenv_stim(currrange,:);
    
   if (nranges==1) || ( nranges>1 && n~=1) % To avoid removing stim artifact for 0 amplitude
 
        if binsize==10
            temp_matrix(:,(pret/binsize)) = 0;
            temp_matrix(:,(pret/binsize)+1) = 0;
            temp_matrix(:,(pret/binsize)+2) = 0;
            
            temp(:,(pret/binsize)) = 0;
            temp(:,(pret/binsize)+1) = 0;
            temp(:,(pret/binsize)+2) = 0;
            
        end
        
        if binsize==5
            
            %             temp(:,(pret/binsize)-2) = 0;
            %             temp(:,(pret/binsize)-1) = 0;
            %             temp(:,(pret/binsize)) = 0;
            %             temp(:,(pret/binsize)+1) = temp(:,(pret/binsize)+8);
            %             temp(:,(pret/binsize)+2) = temp(:,(pret/binsize)+8);
            %             temp(:,(pret/binsize)+3) = temp(:,(pret/binsize)+5);
            %             temp(:,(pret/binsize)+3) = 0;
            
            temp_matrix(:,(pret/binsize)) = 0;
            temp_matrix(:,(pret/binsize)+1) = 0;
            temp_matrix(:,(pret/binsize)+2) = 0;
            temp_matrix(:,(pret/binsize)+3) = 0;
            
            temp(:,(pret/binsize)) = 0;
            temp(:,(pret/binsize)+1) = 0;
            temp(:,(pret/binsize)+2) = 0;
            temp(:,(pret/binsize)+3) = 0;
            
            
        end
    end
    
    
    
    cmd=sprintf('stimhistall%d = temp;',n); eval(cmd);
    cmd=sprintf('stimhist%d = mean(temp,1);',n); eval(cmd);
    
    cmd=sprintf('eeghistall%d = tempeeg;',n); eval(cmd);
    cmd=sprintf('ripamphistall%d = tempripamp;',n); eval(cmd);
    cmd=sprintf('ripenvhistall%d = tempripenv;',n); eval(cmd);
    
    Stimhist_matr = [Stimhist_matr; mean(temp_matrix,1)];
    
end


%% ------------------------------------------------
% PLOT
% -------------------------------------------------


%% Plot Ripple Power for all stimulations on one plot

% ----------------------------------------------
figure; hold on;
%redimscreen_halfvert(0);
orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
for n=1:nranges
    cmd=sprintf('ripenv_hist = ripenvhistall%d;',n); eval(cmd);
    cmd=sprintf('amp = amp%d;',n); eval(cmd);
    subplot(nranges,1,(nranges-n+1)); hold on;
    yplot = mean(ripenv_hist,1);
    
    taxis = [1:size(ripenv_stim,2)]*1000/e.samprate;
    taxis = taxis-pret;
    %taxis = taxis - 0.2*e.samprate;
    %set(gca,'XLim',[-pret-binsize postt+binsize]);
    
    %jbfill(taxis,yplot+sem(ripenv_hist,1),yplot-sem(ripenv_hist,1));
    %plot(taxis, yplot+sem(ripenv_hist,1),'r','Linewidth',2);
    %plot(taxis, yplot-sem(ripenv_hist,1),'r','Linewidth',2);
    
    % Get Axis if 0uA and get Baseline Fr0 from (pret/2)/binsize bins, eg. 100ms if pret is 200ms
    if n==1
        [ylimits] = get(gca,'YLim');
        y1=ylimits(1); y2=ylimits(2);
        basefr0 = mean(yplot(1:round(0.1*e.samprate)));
        basefr0 = roundn(basefr0,-1);
        std0 = std(yplot(1:round(0.1*e.samprate))); std0 = roundn(std0,-1);
        plot(taxis, basefr0*ones(size(taxis)),'k--','Linewidth',2 );
        plot(taxis, (basefr0+5*std0)*ones(size(taxis)),'c-.','Linewidth',2 );
        text(-pret+0.05*e.samprate,0.97*y2,['Power ' num2str(basefr0)], 'FontSize',14,'Fontweight','bold')
    end
    
    
    if n~=1,
        curr_basefr = mean(yplot(1:0.1*e.samprate));
        curr_basefr = roundn(curr_basefr,-1);
        cstd = std(yplot(1:0.1*e.samprate)); cstd = roundn(cstd,-1);
        plot(taxis, basefr0*ones(size(taxis)),'k--','Linewidth',2 );
        plot(taxis, (basefr0+5*std0)*ones(size(taxis)),'c-.','Linewidth',2 );
        plot(taxis, curr_basefr*ones(size(taxis)),'b--','Linewidth',2 );
        plot(taxis, (curr_basefr+5*cstd)*ones(size(taxis)),'g-.','Linewidth',2 );
        text(-pret+0.05*e.samprate,0.97*y2,['Power ' num2str(curr_basefr)],...
            'FontSize',14,'Fontweight','bold','Color','b');
    end
    
    
   plot(taxis, yplot,'r','Linewidth',4);
     
    %jbfill(taxis,yplot+std(ripenv_hist,1),yplot-std(ripenv_hist,1));
    
    
    % Plot Line at 0ms
    [ylimits] = get(gca,'YLim');
    ypts = ylimits(1):0.1:ylimits(2);
    xpts = 0*ones(size(ypts));
    plot(xpts , ypts, 'k-.','Linewidth',2);
    % Plot lines at 100ms and 200ms
    xpts = 100*ones(size(ypts));
    plot(xpts , ypts, 'k:','Linewidth',2);
    xpts = 200*ones(size(ypts));
    plot(xpts , ypts, 'k:','Linewidth',2);
    % Set Yaxis based on 0uA
    %set(gca,'YLim',[y1 y2]);   % Same axis - All Graphs, based on 0uA
    
    
    % Set Ticks
    %set(gca,'ytick',[]);
    if n~=1,
        set(gca,'xtick',[]);
%     else
%         set(gca,'xtick',[-pret]);
    end
    % Set Labels
    ylabel([num2str(amp) ' uA'],'FontSize',16,'Fontweight','bold');
    if n==1
        %ylabel('MU FR');
        xlabel('Time(ms)');
    end
    % Set Title
    if n==nranges,
        title(['Ripple Power aligned to stimulation'  ' - Tets ' num2str(tetu)],...
            'FontSize',16,'Fontweight','bold');
    end
    % Save Graph if asked to
    if saveg1==1,
        orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
        saveas(gcf,['RipPower_StimAmpl_' num2str(day) '_tet' num2str(tet) '_SEPstd' num2str(std)],'fig');
        saveas(gcf,['RipPower_StimAmpl_' num2str(day) '_tet' num2str(tet) '_SEPstd' num2str(std)],'jpg');
    end
    
end












%% Matrix
% -------
Stimhist_matr=Stimhist_matr(:,1:end-1);
figure; hold on;
orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
imagesc(Stimhist_matr);
title(['Multiunit Firing aligned to stimulation-' num2str(binsize) 'ms bins'],...
    'FontSize',16,'Fontweight','bold');
ylabel('Stimulation Amplitude (uA)','FontSize',16,'Fontweight','bold');
xlabel('Time(ms)','FontSize',16,'Fontweight','bold');
axis([0 ((pret+postt)/binsize)+2 0.5 nranges+0.5]);
ypts = 0:1:nranges+1;
xpts = ((pret/binsize))*ones(size(ypts));
% Plot Line at 0 ms - Onset of stimulation
plot(xpts , ypts, 'k--','Linewidth',2);
% Plot lines at 100ms and 200ms
xpts = ((pret+100)/binsize)*ones(size(ypts));
plot(xpts , ypts, 'k:','Linewidth',3);
xpts = ((pret+200)/binsize)*ones(size(ypts));
plot(xpts , ypts, 'k:','Linewidth',3);

set(gca,'xtick',[-2,(pret/binsize)-1,(pret+postt)/binsize-2],'xticklabel',{num2str([-pret,0,postt]')},...
    'FontSize',14,'Fontweight','bold');

% Make Ylabel
for n=1:nranges,
    cmd=sprintf('amp = amp%d;',n); eval(cmd);
    String{n}=num2str(amp);
    
end

set(gca,'ytick',[1:nranges],'yticklabel',String,...
    'FontSize',14,'Fontweight','bold');

if saveg1==1,
    orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
    saveas(gcf,['Image_MultiFR_StimAmpl_' num2str(day) '_tet' num2str(tet) '_SEPstd' num2str(std)],'fig');
    saveas(gcf,['Image_MultiFrStimAmpl_' num2str(day) '_tet' num2str(tet) '_SEPstd' num2str(std)],'jpg');
end



%% Plot Bar Graphs for all stimulations on one plot
% ----------------------------------------------
figure; hold on;
%redimscreen_halfvert(0);
orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
for n=1:nranges
    cmd=sprintf('stimhist = stimhist%d;',n); eval(cmd);
    cmd=sprintf('amp = amp%d;',n); eval(cmd);
    subplot(nranges,1,(nranges-n+1)); hold on;
    yplot = (1000/binsize_plot)*stimhist; %Multiunit fr in "binsize"ms bins
    taxis = [-pret:binsize:postt];
    % taxis = taxis - pret/binsize;
    set(gca,'XLim',[-pret-binsize postt+binsize]);
    bar(taxis, yplot,'r');
    % Get Axis if 0uA and get Baseline Fr0 from (pret/2)/binsize bins, eg. 100ms if pret is 200ms
    if n==1
        [ylimits] = get(gca,'YLim');
        y1=ylimits(1); y2=ylimits(2);
        basefr0 = mean(yplot(1:(pret/2)/binsize));
        basefr0 = round(10*basefr0)./10; % Because roundn does not exist in Lab matlab
        std0 = std(yplot(1:(pret/2)/binsize)); std0 = roundn(std0,-1);
        plot(taxis, (basefr0+3*std0)*ones(size(taxis)),'c-.','Linewidth',2 );
        text(-pret+binsize,0.8*y2,['Fir Rate ' num2str(basefr0) 'Hz'], 'FontSize',14,'Fontweight','bold')
    end
    % Plot Line at 0ms
    ypts = 0:1:y2;
    xpts = 0*ones(size(ypts));
    plot(xpts , ypts, 'k-.','Linewidth',2);
    
    % Plot lines at 100ms and 200ms
    xpts = 100*ones(size(ypts));
    plot(xpts , ypts, 'k:','Linewidth',3);
    xpts = 200*ones(size(ypts));
    plot(xpts , ypts, 'k:','Linewidth',3);
    % Set Yaxis based on 0uA
    set(gca,'YLim',[y1 y2]);   % Same axis - All Graphs, based on 0uA
    % Plot Baseline FR
    plot(taxis, basefr0*ones(size(taxis)),'k--','Linewidth',2 );
    if n~=1,
        curr_basefr = mean(yplot(1:(pret/2)/binsize));
        curr_basefr = round(10*curr_basefr)./10; % Because roundn does not exist in Lab matlab
        cstd = std(yplot(1:(pret/2)/binsize)); cstd = roundn(cstd,-1);
        plot(taxis, curr_basefr*ones(size(taxis)),'b--','Linewidth',2 );
        %plot(taxis, (curr_basefr+3*cstd)*ones(size(taxis)),'g-.','Linewidth',2 );
        text(-pret+binsize,0.8*y2,['Fir Rate ' num2str(curr_basefr) 'Hz'],...
            'FontSize',14,'Fontweight','bold','Color','b');
    end
    % Set Ticks
    set(gca,'ytick',[]);
    if n~=1,
        set(gca,'xtick',[]);
    end
    % Set Labels
    ylabel([num2str(amp) ' uA'],'FontSize',16,'Fontweight','bold');
    if n==1
        %ylabel('MU FR');
        xlabel('Time(ms)');
    end
    % Set Title
    if n==nranges,
        title(['Multiunit Firing aligned to stimulation in ' num2str(binsize) 'ms bins - Tets ' num2str(tetu)],...
            'FontSize',16,'Fontweight','bold');
    end
    % Save Graph if asked to
    if saveg1==1,
        orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
        saveas(gcf,['MultiFR_StimAmpl_' num2str(day) '_tet' num2str(tet) '_SEPstd' num2str(std)],'fig');
        saveas(gcf,['MultiFrStimAmpl_' num2str(day) '_tet' num2str(tet) '_SEPstd' num2str(std)],'jpg');
    end
    
end


%% Plot Bar graph for 0uA and 2nd chosen amplitude
% --------------------------------------------------
figure; hold on;
orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
redimscreen_widevert(0);

[ns]=[1,twoplot_idx]; % 0uA and 2nd amplitude
for i=1:length(ns),
    n=ns(i);
    cmd=sprintf('stimhist = stimhist%d;',n); eval(cmd);
    cmd=sprintf('amp = amp%d;',n); eval(cmd);
    subplot(1,length(ns),i); hold on;
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
        text(-pret+binsize,0.8*y2,['Fir Rate ' num2str(basefr0) 'Hz'], 'FontSize',14,'Fontweight','bold');
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
    plot(taxis, (basefr0+3*std0)*ones(size(taxis)),'c-.','Linewidth',2 );
    if n~=1,
        curr_basefr = mean(yplot(1:(pret/2)/binsize));
        curr_basefr = round(10*curr_basefr)./10; % Because roundn does not exist in Lab matlab
        cstd=std(yplot(1:(pret/2)/binsize));
        cstd=roundn(cstd,-1);
        plot(taxis, curr_basefr*ones(size(taxis)),'b--','Linewidth',2 );
        plot(taxis, (curr_basefr+3*cstd)*ones(size(taxis)),'g-.','Linewidth',2 );
        text(-pret+binsize,0.8*y2,['Fir Rate ' num2str(curr_basefr) 'Hz'],...
            'FontSize',14,'Fontweight','bold','Color','b');
    end
    
    
    %set(gca,'xtick',[]);

    ylabel([ num2str(amp) 'uA'],'FontSize',16,'Fontweight','bold');
    if i==length(ns)
        xlabel('Time(ms)','FontSize',16,'Fontweight','bold');
    end
    if i==1
        title(['Stimulation Calibration- ' num2str(binsize) 'ms bins; Tet ' num2str(tetu)],'FontSize',14,'Fontweight','bold');
    end
end


%% Plot example single-trial firing rate and ripple power aligned to stimulation
% ------------------------------------------------------------------------

if plot_ex==1
    
    figure; hold on;
    orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
    redimscreen_widevert(0);
    
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
                text(-pret+binsize,0.8*y2,['Fir Rate ' num2str(basefr0) 'Hz'], 'FontSize',10,'Fontweight','bold');
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
                    'FontSize',10,'Fontweight','bold','Color','b');
            end
            
            
            %set(gca,'xtick',[]);
            if p==1
                title([num2str(amp) 'uA; ' num2str(binsize) 'ms bins; Tet ' num2str(tetu)],'FontSize',14,'Fontweight','bold');
                %ylabel([ num2str(amp) 'uA'],'FontSize',14,'Fontweight','bold');
            end
            
            if p==5,
                xlabel('Time(ms)','FontSize',14,'Fontweight','bold');
            end
            
            %     if i==1
            %         title(['Stimulation Calibration'],'FontSize',12,'Fontweight','bold');
            %     end
            
            
        end
        
    end
    
end


%% Plot Multiunit rate around ripples
% -------------------------------------

if dorip==1
    %plot(taxis, 2*ripnostim_pre(i,:),['r-'],'Linewidth',2,'Markersize',6);
    
    ypts = 0:1.1*max(yplot);
    xpts = 200*ones(size(ypts));
    plot(xpts , ypts, 'k--','Linewidth',2);
    
    title(['Sleep Multiunit Firing around ripples-Day' num2str(day)],...
        'FontSize',24,'Fontweight','bold');
    %axis([0 800 -800 600]);
    ylabel('Instantaeous Multiunit Firing Rate');
    xlabel('Time(ms)');
    
    %text( 4, 2450,['DetRate(4,6,7) =' num2str(round(DetRateper*100)/100) '%'],'FontSize', 24, 'FontWeight','bold');
    
    if saveg1==1,
        orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
        saveas(gcf,['MultiFRrip_' num2str(day) '_tet' num2str(tet) '_SEPstd' num2str(sd)],'fig');
        saveas(gcf,['MultiFrrip_' num2str(day) '_tet' num2str(tet) '_SEPstd' num2str(sd)],'jpg');
    end
    
end



