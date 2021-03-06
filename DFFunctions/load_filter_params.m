
function Fp = load_filter_params(Fp, varargin)
% Fp = load_filter_params(Fp, varargin)
% Fp.params : cell array of strings
% each string matches a param set
% succesors take priority if vars clash

add_paths = 1;
if ~isempty(varargin)
    assign(varargin{:})
end

epochfilter = '';
tetfilter = '';
cellFilter = '';
timefilter = {};
options = {};
sysDateTime = clock;
fprintf('----filter params----\n');
for s = Fp.params
fprintf('* %s\n', s{1})
switch s{1}
%% ========= XCnormAC ========= 
case 'rmsXCnormAC'

case 'dfa_rmsXCnormAC'
    iterator = 'singleDayAnal';
    filtfunction = 'dfa_rmsXCnormAC';
    datatypes = {'ca1rippleskons', 'lick', 'DIO', 'task', 'pos'};
    options = {};
    
case 'ccPvRXCfromRew'

case 'dfa_ccPvRXCfromRew'    
    iterator = 'singleDayAnal';
    filtfunction = 'dfa_ccPvRXCfromRew';
    datatypes = {'ca1rippleskons', 'lick', 'DIO', 'task'};
    options = {};
    
case 'XCfromRew'    
    win = [-1 1];
    bin = .01;
    
case 'dfa_XCfromRew'
    iterator = 'singleDayAnal';
    filtfunction = 'dfa_XCfromRew';
    datatypes = {'ca1rippleskons', 'lick', 'DIO', 'task'};
    options = {};
    

%% ========= XPtrigSWR =========    
case 'XPtrigSWR'
    win = [-1 1];
    bin = .001;
    eventType = 'ca1rippleskons';    
    
case 'dfa_XPtrigSWR'    
    iterator = 'singleDayAnal'; 
    filtfunction = 'dfa_XPtrigSWR';
    datatypes = {'ca1rippleskons', 'lick'};
    options = {'win', win, 'bin', bin, 'eventType', eventType}; 
    
%% ========= rewTrigSWRXP =========
case 'rewTrigSWRXP'
    win = [-2 20];
    bin = .001;
    eventType = 'ca1rippleskons';
    minILIthresh = .01;
    
case 'dfa_rewTrigSWRXP'
    iterator = 'singleDayAnal';
    filtfunction = 'dfa_rewTrigSWRXP';
    datatypes = {'ca1rippleskons', 'lick', 'DIO', 'task'};
    options = {'win', win, 'bin', bin, 'eventType', eventType, ...
        'minILIthresh', minILIthresh};     
    
%% ========= common =========    
case 'valid_ntrodes'
    %             ntAreas = {'ca1', 'mec', 'ref'}; %, 'por', 'v2l', 'sub'};
    areas = {{'ca1', 'd'}, {'mec', 'deep'}, {'mec', 'supf'}};
    tetfilter = '(isequal($valid,''yes'') && (isequal($area,''ca1'') || isequal($area,''mec'')))';
%% ========= XP pairwise phase corr =========    
case 'dfa_phaseXcorr'
    bin = 0.01; % fraction of full rotations
    tmax = 2; % full rotations
    
    iterator = 'iter_multiEpoch_multiCluster';
    filtfunction = 'dfa_phaseXcorr';

    cellpairfilter = {'allcomb', ...
        '($numspikes > 100) && (isequal($area, ''ca1'')) && (all(cellfun(''isempty'',(arrayfun(@(x) strfind(x,''mua''), $tags, ''un'', 0)))))', ...
        '($numspikes > 100) && (isequal($area, ''ca1'')) && (all(cellfun(''isempty'',(arrayfun(@(x) strfind(x,''mua''), $tags, ''un'', 0)))))'};

    eventType = 'lick';
    datatypes = {'spikes', 'lick'};
    options = {'eventName',eventType, 'bin', bin, 'tmax', tmax};
%% ========= XP-SWR mod REACTIVATION =========
case 'dfa_reactivationPSTH'
    bin = 0.1; % seconds
    win = [-2 2];
    perPC = 1;
    % cellfilter
    cellFilter = '(isequal($area, ''ca1'')) && ($numspikes > 100) && (all(cellfun(''isempty'',(arrayfun(@(x) strfind(x,''mua''), $tags, ''un'', 0)))))';
    % ripple filter
%     eventType = 'ca1rippleskons';
    consensus_numtets = 2;   % minimum # of tets for consensus event detection
    minstdthresh = 2;        % STD. how big your ripples are
    exclusion_dur = .5;  % seconds within which consecutive events are eliminated / ignored
    minvelocity = 0;
    maxvelocity = 4;
    rippleFilter = {{'getconstimes', '($cons == 1)', ...
        'ca1rippleskons', 1,'consensus_numtets',consensus_numtets, ...
        'minstdthresh', minstdthresh,'exclusion_dur',exclusion_dur, ...
        'minvelocity', minvelocity,'maxvelocity',maxvelocity}};    
    % lick burst filter
    maxIntraBurstILI = 0.25;  % max burst ili threshold in seconds
    minBoutLicks = 3; % filter out bouts with less than boutNum licks
    minLickTimeFromSwr = 1; % time duration from closest lick
    
    % template filter
%     templateFilter = {'get2dstate','($velocity>4)'};
      
    iterator = 'singleDayAnal'; % city
    datatypes = {'cellinfo'};
    options = {'win', win, 'bin', bin, 'perPC', perPC, 'maxIntraBurstILI', maxIntraBurstILI, ...
        'minBoutLicks', minBoutLicks, 'cellfilter', cellFilter, ...
        'rippleFilter', rippleFilter};

%% ========= SWR iLB vs eLB =========
case 'pctILB'
    eventType = 'ca1rippleskons';
        % lick burst filter
    maxILIthresh = 1; % max burst ili threshold in seconds
    minBoutLicks = 2; % filter out bouts with less than boutNum licks
    minILIthresh = .06; 
    stdthresh = [2 3 4];
    
case 'dfa_pctILB'   
    iterator = 'singleDayAnal'; %'singleepochanal';
    filtfunction = 'dfa_pctILB';
    datatypes = {'ca1rippleskons','task', 'lick', 'DIO'};
    options = {'eventType', eventType, 'maxILIthresh', maxILIthresh,...
        'minBoutLicks', minBoutLicks, 'stdthresh', stdthresh, 'minILIthresh', minILIthresh};
        
%% ========= XP-SWR mod =========
case 'XPtrigAvgRip'
    win = [-2 2];
    eventType = 'ca1rippleskons';
    swin = 1500;
    stp = 375;
    wienerFiltDiv = 10000;
    circperms = 25;
    
case 'dfa_XPtrigAvgRip'
    iterator = 'singleDayAnal'; %'singleepochanal';
    filtfunction = 'dfa_XPtrigAvgRip';
    datatypes = {'ca1rippleskons', 'lick', 'pos'};
    options = {'win', win, 'eventType', eventType, 'swin', swin, 'stp', stp, ...
        'wienerFiltDiv', wienerFiltDiv, 'circperms', circperms};     
    

case 'wXPTrigSWR'
    %     maxTimeSinceRew = 5;
    bin = .01;
    tmax = .5;
    eventType = 'ca1rippleskons';
    minILIthresh = .080; % seconds
    maxILIthresh = .250; % seconds
    minBoutLicks = 3;
    % xcorr / excorr
    excShortBin = bin; 
    excLongBin = .500;
    rmsmincounts = 1; % min bin count within rmstamax. otherwise nan
    rmstmax = .25; % seconds
    % shuf
    compute_shuffle = 1;
    numshuffs = 1000;
    maxShift = tmax/2; %ms integer. timemod shift


case 'dfa_lickswrcorr'
    % make sure to include a 'ripples' timefilter in the paramset
    % func
    iterator = 'singleDayAnal'; %'singleepochanal';
    filtfunction = 'dfa_lickswrcorr';
    datatypes = {'ca1rippleskons','task', 'lick', 'DIO'};
    options = {'bin', bin, 'tmax', tmax, 'eventType', eventType, ...
        'excShortBin', excShortBin, 'excLongBin', excLongBin, ...
        'rmsmincounts', rmsmincounts, 'rmstmax', rmstmax, 'compute_shuffle', ...
        compute_shuffle, 'numshuffs', numshuffs, 'shuffOffset', maxShift, ...
        'minILIthresh', minILIthresh, 'maxILIthresh', maxILIthresh, ...
        'minBoutLicks', minBoutLicks};         
    
case 'XPprepostSWR'    
        %     maxTimeSinceRew = 5;
    bin = .01;
    tmax = .5;
    eventType = 'ca1rippleskons';
    minILIthresh = .06; % seconds
    maxILIthresh = 1; % seconds
    minBoutLicks = 2;
    % xcorr / excorr
    excShortBin = bin; 
    excLongBin = .500;
    rmsmincounts = 1; % min bin count within rmstamax. otherwise nan
    rmstmax = .25; % seconds
    % shuf
    compute_shuffle = 1;
    numshuffs = 1000;
    maxShift = tmax/2; %ms integer. timemod shift
    
case 'dfa_XPprepostSWR'
    % make sure to include a 'ripples' timefilter in the paramset
    % func
    iterator = 'singleDayAnal'; %'singleepochanal';
    filtfunction = 'dfa_XPprepostSWR';
    datatypes = {'ca1rippleskons','task', 'lick', 'DIO'};
    options = {'bin', bin, 'tmax', tmax, 'eventType', eventType, ...
        'excShortBin', excShortBin, 'excLongBin', excLongBin, ...
        'minILIthresh', minILIthresh, 'maxILIthresh', maxILIthresh, ...
        'rmsmincounts', rmsmincounts, 'rmstmax', rmstmax, 'compute_shuffle', ...
        compute_shuffle, 'numshuffs', numshuffs, 'shuffOffset', maxShift, ...
        'minBoutLicks', minBoutLicks};     
    
  
%% ========= 'dfa_eventTrigSpiking'=========
case 'lickboutlicks'
%     minILIthresh = .06; % seconds
%     maxILIthresh = .250; % seconds
%     minBoutLicks = 3;
%     timefilter{end+1} = {'getLickBoutLicks', '($LBlick == 1)', ...
%         'minILIthresh', minILIthresh, 'maxILIthresh', maxILIthresh, ...
%         'minBoutLicks', minBoutLicks};

case 'dfa_eventTrigSpiking'           
    iterator = 'singleDayCellAnal';
    filtfunction = 'dfa_eventTrigSpiking';
    datatypes = {'spikes', 'cellinfo', eventType};
    options = {'win', win, 'bin', bin, 'wbin', wbin, ...
        'eventType', eventType};

case 'wtrackSWRTrigSpiking'
    % for dfa_eventTrigSpiking
    win = [1.5 1.5];
    bin = .001;
    wbin = .01; % seconds. wider psth
%             smbins = 10; % bins. smooth across x bins (wbin x smbins = range of influence)
    eventType = 'ca1rippleskons';
    
case 'PeventSpiking'% formely wtrackLickTrigSpiking
    win = [1.5 1.5];
    bin = .001;
    wbin = .01; % seconds. wider psth
%             smbins = 10; % bins. smooth across x bins (wbin x smbins = range of influence)
    eventType = 'Pevent';
    expvars = {'all'}; %, 'wetLickBursts', 'dryLickBursts'};

case 'wtrackLickTrigSpiking' % formerly wtrackLickTrigSpiking
    % for dfa_eventTrigSpiking
    win = [1.5 1.5];
    bin = .001;
    wbin = .01; % seconds. wider psth
%             smbins = 10; % bins. smooth across x bins (wbin x smbins = range of influence)
    eventType = 'lick';
    expvars = {'all'}; %, 'wetLickBursts', 'dryLickBursts'};

%% ========= 'dfa_eventTrigLFP'=========
case 'dfa_eventTrigLFP'
    iterator = 'multitetrodeanal';
    filtfunction = 'dfa_eventTrigLFP';
    datatypes = {eventType, LFPtypes{:}};
    options = {'eventType', eventType, 'LFPtypes', LFPtypes, 'win', win};

case 'wtrackSWRTrigLFP'
    % for dfa_eventTrigLFP
    eventType = 'ca1rippleskons';
    win = [1.5 1.5];
    LFPtypes = {'eeg'};

case 'wtrackLickTrigLFP'
    % for dfa_eventTrigLFP
    eventType = 'lick';
    win = [1.5 1.5];
    LFPtypes = {'eeg'};

%% ========= ripPos    
case 'ripPos'    
    eventType = 'ca1rippleskons';
    
case 'dfa_ripPos'
    iterator = 'singleDayAnal'; %'singleepochanal';
    filtfunction = 'dfa_ripPos';
    datatypes = {'ca1rippleskons','pos', 'lick'};
    options = {'eventType', eventType};     
    
    
%% ========= DAY EPOCH TETRODE CELL RIPPLE filters
case 'day1'
    days = 1;
    
case 'ripples>3'
    eventType = 'ca1rippleskons';
    consensus_numtets = 2;   % minimum # of tets for consensus event detection
    minstdthresh = 3;        % STD. how big your ripples are
    exclusion_dur = 0;  % seconds within which consecutive events are eliminated / ignored
    minvelocity = 0;
    maxvelocity = 4;
    timefilter{end+1} = {'getconstimes', '($cons == 1)', ...
        'ca1rippleskons', 1,'consensus_numtets',consensus_numtets, ...
        'minstdthresh', minstdthresh,'exclusion_dur',exclusion_dur, ...
        'minvelocity', minvelocity,'maxvelocity',maxvelocity};
    
case 'ripples>5'
    eventType = 'ca1rippleskons';
    consensus_numtets = 2;   % minimum # of tets for consensus event detection
    minstdthresh = 5;        % STD. how big your ripples are
    exclusion_dur = 0;  % seconds within which consecutive events are eliminated / ignored
    minvelocity = 0;
    maxvelocity = 4;
    timefilter{end+1} = {'getconstimes', '($cons == 1)', ...
        'ca1rippleskons', 1,'consensus_numtets',consensus_numtets, ...
        'minstdthresh', minstdthresh,'exclusion_dur',exclusion_dur, ...
        'minvelocity', minvelocity,'maxvelocity',maxvelocity};    

    
case 'ripples>2'
    eventType = 'ca1rippleskons';
    consensus_numtets = 2;   % minimum # of tets for consensus event detection
    minstdthresh = 2;        % STD. how big your ripples are
    exclusion_dur = 0;  % seconds within which consecutive events are eliminated / ignored
    minvelocity = 0;
    maxvelocity = 4;
    timefilter{end+1} = {'getconstimes', '($cons == 1)', ...
        'ca1rippleskons', 1,'consensus_numtets',consensus_numtets, ...
        'minstdthresh', minstdthresh,'exclusion_dur',exclusion_dur, ...
        'minvelocity', minvelocity,'maxvelocity',maxvelocity};    

case 'wetLickBursts'
    timefilter{end+1} = {'getWetLickBouts', ...
        '($wetLB == 1)', 'minILIthresh', .01, 'maxILIthresh', 1, ...
        'minBoutLicks', 2};
    
case 'dryLickBursts'
    timefilter{end+1} = {'getWetLickBouts', ...
        '($dryLB == 1)'};   
    
case 'sleep'
    env = 'sleep';
    epochfilter = sprintf('(isequal($environment,''%s''))', ...
        env);
    
case 'wtrack'
    env = 'wtrack';
    epochfilter = sprintf('(isequal($environment,''%s''))', ...
        env);
    
case 'wtrackdays'
    env = 'wtrack';
    daytype = 'wtrack';
    epochfilter = sprintf(...
        '(isequal($daytype,''%s'')) && (isequal($environment,''%s''))', ...
        daytype, env);
        
case 'exemplar_wepochs'
    env = 'wtrack';
    epochfilter = sprintf(...
        '(isequal($exemplar,1)) && (isequal($environment,''wtrack''))');
    
case 'sleepwtrackdays'
    env = 'sleep';
    daytype = 'wtrack';
    epochfilter = sprintf('(isequal($daytype,''%s'')) && (isequal($environment,''%s''))', ...
        daytype, env);
    
case 'openfield'
    env = 'openfield';
    epochfilter = sprintf('(isequal($environment,''%s''))', ...
        env);
    %         case 'all_epoch_types'
    %             epochType = {'sleep','run', 'run', 'run','run', 'run'};
    %             epochEnvironment = {'sleep', 'wtrack', 'wtrackrotated', 'openfield', ...
    %                 'sixarmtrack_right', 'sixarmtrack_left'};
    %             eptypeEnv = [epochType; epochEnvironment];
    %             epochfilter = sprintf(...
    %                 ' || ((isequal($type,''%s'')) && (isequal($environment,''%s'')))', ...
    %                 eptypeEnv{:});
    %             epochfilter = epochfilter(4:end); %trim first ||
case 'ca1ntrodes'
    tetfilter = '(isequal($valid,''yes'') && isequal($area,''ca1''))';

case 'same tet for eeg'
    eegfilter = {'geteegtet', 'theta','sametet', 1};
case 'ca1SU'
    cellFilter = '(isequal($area, ''ca1'')) && ($numspikes > 100) && (all(cellfun(''isempty'',(arrayfun(@(x) strfind(x,''mua''), $tags, ''un'', 0)))))';
case 'nonMU_cells'
    cellFilter = '($numspikes > 100) && (all(cellfun(''isempty'',(arrayfun(@(x) strfind(x,''mua''), $tags, ''un'', 0)))))';
case '>100spikes_cells'
    cellFilter = '($numspikes > 100)';
case '2s_win'
    window = [2 2];         % size of psth window (in sec)
    %make sure this matches the binsize in plot params
    binsize = .001;             % size of bins (in sec)
    frbinsize = 0.01;
    time = (-window(1)-0.5*binsize):binsize:(window(2)+0.5*binsize);
case '<1cm/s'
    timefilter{end+1} = {'get2dstate','($velocity<1)'};
case '<4cm/s'
    timefilter{end+1} = {'get2dstate','($velocity<4)'};
case '>4cm/s'
    timefilter{end+1} = {'get2dstate','($velocity>4)'};

case 'correcttrials'
    timefilter{end+1} = {'getWtracktrialstate','($correct==1)'};

case 'errortrials'
    timefilter{end+1} = {'getWtracktrialstate','($correct==0)'};

case 'outbound'
    timefilter{end+1} = {'getWtracktrialstate','($outbound==1)'};

case 'inbound'
    timefilter{end+1} = {'getWtracktrialstate','($inbound==1)'};

case 'proximalWell'
    timefilter{end+1} = {'getnearwells','($nearwell == 1)'};

case 'distalWell'
    timefilter{end+1} = {'getnearwells','($nearwell == 0)'};            

case 'excludeNoise'
    timefilter{end+1} = {'excludenoiseevents', '($noise == 0)', 'ca1noisekons', ...
        'exclpad', .5, 'stdthresh', 1, 'excludeman', 1}; %15

case 'firstToLastWellVisit'
    timefilter{end+1} = {'getpriortofirstwell', '($prefirst == 0)'};
    timefilter{end+1} = {'getpostlastwell', '($postlast == 0)'};
    
case 'excludePriorFirstWell'
    timefilter{end+1} = {'getpriortofirstwell', '($prefirst == 0)'};

case 'excludeAfterLastWell'
    timefilter{end+1} = {'getpostlastwell', '($postlast == 0)'};
case 'referenced'
    uselfptype = 'eeg';

case 'unreferenced'
    uselfptype = 'eeggnd';

case 'lickbouts'
    maxIntraBurstILI = 0.25;  % max burst ili threshold in seconds
    minBoutLicks = 3; %filter out bouts with less than boutNum licks
    % timefilt1: lick bouts, 
    timefilter{end+1} = {'getLickBout', '($lickBout == 1)', ...
        'maxIntraBurstILI', maxIntraBurstILI, ...
        'minBoutLicks', minBoutLicks};

case 'nolickbouts'
    % timefilt: immoble not lick bout
    maxIntraBurstILI = 0.25;  % max burst ili threshold in seconds
    minBoutLicks = 3; % filter out bouts with less than boutNum licks
    minTimeFromLick = .5; % time duration from closest lick
    timefilter{end+1} = {'getLickBout', ...
        sprintf('($lickBout == 0) & ($timeFromLick > %d)',minTimeFromLick), ...
        'maxIntraBurstILI', maxIntraBurstILI, 'minBoutLicks', minBoutLicks};

case 'licks'
    eventType = 'lick';            
    maxIntraBurstILI = 0.01;  % max burst ili threshold in seconds
    minBoutLicks = 1; %filter out bouts with less than boutNum licks
    % timefilt1: lick bouts, 
    timefilter{end+1} = {'getLickBout', '($lickBout == 1)', ...
        'maxIntraBurstILI', maxIntraBurstILI, ...
        'minBoutLicks', minBoutLicks};

    %% filter function specific params

case 'swrTrigSUmod'
    eventType = 'ca1rippleskons';


case 'lickSpikeXC'
    eventType = 'lick';

case 'lickTrigSUmod'
    eventType = 'lick';

case 'swrlickmod'
    filtfunction = 'swrlickmod';

case 'wavelets4-300Hz'
    waveSet = '4-300Hz';
%             wp = getWaveParams(waveSet);
case '4-300Hz_focusSWR'
    waveSet = '4-300Hz_focusSWR';
%             wp = getWaveParams(waveSet);
case '4-350Hz'
    waveSet = '4-350Hz';

case 'reactivationPLTH'
    iterator = 'singleepochanal';
    datatypes = {'spikes', 'pos', 'licks'};

case 'dfa_lickBoutSpikeCorr'

    cellpairfilter = {'allcomb', ...
        '($numspikes > 100) && (all(cellfun(''isempty'',(arrayfun(@(x) strfind(x,''mua''), $tags, ''un'', 0)))))', ...
        '($numspikes > 100) && (all(cellfun(''isempty'',(arrayfun(@(x) strfind(x,''mua''), $tags, ''un'', 0)))))';};
    filtfunction = 'dfa_lickBoutSpikeCorr';
    eventType = 'lick';
    iterator = 'singlecellanal';
    datatypes = {'spikes'};
    options = {'savefigas', 'png', 'eventName',eventType};

case 'dfa_lickphaseSUclustering'
    eventType = 'lick';
    filtfunction = 'dfa_lickphaseSUclustering';
    iterator = 'singlecellanal';
    datatypes = {'lick', 'spikes', 'cellinfo', 'DIO', 'task'};
    options = {'savefigas', 'png', 'eventName',eventType};

case 'lickSpikesXC'
    eventType = 'lick';
    tmax = 1;
    bin = .01;
    numShufs = 200;

case 'dfa_lickXCorrSpikes'
    filtfunction = 'dfa_lickXCorrSpikes';
    iterator = 'singleDayCellAnal';
    datatypes = {'lick', 'spikes', 'cellinfo'};
    options = {'savefigas', 'png', 'bin', bin, 'tmax',tmax,'eventType',eventType, ...
        'numShufs', numShufs};

case 'dfa_licktrigspiking'
    window = [1 1];
    binsize = .001;
    eventType = 'lick';
    options = {'savefigas', 'png', 'binsize', binsize, 'window', window, 'eventName', ...
        eventType};
    filtfunction = 'dfa_licktrigspiking';
    iterator = 'singlecellanal';
    datatypes = {'lick', 'task', 'DIO', 'spikes', 'cellinfo'};

case 'lickTrigSpikingMod'
    respwin = [0 200]; % response period in ms rel to swr on
    basewin = [-300 -100]; % baseline period in ms rel to swr on
    minNumSwr = 10;
    % minNumSwrSpikes = 10;
    nshuffs = 1000;
    shuffms = 700;
    dmatIdx = {'lickbout'};


case 'savefigs'
    savefigs = 1;
    pausefigs = 0;
case 'pausefigs'
    savefigs = 0;
    pausefigs = 1;

case 'dfa_plotDataChunks'
    splitSize = 30; % seconds
    Yoffset = 600; 
    eventType = 'ca1rippleskons';
    consensus_numtets = 2;   % minimum # of tets for consensus event detection
    minstdthresh = 2;        % STD. how big your ripples are
    exclusion_dur = 0;  % seconds within which consecutive events are eliminated / ignored
    minvelocity = 0;
    maxvelocity = 4;
    skipNoRipChunks = 0;
    skipExist = 0;
    centerEvents = 0;
    centerOn = 'ca1ripplekons';
    savefigas = {'mfig', 'png'};

    timefilter{end+1} = {'getconstimes', '($cons == 1)', ...
        'ca1rippleskons', 1,'consensus_numtets',consensus_numtets, ...
        'minstdthresh', minstdthresh,'exclusion_dur',exclusion_dur, ...
        'minvelocity', minvelocity,'maxvelocity',maxvelocity};            

    options = {'splitSize', splitSize, 'Yoffset', Yoffset, 'minstdthresh',...
        minstdthresh, 'maxvelocity', maxvelocity, 'skipNoRipChunks', ...
        skipNoRipChunks, 'skipExist', skipExist, 'centerEvents', centerEvents, ...
        'centerOn', centerOn, 'savefigas', {'png', 'mfig'}};

    filtfunction = 'dfa_plotDataChunks';
    iterator = 'multitetrodeanal';
    datatypes = {'ca1rippleskons', 'pos', 'linpos', 'eeg', 'spikes', ...
        'tetinfo', 'cellinfo', 'task', 'ca1noisekons', 'ripple', 'theta', 'DIO'};

case 'dfa_getPeriEventVelocity'
    win = [2 2];
    eventType = 'rippleskons';
    eventSourceArea = 'ca1';
    eventDataLabel = [eventSourceArea eventType];
    consensus_numtets = 2;   % minimum # of tets for consensus event detection
    minstdthresh = 3;        % STD. how big your ripples are
    exclusion_dur = 0;  % seconds within which consecutive events are eliminated / ignored
    minvelocity = 0;
    maxvelocity = 4;

    timefilter{end+1} = {'getconstimes', '($cons == 1)', ...
        'ca1rippleskons', 1,'consensus_numtets',consensus_numtets, ...
        'minstdthresh', minstdthresh,'exclusion_dur',exclusion_dur, ...
        'minvelocity', minvelocity,'maxvelocity',maxvelocity};            
    options = {'eventtype', [eventSourceArea eventType], 'win', win};
    filtfunction = 'dfa_getPeriEventVelocity';
    iterator = 'singleepochanal';
    datatypes = {[eventSourceArea eventType], 'pos'};



case 'dfa_riptriglfp'
%             LFPtypes = {'eeggnd', 'eeg', 'theta', };
    eventSourceArea = 'ca1';
    eventType = 'rippleskons';
    LFPtypes = {'eeggnd', 'eeg'};%, 'theta', 'ripple'}; %'lowgamma', 'fastgamma', 'ripple'};
    LFPrangesHz = {'1-400', '1-400'};%, '6-9', '140-250'}; %'20-50', '65-140',};
    srate = 1500;
    win = [1.5 1.5];
    time = -win(1):1/srate:win(2);
    TF = 1;
    eventType = 'rippleskons';
    eventSourceArea = 'ca1';
    eventDataLabel = [eventSourceArea eventType];
    consensus_numtets = 2;   % minimum # of tets for consensus event detection
    minstdthresh = 3;        % STD. how big your ripples are
    exclusion_dur = 0;  % seconds within which consecutive events are eliminated / ignored
    minvelocity = 0;
    maxvelocity = 4;
    timefilter{end+1} = {'getconstimes', '($cons == 1)', ...
        'ca1rippleskons', 1,'consensus_numtets',consensus_numtets, ...
        'minstdthresh', minstdthresh,'exclusion_dur',exclusion_dur, ...
        'minvelocity', minvelocity,'maxvelocity',maxvelocity};

    options = {'eventtype', [eventSourceArea eventType], 'LFPtypes', LFPtypes, ...
        'win', win};

    filtfunction = 'dfa_riptriglfp';
    iterator = 'multitetrodeanal';
    datatypes = {[eventSourceArea eventType], LFPtypes{:}};
    %             params{end+1} = {'allLFPtypes', '1s_win'};
    %             F = createfilter('animal', Fp.animals, 'epochs', Fp.epochfilter, ...
    % 'eegtetrodes', Fp.tetfilter, 'excludetime', Fp.timefilter, 'iterator', Fp.iterator);

case 'thetaphasemod'
    iterator = 'singlecelleeganal';
    eegfilter = {'geteegtet', 'theta','sametet', 1};
    datatypes = {'spikes','theta'};

case 'dfa_occNormFiring'
    cellFilter = ...
        '($numspikes > 100) && (all(cellfun(''isempty'',(arrayfun(@(x) strfind(x,''mua''), $tags, ''un'', 0)))))';
    timefilter{end+1} = {'get2dstate', '(abs($velocity) >= 4)'};
    iterator = 'singlecellanal';
    filtfunction = 'dfa_occNormFiring';
    datatypes = {'spikes', 'linpos', 'pos', 'task'};

case 'dfa_riptrigspiking'
    eventType = 'rippleskons';
    srate = 1500;
    win = [1 1];
    binsize = .001;
    eventType = 'swr';
    time = -win(1):1/srate:win(2);
    TF = 1;
    eventSourceArea = 'ca1';
    eventDataLabel = [eventSourceArea eventType];
    consensus_numtets = 1;   % minimum # of tets for consensus event detection
    minstdthresh = 2;        % STD. how big your ripples are
    exclusion_dur = 0;  % seconds within which consecutive events are eliminated / ignored
    minvelocity = 0;
    maxvelocity = 4;
    timefilter{end+1} = {'getconstimes', '($cons == 1)', ...
        'ca1rippleskons', 1,'consensus_numtets',consensus_numtets, ...
        'minstdthresh', minstdthresh,'exclusion_dur',exclusion_dur, ...
        'minvelocity', minvelocity,'maxvelocity',maxvelocity};

    %make sure this matches the binsize in plot params
    binsize = .001;             % size of bins (in sec)
    frbinsize = 0.01;
    time = (-win(1)-0.5*binsize):binsize:(win(2)+0.5*binsize);
    welldist = [];

    options = {'TF', TF, 'window', win, 'binsize', binsize, ...
        'frbinsize', frbinsize, 'minthresh', minstdthresh, 'maxvelocity', ...
        maxvelocity, 'minvelocity', minvelocity, 'consensus_numtets', ...
        consensus_numtets, 'welldist', welldist};

    iterator = 'singlecellanal';
    filtfunction = 'dfa_riptrigspiking';
    datatypes = {'spikes', 'ca1rippleskons', 'pos', 'task'};

case 'riptrigspiking_corr'
    iterator = 'singlecellanal';
    filtfunction = 'calcxcorrmeasures';
    datatypes = {'spikes', 'ca1rippleskons', 'pos', 'task'};

case 'dfa_suCoactiveXcorr'
    iterator = 'singlecellanal';
    filtfunction = 'dfa_suCoactiveXcorr';

    cellpairfilter = {'allcomb', ...
        '($numspikes > 100) && (isequal($area, ''ca1'')) && (all(cellfun(''isempty'',(arrayfun(@(x) strfind(x,''mua''), $tags, ''un'', 0)))))', ...
        '($numspikes > 100) && (isequal($area, ''ca1'')) && (all(cellfun(''isempty'',(arrayfun(@(x) strfind(x,''mua''), $tags, ''un'', 0)))))'};

    eventType = 'swr';
    iterator = 'singlecellanal';
    datatypes = {'spikes'};
    options = {'savefigas', 'png', 'eventName',eventType};

case 'mua_calcxcorrmeasures'
    % requires 'tetrodepairs', Fp.tetpairfilter,
    tetpairfilter = {'(isequal($area, ''mec''))', ...
        '(isequal($area, ''mec''))'};
    iterator = 'singletetrodeanal';
    filtfunction = 'mua_calcxcorrmeasures';
    datatypes = {'spikes', 'linpos'};

case 'dfa_perripspikingcorr'
    % requires 'tetrodepairs', Fp.tetpairfilter,
    tetpairfilter = {'(isequal($area, ''mec''))', ...
        '(isequal($area, ''mec''))'};
    iterator = 'singletetrodeanal';
    filtfunction = 'dfa_perripspikingcorr';
    datatypes = {'spikes', 'ca1rippleskons'};
otherwise
    error('invalid filter param: %s\n', s{1});
%             filtfunction = s{1};
end

end
fprintf('---------------------\n');
w = whos;
for a = 1:length(w)
    Fp.(w(a).name) = eval(w(a).name);
end
if add_paths
    try
        Fp.paths = make_paths(Fp.filtfunction);
    catch
        disp('could not generate paths');
    end
end
end