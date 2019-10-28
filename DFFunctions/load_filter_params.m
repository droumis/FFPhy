
function Fp = load_filter_params(Fp, varargin)
% ex. filter_params = load_filter_params(Fp, 'add_params', {'ripples', 'wtrack'});
% loads filter parameter variables from group keyword as a convenience

% params input can either be a struct with a valid 'filtfunction' field
% or a cell of strings, each specifying some param group
% if Fp is a struct, include additional filter groups as a cell array 
% of strings with 'add_params' varargin
% succesors will take priority when var names clash

% Author: Demetris Roumis June 2019

add_paths = 1;
set_filt_func = 1;
if ~isempty(varargin)
    assign(varargin{:})
end

epochfilter = '';
tetfilter = '';
cellfilter = '';
timefilter = {};
options = {};
sysDateTime = clock;
for s = Fp.params
    fprintf('%s\n', s{:})
    switch s{1}
        %% EPOCH TETRODE CELL RIPPLE filters
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
            epochfilter = sprintf('(isequal($daytype,''%s'')) && (isequal($environment,''%s''))', ...
                daytype, env);
        case 'exemplar_wepochs'
            env = 'wtrack';
            epochfilter = sprintf('(isequal($exemplar,1)) && (isequal($environment,''wtrack''))');
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
        case 'valid_ntrodes'
%             ntAreas = {'ca1', 'mec', 'ref'}; %, 'por', 'v2l', 'sub'};
            tetfilter = '(isequal($valid,''yes'') && (isequal($area,''ca1'') || isequal($area,''mec'')))';
%             tetfilter = tetfilter(4:end); %trim first ||
        case 'same tet for eeg'
            eegfilter = {'geteegtet', 'theta','sametet', 1};
        case 'ca1SU'
            cellfilter = '(isequal($area, ''ca1'')) && ($numspikes > 100) && (all(cellfun(''isempty'',(arrayfun(@(x) strfind(x,''mua''), $tags, ''un'', 0)))))';
        case 'nonMU_cells'
            cellfilter = '($numspikes > 100) && (all(cellfun(''isempty'',(arrayfun(@(x) strfind(x,''mua''), $tags, ''un'', 0)))))';
        case '>100spikes_cells'
            cellfilter = '($numspikes > 100)';
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
            timefilter{end+1} = {'getnearwells','($nearwell == 0)'};
        
        case 'distalWell'
            timefilter{end+1} = {'getnearwells','($nearwell == 1)'};            
            
        case 'excludeNoise'
            timefilter{end+1} = {'excludenoiseevents', '($noise == 0)', 'ca1noisekons', ...
                'exclpad', 1, 'stdthresh', 15, 'excludeman', 1}; %15
            
        case 'excludePriorFirstWell'
            timefilter{end+1} = {'getpriortofirstwell', '($prefirst == 0)'};
            
        case 'excludeAfterLastWell'
            timefilter{end+1} = {'getpostlastwell', '($postlast == 0)'};
        case 'referenced'
            uselfptype = 'eeg';
        case 'unreferenced'
            uselfptype = 'eeggnd';
        case 'ripples'
            TF = 1;
            eventtype = 'rippleskons';
            eventSourceArea = 'ca1';
            eventDataLabel = [eventSourceArea eventtype];
            consensus_numtets = 2;   % minimum # of tets for consensus event detection
            minstdthresh = 3;        % STD. how big your ripples are
            exclusion_dur = 0;  % seconds within which consecutive events are eliminated / ignored
            minvelocity = 0;
            maxvelocity = 4;
            
            welldist = [];
            timefilter{end+1} = {'getconstimes', '($cons == 1)', ...
                'ca1rippleskons', 1,'consensus_numtets',consensus_numtets, ...
                'minstdthresh', minstdthresh,'exclusion_dur',exclusion_dur, ...
                'minvelocity', minvelocity,'maxvelocity',maxvelocity};
            %     F = setfilterfunction(F, Fp.filtfunction, Fp.datatypes,...
            %         'TF',Fp.TF,'window',Fp.window, ...
            %         'binsize',Fp.binsize,'frbinsize',Fp.frbinsize,'minthresh', ...
            %         Fp.minstdthresh,'maxvelocity',Fp.maxvelocity,'minvelocity', ...
            %         Fp.minvelocity, 'consensus_numtets',Fp.consensus_numtets,'welldist', ...
            %         Fp.welldist);
            %% filter function specific params
        case 'dfa_reactivationPSTH'
            bin = 0.01; % seconds
            win = [-2 2];
            
            consensus_numtets = 2;   % minimum # of tets for consensus event detection
            minstdthresh = 3;        % STD. how big your ripples are
            exclusion_dur = 0;  % seconds within which consecutive events are eliminated / ignored
            minvelocity = 0;
            maxvelocity = 4;
            maxIntraBurstILI = 0.25;  % max burst ili threshold in seconds
            minBoutLicks = 3; % filter out bouts with less than boutNum licks
            minTimeFromLick = .5; % time duration from closest lick
            
            iterator = 'singleepochanal';
            datatypes = {'spikes', 'lick', 'ca1rippleskons', 'cellinfo'};
            options = {'cellfilter', cellfilter, 'win', win, 'bin', bin, ...
                'swrTimeFilter', ...
                [{{'getconstimes', '($cons == 1)', ...
                'ca1rippleskons', 1,'consensus_numtets',consensus_numtets, ...
                'minstdthresh', minstdthresh,'exclusion_dur',exclusion_dur, ...
                'minvelocity', minvelocity,'maxvelocity',maxvelocity}}, ...
                {{'getpriortofirstwell', '($prefirst == 0)'}}], ...
                'burstTimeFilter', ...
                {{'getLickBout', '($lickBout == 1)', 'maxIntraBurstILI', maxIntraBurstILI, ...
                'minBoutLicks', minBoutLicks}}};
            
            
        case 'swrlickmod'
            filtfunction = 'swrlickmod';
            
        case 'wavelets4-300Hz'
            waveSet = '4-300Hz';
            wp = getWaveParams(waveSet);
            
        case 'reactivationPLTH'
            iterator = 'singleepochanal';
            datatypes = {'spikes', 'pos', 'licks'};

        case 'dfa_lickBoutSpikeCorr'

            cellpairfilter = {'allcomb', ...
                '($numspikes > 100) && (all(cellfun(''isempty'',(arrayfun(@(x) strfind(x,''mua''), $tags, ''un'', 0)))))', ...
                '($numspikes > 100) && (all(cellfun(''isempty'',(arrayfun(@(x) strfind(x,''mua''), $tags, ''un'', 0)))))';};
            filtfunction = 'dfa_lickBoutSpikeCorr';
            eventName = 'lick';
            iterator = 'singlecellanal';
            datatypes = {'spikes'};
            options = {'savefigas', 'png', 'eventName',eventName};
            
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
                
        case 'dfa_lickphaseSUclustering'
            eventName = 'lick';
            filtfunction = 'dfa_lickphaseSUclustering';
            iterator = 'singlecellanal';
            datatypes = {'lick', 'spikes', 'cellinfo', 'DIO', 'task'};
            options = {'savefigas', 'png', 'eventName',eventName};
            
        case 'dfa_lickXCorrSpikes'
            eventName = 'lick';
            filtfunction = 'dfa_lickXCorrSpikes';
            iterator = 'singlecellanal';
            tmax = 1;
            bin = .02;
            datatypes = {'lick', 'spikes', 'cellinfo'};
            options = {'savefigas', 'png', 'bin', bin, 'tmax',tmax,'eventName',eventName};
            
        case 'dfa_licktrigspiking'
            window = [1 1];
            binsize = .001;
            eventName = 'lick';
            options = {'savefigas', 'png', 'binsize', binsize, 'window', window, 'eventName', ...
                eventName};
            filtfunction = 'dfa_licktrigspiking';
            iterator = 'singlecellanal';
            datatypes = {'lick', 'task', 'DIO', 'spikes', 'cellinfo'};
            
        case 'dfa_lickswrcorr'
           % /home/droumis/Src/Matlab/filterframework_dr/DFFunctions/dfa_lickswrcorr.m
            bin = .01;
            tmax = 1;
            sw1 = bin*2;
            sw2 = .250;
            sigperc = .975;
            eventtype = 'ca1rippleskons';
            consensus_numtets = 2;   % minimum # of tets for consensus event detection
            minstdthresh = 2;        % STD. how big your ripples are
            exclusion_dur = 0;  % seconds within which consecutive events are eliminated / ignored
            minvelocity = 0;
            maxvelocity = 4;
            timefilter{end+1} = {'getconstimes', '($cons == 1)', ...
                'ca1rippleskons', 1,'consensus_numtets',consensus_numtets, ...
                'minstdthresh', minstdthresh,'exclusion_dur',exclusion_dur, ...
                'minvelocity', minvelocity,'maxvelocity',maxvelocity};            
            options = {'savefigas', 'png', 'bin', bin, 'tmax', tmax, 'savefigs', savefigs, ...
                'pausefigs', pausefigs};
            filtfunction = 'dfa_lickswrcorr';
            iterator = 'singleepochanal';
            datatypes = {'ca1rippleskons','task', 'lick', 'DIO'};
        case 'savefigs'
            savefigs = 1;
            pausefigs = 0;
        case 'pausefigs'
            savefigs = 0;
            pausefigs = 1;
        case 'dfa_plotDataChunks'
            splitSize = 30; % seconds
            Yoffset = 600; 
            eventtype = 'ca1rippleskons';
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
            eventtype = 'rippleskons';
            eventSourceArea = 'ca1';
            eventDataLabel = [eventSourceArea eventtype];
            consensus_numtets = 2;   % minimum # of tets for consensus event detection
            minstdthresh = 3;        % STD. how big your ripples are
            exclusion_dur = 0;  % seconds within which consecutive events are eliminated / ignored
            minvelocity = 0;
            maxvelocity = 4;
            
            timefilter{end+1} = {'getconstimes', '($cons == 1)', ...
                'ca1rippleskons', 1,'consensus_numtets',consensus_numtets, ...
                'minstdthresh', minstdthresh,'exclusion_dur',exclusion_dur, ...
                'minvelocity', minvelocity,'maxvelocity',maxvelocity};            
            options = {'eventtype', [eventSourceArea eventtype], 'win', win};
            filtfunction = 'dfa_getPeriEventVelocity';
            iterator = 'singleepochanal';
            datatypes = {[eventSourceArea eventtype], 'pos'};
            
        case 'dfa_riptriglfp'
%             LFPtypes = {'eeggnd', 'eeg', 'theta', };
            eventSourceArea = 'ca1';
            eventtype = 'rippleskons';
            LFPtypes = {'eeggnd', 'eeg'};%, 'theta', 'ripple'}; %'lowgamma', 'fastgamma', 'ripple'};
            LFPrangesHz = {'1-400', '1-400'};%, '6-9', '140-250'}; %'20-50', '65-140',};
            srate = 1500;
            win = [2 2];
            time = -win(1):1/srate:win(2);
            TF = 1;
            eventtype = 'rippleskons';
            eventSourceArea = 'ca1';
            eventDataLabel = [eventSourceArea eventtype];
            consensus_numtets = 2;   % minimum # of tets for consensus event detection
            minstdthresh = 3;        % STD. how big your ripples are
            exclusion_dur = 0;  % seconds within which consecutive events are eliminated / ignored
            minvelocity = 0;
            maxvelocity = 4;
            timefilter{end+1} = {'getconstimes', '($cons == 1)', ...
                'ca1rippleskons', 1,'consensus_numtets',consensus_numtets, ...
                'minstdthresh', minstdthresh,'exclusion_dur',exclusion_dur, ...
                'minvelocity', minvelocity,'maxvelocity',maxvelocity};
            
            options = {'eventtype', [eventSourceArea eventtype], 'LFPtypes', LFPtypes, ...
                'win', win};
            
            filtfunction = 'dfa_riptriglfp';
            iterator = 'multitetrodeanal';
            datatypes = {[eventSourceArea eventtype], LFPtypes{:}};
            %             params{end+1} = {'allLFPtypes', '1s_win'};
            %             F = createfilter('animal', Fp.animals, 'epochs', Fp.epochfilter, ...
            % 'eegtetrodes', Fp.tetfilter, 'excludetime', Fp.timefilter, 'iterator', Fp.iterator);
            
        case 'thetaphasemod'
            iterator = 'singlecelleeganal';
            eegfilter = {'geteegtet', 'theta','sametet', 1};
            datatypes = {'spikes','theta'};
            
        case 'dfa_occNormFiring'
            cellfilter = ...
                '($numspikes > 100) && (all(cellfun(''isempty'',(arrayfun(@(x) strfind(x,''mua''), $tags, ''un'', 0)))))';
            timefilter{end+1} = {'get2dstate', '(abs($velocity) >= 4)'};
            iterator = 'singlecellanal';
            filtfunction = 'dfa_occNormFiring';
            datatypes = {'spikes', 'linpos', 'pos', 'task'};
            
        case 'dfa_riptrigspiking'
            eventtype = 'rippleskons';
            srate = 1500;
            win = [1 1];
            binsize = .001;
            eventName = 'swr';
            time = -win(1):1/srate:win(2);
            TF = 1;
            eventtype = 'rippleskons';
            eventSourceArea = 'ca1';
            eventDataLabel = [eventSourceArea eventtype];
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
            
            eventName = 'swr';
            iterator = 'singlecellanal';
            datatypes = {'spikes'};
            options = {'savefigas', 'png', 'eventName',eventName};
            
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
            filtfunction = s{1};
    end
end
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