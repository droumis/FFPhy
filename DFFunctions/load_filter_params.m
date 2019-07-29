
function filter_params = load_filter_params(params, varargin)
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
% default_params = {'all_epoch_types'};

if isa(params, 'string')
    params = {params};
elseif isa(params, 'struct')
    filter_params = params;
    try
        params = {filter_params.filtfunction};
    catch
        params = {};
    end
end
% 
% if ~isempty(default_params)
%     params = {default_params{:}, params{:}};
% end

add_params = {};
if ~isempty(varargin)
    assign(varargin{:})
end
if ~isempty(add_params)
    % add more params on top of filtfunction defaults
    params = {params{:}, add_params{:}};
end

epochfilter = '';
tetfilter = '';
cellfilter = '';
timefilter = {};
options = {};
sysDateTime = clock;
for s = params
    fprintf('%s\n', s{:})
    switch s{1}
        %% EPOCH TETRODE CELL RIPPLE filters
        case 'sleep'
            epochEnvironment = 'sleep';
            epochfilter = sprintf('(isequal($environment,''%s''))', ...
                epochEnvironment);
        case 'wtrack'
            epochEnvironment = 'wtrack';
            epochfilter = sprintf('(isequal($environment,''%s''))', ...
                epochEnvironment);
        case 'wtrackdays'
            epochEnvironment = 'wtrack';
            daytype = 'wtrack';
            epochfilter = sprintf('(isequal($daytype,''%s'')) && (isequal($environment,''%s''))', ...
                daytype, epochEnvironment);
        case 'openfield'
            epochEnvironment = 'openfield';
            epochfilter = sprintf('(isequal($environment,''%s''))', ...
                epochEnvironment);
            %         case 'all_epoch_types'
            %             epochType = {'sleep','run', 'run', 'run','run', 'run'};
            %             epochEnvironment = {'sleep', 'wtrack', 'wtrackrotated', 'openfield', ...
            %                 'sixarmtrack_right', 'sixarmtrack_left'};
            %             eptypeEnv = [epochType; epochEnvironment];
            %             epochfilter = sprintf(...
            %                 ' || ((isequal($type,''%s'')) && (isequal($environment,''%s'')))', ...
            %                 eptypeEnv{:});
            %             epochfilter = epochfilter(4:end); %trim first ||
        case 'nonref_ntrodes'
            ntAreas = {'ca1', 'mec', 'por', 'v2l', 'sub'};
            tetfilter = sprintf(' || (isequal($area,''%s''))', ntAreas{:});
            tetfilter = tetfilter(4:end); %trim first ||
        case 'same tet for eeg'
            eegfilter = {'geteegtet', 'theta','sametet', 1};
        case 'nonMU_cells'
            cellfilter = ...
                '($numspikes > 100) && (all(cellfun(''isempty'',(arrayfun(@(x) strfind(x,''mua''), $tags, ''un'', 0)))))';
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
                'exclpad', 1, 'stdthresh', 15}; %15
            
        case 'excludePriorFirstWell'
            timefilter{end+1} = {'getpriortofirstwell', '($prefirst == 0)'};
            
        case 'excludeAfterLastWell'
            timefilter{end+1} = {'getpostlastwell', '($postlast == 0)'};            
            
        case 'ripples'
            TF = 1;
            eventtype = 'rippleskons';
            eventSourceArea = 'ca1';
            eventDataLabel = [eventSourceArea eventtype];
            consensus_numtets = 2;   % minimum # of tets for consensus event detection
            minstdthresh = 2;        % STD. how big your ripples are
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
        case 'wavelets4-300Hz'
            waveSet = '4-300Hz';

            
        case 'behavestate'
            
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
    filter_params.(w(a).name) = eval(w(a).name);
end
if add_paths
    try
        filter_params.paths = make_paths(filter_params.filtfunction);
    catch
        disp('could not generate paths');
    end
end
end