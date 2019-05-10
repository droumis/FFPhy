
function filter_params = load_filter_params(param_set)
%% param_set : cell array of strings
% loads variables for each case of the input param_set to caller workspace
% i.e.:: load_filter_params({'occnormfiring', 'riptrighist'})
% succesors will take priority when var names clash

if ~isa(param_set,'cell')
    param_set = {param_set};
end
for s = param_set
    switch s{1}
        %% EPOCH TETRODE CELL RIPPLE filters
        case 'wtrack'
            epochEnvironment = 'wtrack';
            epochfilter = sprintf('(isequal($environment,''%s''))', ...
                epochEnvironment);
        case 'all_epoch_types'
            epochType = {'sleep','run', 'run', 'run','run', 'run'};
            epochEnvironment = {'sleep', 'wtrack', 'wtrackrotated', 'openfield', ...
                'sixarmtrack_right', 'sixarmtrack_left'};
            eptypeEnv = [epochType; epochEnvironment];
            epochfilter = sprintf(...
                ' || ((isequal($type,''%s'')) && (isequal($environment,''%s'')))', ...
                eptypeEnv{:});
            epochfilter = epochfilter(4:end); %trim first ||
        case 'nonref_ntrodes'
            ntAreas = {'ca1', 'mec', 'por', 'v2l', 'sub'};
            tetfilter = sprintf(' || (isequal($area,''%s''))', ntAreas{:});
            tetfilter = tetfilter(4:end); %trim first ||
        case 'nonMU_cells'
            cellfilter = ...
                '($numspikes > 100) && (all(cellfun(''isempty'',(arrayfun(@(x) strfind(x,''mua''), $tags, ''un'', 0)))))';
        case '>100spikes_cells'
            cellfilter = '($numspikes > 100)';
        case 'ripples'
            TF = 1;
            eventtype = 'rippleskons';
            eventSourceArea = 'ca1';
            eventDataLabel = [eventSourceArea eventtype];
            consensus_numtets = 1;   % minimum # of tets for consensus event detection
            minstdthresh = 2;        % STD. how big your ripples are
            exclusion_dur = 0.5;  % seconds within which consecutive events are eliminated / ignored
            minvelocity = 0;
            maxvelocity = 4;
            window = [0.5 0.5];         % size of psth window (in sec)
            %make sure this matches the binsize in plot params
            binsize = .001;             % size of bins (in sec)
            frbinsize = 0.02;
            time = (-window(1)-0.5*binsize):binsize:(window(2)+0.5*binsize);
            welldist = [];
            timefilter{1} = {'getconstimes', '($cons == 1)', ...
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
            
        case 'ratemaps'
            timefilter{1} = {'get2dstate', '(abs($velocity) >= 4)'};
            iterator = 'singlecellanal';
            filtfunction = 'dfa_occNormFiring';
            datatypes = {'spikes', 'linpos', 'pos', 'task'};
            
        case 'riptrigspiking'
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
    end
end
w = whos;
for a = 1:length(w)
    filter_params.(w(a).name) = eval(w(a).name);
end
end