
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
        case 'occnormfiring_openfield'
            filtfunction = 'dfa_occNormFiring';
            datatypes = {'spikes', 'linpos', 'pos', 'task'};
            epochType = {'run'};
            epochEnvironment = {'openfield'};
            ntAreas = {'ca1', 'mec', 'por', 'v2l', 'sub'};
            runOccnormfiring = 1;
            cellfilter = ...
            '($numspikes > 100) && (all(cellfun(''isempty'',(arrayfun(@(x) strfind(x,''mua''), $tags, ''un'', 0)))))';
            eptypeEnv = [epochType; epochEnvironment];
            epochfilter = sprintf(...
                ' || ((isequal($type,''%s'')) && (isequal($environment,''%s'')))', ...
                eptypeEnv{:});
            epochfilter = epochfilter(4:end); %trim first ||
            tetfilter = sprintf(' || (isequal($area,''%s''))', ntAreas{:}); 
            tetfilter = tetfilter(4:end); %trim first ||
            timefilter{1} = {'get2dstate', '(abs($velocity) >= 4)'};
            iterator = 'singlecellanal';
            
        case 'riptrigspiking'
            epochType = {'run'};
            epochEnvironment = {'wtrack', 'wtrackrotated','openfield'};
            ntAreas = {'ca1', 'mec', 'por', 'v2l', 'sub'};
            eptypeEnv = [epochType; epochEnvironment];
            TF = 1; 
            eventtype = 'rippleskons';
            eventSourceArea = 'ca1';
            eventDataLabel = [eventSourceArea eventtype];

            window = [0.5 0.5];         % size of psth window (in sec)
             %make sure this matches the binsize in plot params
            binsize = .001;             % size of bins (in sec)
            frbinsize = 0.02;
            time = (-window(1)-0.5*binsize):binsize:(window(2)+0.5*binsize);
            welldist = [];
            
            cellfilter = ...
'($numspikes>100) && (all(cellfun(''isempty'',(arrayfun(@(x) strfind(x,''mua''),$tags,''un'',0)))))';
            epochfilter = sprintf(...
' || ((isequal($type,''%s'')) && (isequal($environment,''%s'')))',eptypeEnv{:});
            epochfilter = epochfilter(4:end); %trim first ||
            tetfilter = sprintf(' || (isequal($area,''%s''))', ntAreas{:}); 
            tetfilter = tetfilter(4:end); %trim first ||
            iterator = 'singlecellanal';
            filtfunction = 'dfa_riptrigspiking';
            datatypes = {'spikes', 'ca1rippleskons', 'pos', 'task'};
            
            consensus_numtets = 1;   % minimum # of tets for consensus event detection
            minstdthresh = 2;        % STD. how big your ripples are
            exclusion_dur = 0.5;  % seconds within which consecutive events are eliminated / ignored
            minvelocity = 0;
            maxvelocity = 4;
            timefilter{1} = {'getconstimes', '($cons == 1)', ...
                'ca1rippleskons', 1,'consensus_numtets',consensus_numtets, ...
                'minstdthresh', minstdthresh,'exclusion_dur',exclusion_dur, ...
                'minvelocity', minvelocity,'maxvelocity',maxvelocity};

        case 'riptrigspiking_withMU'
            epochType = {'run'};
            epochEnvironment = {'wtrack'};
            ntAreas = {'ca1', 'mec', 'por', 'v2l', 'sub'};
            eptypeEnv = [epochType; epochEnvironment];
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
            
            cellfilter = '($numspikes>100)';
            epochfilter = sprintf(...
' || ((isequal($type,''%s'')) && (isequal($environment,''%s'')))',eptypeEnv{:});
            epochfilter = epochfilter(4:end); %trim first ||
            tetfilter = sprintf(' || (isequal($area,''%s''))', ntAreas{:}); 
            tetfilter = tetfilter(4:end); %trim first ||
            iterator = 'singlecellanal';
            filtfunction = 'dfa_riptrigspiking';
            datatypes = {'spikes', 'ca1rippleskons', 'pos', 'task'};
            timefilter{1} = {'getconstimes', '($cons == 1)', ...
                'ca1rippleskons', 1,'consensus_numtets',consensus_numtets, ...
                'minstdthresh', minstdthresh,'exclusion_dur',exclusion_dur, ...
                'minvelocity', minvelocity,'maxvelocity',maxvelocity};
            
        case 'riptrigspiking_corr'
            epochType = {'sleep','run', 'run'};
            epochEnvironment = {'sleep', 'wtrack', 'openfield'};
            ntAreas = {'ca1', 'mec', 'por', 'v2l', 'sub'};
            eptypeEnv = [epochType; epochEnvironment];
            TF = 1; 
            eventtype = 'rippleskons';
            eventSourceArea = 'ca1';
            eventDataLabel = [eventSourceArea eventtype];
            consensus_numtets = 2;   % minimum # of tets for consensus event detection
            minstdthresh = 2;        % STD. how big your ripples are
            exclusion_dur = 0.5;  % seconds within which consecutive events are eliminated / ignored
            minvelocity = 0;
            maxvelocity = 4;
            
            cellpairfilter = {'allcomb', '(isequal($area, ''mec''))', ...
            '(isequal($area, ''mec''))'};
            epochfilter = sprintf(...
' || ((isequal($type,''%s'')) && (isequal($environment,''%s'')))',eptypeEnv{:});
            epochfilter = epochfilter(4:end); %trim first ||
%             tetfilter = sprintf(' || (isequal($area,''%s''))', ntAreas{:}); 
%             tetfilter = tetfilter(4:end); %trim first ||
            iterator = 'singlecellanal';
            filtfunction = 'calcxcorrmeasures';
            datatypes = {'spikes', 'ca1rippleskons', 'pos', 'task'};
            timefilter{1} = {'getconstimes', '($cons == 1)', ...
                'ca1rippleskons', 1,'consensus_numtets',consensus_numtets, ...
                'minstdthresh', minstdthresh,'exclusion_dur',exclusion_dur, ...
                'minvelocity', minvelocity,'maxvelocity',maxvelocity};
            
            
        case 'mua_calcxcorrmeasures'
            consensus_numtets = 2;   % minimum # of tets for consensus event detection
            minstdthresh = 2;        % STD. how big your ripples are
            exclusion_dur = 0.5;  % seconds within which consecutive events are eliminated / ignored
            minvelocity = 0;
            maxvelocity = 4;
            epochType = {'run'};
            epochEnvironment = {'wtrack'};
            epochfilter = '(isequal($environment, ''wtrack''))';
            tetpairfilter = {'(isequal($area, ''mec''))', '(isequal($area, ''mec''))'};
            iterator = 'singletetrodeanal';
            timefilter{1} = {'getconstimes', '($cons == 1)', ...
                'ca1rippleskons', 1,'consensus_numtets',consensus_numtets, ...
                'minstdthresh', minstdthresh,'exclusion_dur',exclusion_dur, ...
                'minvelocity', minvelocity,'maxvelocity',maxvelocity};
            filtfunction = 'mua_calcxcorrmeasures';
            datatypes = {'spikes', 'linpos'};
    end
end
w = whos;
for a = 1:length(w)
filter_params.(w(a).name) = eval(w(a).name);
end
end