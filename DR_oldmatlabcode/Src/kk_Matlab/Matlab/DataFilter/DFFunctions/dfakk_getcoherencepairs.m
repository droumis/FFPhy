function [out] = dfakk_getcoherencepairs(index, excludeperiods, eeg, tetinfo, freqrange, tetfilter1, tetfilter2, varargin)


% default
windowsize = 0.5;
params.tapers = [3 5];
params.Fs = 1500;
params.fpass = [0 400];
params.err = 0;
params.trialave = 1;

% assign the options


for option = 1:2:length(varargin)-1
    switch varargin{option}          
        case 'windowsize'             % in seconds
            windowsize = varargin{option+1};
        case 'params'
            params = varargin{option+1};
        otherwise
            error(['Option ''', varargin{option}, ''' not defined']);
    end
end


if isempty(tetfilter1) || isempty(tetfilter2)
    error('need to specify gammatetfilter!')
end


%% First receive valid periods for the epoch's day. %%%

day = index(1,1);
epoch = index(1,2);

    tetlist1 =  evaluatefilter(tetinfo{day}{epoch},tetfilter1); 
    tetlist2 =  evaluatefilter(tetinfo{day}{epoch},tetfilter2);
    alltet = unique([tetlist1 ; tetlist2])';

if isempty(tetlist1) || isempty(tetlist2) || (length(unique(alltet)) < 2)
    disp(sprintf('day %d epoch %d : no valid pairs',round(day),round(epoch)));
    out.pairs = [];
    out.freqs = [];
    out.coherograms = [];
    return
end

g = eeg{day}{epoch}{tetlist1(1)};
times = geteegtimes(g);         

%obtain includetimes
includetimes = ~isExcluded(times, excludeperiods);
includeperiods = vec2list(includetimes,times);

% (if there are no valid times this epoch, report no outputs) 
if sum(includetimes)==0
    disp(sprintf('day %d epoch %d: no includetimes',round(day),round(epoch)));
    out.pairs = [];
    out.freqs = [];
    out.coherograms = [];

    return
end

%%%%%%%%


%% Second, iterate through each valid period, chop into windows, and collect into eegwindows.

%initialize output vector
eegwindows = {};
for tet = alltet
    eegwindows{tet} = [];
end
    
for p=1:size(includeperiods,1)
    
    duration = includeperiods(p,2)-includeperiods(p,1);
    numwindows = floor(duration / windowsize);
    
    if numwindows
        
        starttime = includeperiods(p,1);
        for tet = alltet
            
            times = geteegtimes(eeg{day}{epoch}{tet});
            data = eeg{day}{epoch}{tet}.data;
            winsize_samp = windowsize * round(eeg{day}{epoch}{tet}.samprate);
            startind = lookup(starttime,times);
            
            for w=1:numwindows
                eegwindows{tet} = [eegwindows{tet}  data( (startind + (w-1)*winsize_samp):(startind + w*winsize_samp - 1) )];
            end
        end
    end
end

%% Thirdly, set up all unique 1-2 pairs and the calculate coherences.

pairs = combvec(tetlist1',tetlist2')';
% remove self-pairs and repeat pairs
dummy = size(pairs,1);
for p = 1:size(pairs,1)
    ind = dummy - p + 1;
    if pairs(ind,1)==pairs(ind,2) || logical(rowfind(fliplr(pairs(ind,:)),pairs(1:(ind-1),:)))
        pairs(ind,:) = [];
    end
end
pairs = [repmat([day epoch],size(pairs,1),1) pairs]; 
% calculate coherences
for p = 1:size(pairs,1)
   [C(p,:),~,~,~,~,freqs]=coherencyc(eegwindows{pairs(p,3)},eegwindows{pairs(p,4)},params);
   disp('coherence calculated!')
end


%% Lastly send to output.

out.pairs = pairs;
out.freqs = freqs
out.coherograms = C;


end