function [out] = dfakk_getgammagndtrigspiking2(index, excludeperiods, lowgammagnd, gammal, spikes, tetinfo,varargin)
% [out] = dfakk_getgammatrigspiking(index, excludeperiods, lowgamma, gammal, spikes, tetinfo, varargin)

% adapted from dfakk_getriptrigspiking

% 

% Note that gamma times are inherited from kk_getgammatimes -- thus
% kk_getgammatimes is where you should specify minthresh (i.e. 2 SD)

% This function differs fundamentally from ripple PSTH (i.e. dfakk_getriptrigspiking)

%
%   index [day epoch tetrode cell tetrode cell]
%
%   options are
%	'minthresh',
%		     specifies the minimum threshold of a valid ripple event

%   'gammatetfilter' -- these tetrodes will be used to extract gamma phase
%                    -- note that Carr--Frank-2012 uses CA3 w/ >= 2
%                    principal cellsc

%   out = out.out   An R x C sized matrix where each entry is the number of
%                   spikes cell C fired during ripple R
%         out.index [D E T C], gives the identity of the cells for each
%                   column in out.out
%         out.times [starttime endtime], givest the starttime and endtime
%                   of the gammal for each row in out.out

% assign the options
ripple_flag = '';
minthresh = 0;
binsize = 0.005;
freqrange = [30 50];

for option = 1:2:length(varargin)-1
    switch varargin{option}
        case 'gammatetfilter'
            gammatetfilter = varargin{(option+1)};   % tetrodes used to extract the phase
        case 'local'
            local_flag = varargin{(option+1)};   % see dfskk_gammatrigspiking for explanation of this flag
        case 'freqrange'
            freqrange = varargin{(option+1)};   % range for coherence calculation for pairs regions (for 'local' option)
                                                % in format [lowfreq , highfreq]
        case 'ripple'
            ripple_flag = varargin{(option+1)};
        otherwise
            error(['Option ''', varargin{option}, ''' not defined']);
    end
end

%% First find de novo here valid lowgamma periods for each gammatetfilter tetrode.
    % Then take intersection of valid lowgamma periods with includetimes.

day = index(1,1);
epoch = index(1,2);

if ~isempty(gammatetfilter)
    tetlist =  evaluatefilter(tetinfo{day}{epoch},gammatetfilter)'; 
else
    error('need to specify gammatetfilter!')
end

if isempty(tetlist)
    disp(sprintf('no valid gammatetfilter tetrodes this epoch'));
    out.indices_coherence = [];
    out.fields_indices_coherence = '[<cellindex> , <reftet num> , <coherence b/t cell tetrode and reftet> ]';
    out.spikephases = [];
    out.numspikes_total = [];
    out.numspikes_gammal = [];
    return
end


g = gammal{day}{epoch}{tetlist(1)};
times = g.timerange(1):0.001:g.timerange(end);         % 1 ms step vector

%obtain includetimes
includetimes = ~isExcluded(times, excludeperiods);     % 1s and 0s: ones = included, zeros = excluded

% (if there are no valid periods this epoch, report no outputs) 
if sum(includetimes)==0
    disp(sprintf('no includetimes, day %d epoch %d',day,epoch));
    disp(sprintf('%d tetrodes this epoch',length(tetlist)));
    out.indices_coherence = [];
    out.fields_indices_coherence = '[<cellindex> , <reftet num> , <coherence b/t cell tetrode and reftet> ]';
    out.spikephases = [];
    out.numspikes_total = [];
    out.numspikes_gammal = [];
    return
end

% if doing ripple-only (i.e. if doing theta), get all 2 SD periods for each gammatetfilter tetrode, take intersection

gammaperiods = {};


for tet = tetlist
    if ~strcmp(ripple_flag,'rip')
        eventperiods = [gammal{day}{epoch}{tet}.starttime gammal{day}{epoch}{tet}.endtime];
        eventvec = list2vec(eventperiods,times);
        % take intersection with permissible times
        gammaperiods{tet} = vec2list(eventvec & includetimes,times);
        % report how much of epoch is 2 SD gamma
        percentepoch = round(100*sum(gammaperiods{tet}(:,2)-gammaperiods{tet}(:,1))/(times(end)-times(1)));
        disp(sprintf('d %d e %d t %d, %d%% gamma',day,epoch,tet,percentepoch));
    % if ripple-only, then take gamma regardless of power
    else
        gammaperiods{tet} = vec2list(includetimes,times);
    end
end

%%%%%%%%


%% Then assign spikes to gamma phases, then output.

indices_coherence  = []  ;   % [<cellindex> , <reftet num> , <coherence b/t cell tetrode and reftet> ]
                             % [ day epoch tet cell reftet coherence ] 
% iterate over units
for i = 1:size(index,1)
   
        % if CA1 (same region as gammafiltertet), choose the cell's own
        % tetrode ("local") to get gamma periods
        if any(index(i,3) == tetlist)
            reftet = index(i,3);
            indices_coherence(i,:) = [ index(i,:) , index(i,3) ,  1];  % for reporting purposes
        % if some other region (say CA3), identify CA1 tetrode that is most coherent to get events    
        else
            % check if you've previously calculated-- if so, just retrieve
            if ~isempty(indices_coherence) && ~isempty(find(indices_coherence(:,3)==index(i,3)))
                ind = find(indices_coherence(:,3)==index(i,3));
                    ind = ind(1);
                reftet = indices_coherence(ind,5);
                % for reporting later
                indices_coherence(i,:) = [index(i,:)  reftet indices_coherence(ind,6)];
            % if not previously calculated, then calculate all coherences
            % to the cell's tetrode
            else
                gamma1 = lowgammagnd{index(i,1)}{index(i,2)}{index(i,3)};
                coherences = [];
                for kk = 1:length(tetlist)
                    gamma2 = lowgammagnd{index(i,1)}{index(i,2)}{tetlist(kk)};
                    coherences = [coherences ; tetlist(kk) meancoherence(gamma1,gamma2,freqrange)];
                end
                coherences = sortrows(coherences,-2);
                reftet = coherences(1,1);
                indices_coherence(i,:) = [index(i,:) reftet coherences(1,2)];
            end
     
            disp(sprintf('tet %d: assigned tet %d for ref (coherence %d)',index(i,3), reftet, coherences(1,2)));

        end   

        filtertimes = gammaperiods{reftet};
        
    % get spikes and their phases
    if ~isempty(spikes{index(i,1)}{index(i,2)}{index(i,3)}{index(i,4)})
        spiketimes = spikes{index(i,1)}{index(i,2)}{index(i,3)}{index(i,4)}.data(:,1);
            includedspike_times = spiketimes(logical(isExcluded(spiketimes,filtertimes)));
        filteegtimes=geteegtimes(lowgammagnd{index(i,1)}{index(i,2)}{reftet});
        includedspike_indices = lookup(includedspike_times,filteegtimes);
        reftetphases = lowgammagnd{index(i,1)}{index(i,2)}{reftet}.data(:,2);
    
        spikephases{i} = reftetphases(includedspike_indices);                            % list of phase times of spikes
        numspikes_gammal(i,:) = length(includedspike_indices);             % # spikes firing during gammal period
        numspikes_total(i,:) = length(spiketimes);                                          % # all spikes fired in the epoch 
    else
        disp(sprintf('d%de%dt%dc%d, cell in iterator but empty spike data cell this epoch..',index(i,1),index(i,2),index(i,3),index(i,4)))
        continue
    end
    
end


%% Lastly send to output.

if exist('spikephases')
    out.indices_coherence = indices_coherence;
    out.fields_indices_coherence = '[<cellindex> , <reftet num> , <coherence b/t cell tetrode and reftet> ]';
    out.spikephases = spikephases;
    out.numspikes_total = numspikes_total;
    out.numspikes_gammal = numspikes_gammal;
else
    out.indices_coherence = [];
    out.fields_indices_coherence = '[<cellindex> , <reftet num> , <coherence b/t cell tetrode and reftet> ]';
    out.spikephases = [];
    out.numspikes_total = [];
    out.numspikes_gammal = [];
end

end


