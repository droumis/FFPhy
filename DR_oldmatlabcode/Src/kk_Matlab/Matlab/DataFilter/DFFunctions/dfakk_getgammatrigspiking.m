function [out] = dfakk_getgammatrigspiking(index, excludeperiods, lowgamma, gammal, spikes, tetinfo,varargin)
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
        otherwise
            error(['Option ''', varargin{option}, ''' not defined']);
    end
end

%% First receive valid lowgamma periods for the epoch's day. %%%

day = index(1,1);
epoch = index(1,2);

if ~isempty(gammatetfilter)
    tetlist =  evaluatefilter(tetinfo{day}{epoch},gammatetfilter); 
else
    error('need to specify gammatetfilter!')
end

if isempty(tetlist)
    disp(sprintf('no valid gammatetfilter tetrodes this epoch'));
    out.cellindices = [];
    out.spikephases = [];
    out.numspikes_total = [];
    out.numspikes_gammal = [];
    out.gammal_durations = [];
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
    out.cellindices = [];
    out.spikephases = [];
    out.numspikes_total = [];
    out.numspikes_gammal = [];
    out.gammal_durations = [];
    return
end

% reconstitute valid gammal into riptimes [start end] list
glstart = times([0 diff(includetimes')]==1);
glend = times(diff(includetimes')==-1);

% throw out last gl event if interrupted by end of epoch
if (glend(end)-glstart(end) < 0)
    glstart = glstart(1:(end-1));
end
% throw out first gl event if occurs at beginning of epoch
if (glend(1)-glstart(1) < 0)
   glend = glend(2:end); 
end

gltimes = [glstart' glend'];

%%%%%%%%


%% Secondly, iterate through each valid gammal period, and for each
    % gammatetfilter tetrode, get # clustered cells, filtered eeg.
    % If specified, select the dominant tetrode as the phase tetrode ("phasetet")
    
    
%initialize output vector
gammalperiods = struct;
    
for p=1:size(gltimes,1)
    
    % initialize
    gammalperiods(p).startgl = gltimes(p,1);            % clock time (fit in 1 ms vector)
    gammalperiods(p).endgl = gltimes(p,2);
    gammalperiods(p).sortedtetrodes = [];      % [<tetrode> , <numcells> , <proportion active>] 
    gammalperiods(p).phasetet = [];
    gammalperiods(p).instantaneous_phase = [];              
    
    % iterate through tetfilter tetrodes (i.e. CA3)
        % 1. Find all the tetrodes that *participate at all* in the gamma episode
                % "participation" = any part of any tet's extracted episodes that
                % is in the gamma episode
        % 2. Then identify the tetrode w/ the longest duration of power exceeding 2 SD
                % This will be the tetrode for which you retrieve phase
                % (phasetet)
        % (3. This shouldn't happen, but .numcells is used to break ties.)
        
    tet_numcells_activeprop = [];
        
    for tet=tetlist'
        tet_events = [ gammal{day}{epoch}{tet}.starttime   gammal{day}{epoch}{tet}.endtime ];
        % check if gamma episode overlaps with any of the tetrodes' extracted events
        if isExcluded(gltimes(p,1),tet_events) || isExcluded(gltimes(p,2),tet_events)
            tet_numcells_activeprop = [tet_numcells_activeprop ; tet tetinfo{day}{epoch}{tet}.numcells];
        end
    end
    
    % error checking
    if isempty(tet_numcells_activeprop)
       disp(sprintf('a gamma event fails to be claimed by any tetrodes? day %d epoch %d',day,epoch));
       continue
    end
    
    % if there is a tie, then iterate through tets and collect eeg
    for kk = 1:size(tet_numcells_activeprop,1)
        
        tet = tet_numcells_activeprop(kk,1);
        
        filteeg = lowgamma{day}{epoch}{tet};
            filteegtimes = geteegtimes(filteeg);
        startind = lookup(gltimes(p,1),filteegtimes);
        endind = lookup(gltimes(p,2),filteegtimes);
        mag_envelope = filteeg.data(startind:endind,3);
        
        % retrieve epoch baseline and std magnitudes
        baseline_mag = gammal{day}{epoch}{tet}.baseline;
        std_mag = gammal{day}{epoch}{tet}.std;
        
        % compute proportion of gammal period for which the tetrode is above 2 SD
        proportion_active =  sum(mag_envelope > (baseline_mag + 2 * std_mag) ) / (endind - startind + 1) ;
        % report % active
        tet_numcells_activeprop(kk,3) = proportion_active;
        
        % if local specified, then output list of times -within- gltimes in which tetrode is above 2 SD 
        if local_flag
            activetimes = vec2list(mag_envelope > (baseline_mag + 2 * std_mag) , filteegtimes);
            gammalperiods(p).activetimes{tet_numcells_activeprop(kk,1)} = activetimes;
        end
        
    end
    
    % install sortedtetrodes
    sortedtetrodes = sortrows(tet_numcells_activeprop,[-3 -2]);
    gammalperiods(p).sortedtetrodes = sortedtetrodes;

    % install phasetet
    if local_flag
        % all participating tetrodes
        gammalperiods(p).phasetet = gammalperiods(p).sortedtetrodes(:,1);
    elseif ~local_flag
        % just the most active one
        gammalperiods(p).phasetet = gammalperiods(p).sortedtetrodes(1,1);
    end
    
    % lastly retrieve lowgamma inst. phase from phasetet
    
    % If 'local' specified, then collect ALL participating tetrodes
    % instantaneous phase -- in the same row as 
    if ~local_flag
        filteeg = lowgamma{day}{epoch}{phasetet};
        eegtimes = geteegtimes(filteeg);
        startind = lookup(gltimes(p,1),eegtimes);
        endind = lookup(gltimes(p,2),eegtimes);
        %install
        gammalperiods(p).instantaneous_phase = filteeg.data(startind:endind,2);
    elseif local_flag
        for jj = 1:size(sortedtetrodes,1)
            filteeg = lowgamma{day}{epoch}{sortedtetrodes(jj,1)};
            eegtimes = geteegtimes(filteeg);
            startind = lookup(gltimes(p,1),eegtimes);
            endind = lookup(gltimes(p,2),eegtimes);
            %install
            gammalperiods(p).instantaneous_phase(jj,:) = filteeg.data(startind:endind,2);        
        end
    end
end

list_of_durations = gltimes(:,2)-gltimes(:,1);


%% Thirdly, assign spikes to phases and output.

% Retrieve the gammal phases of spikes.

% (keep track of coherences between tetrodes in this epoch)

% [ tetrode1 tetrode2 coherence ]
mostcoherent_table = [];

% iterate over units


for i = 1:size(index,1)
    
    filtertimes = [];
    
    if local_flag
        
        % if CA1 (same region as gammafiltertet), choose the cell's own tetrode ("local") to get event times
        if any(index(i,3) == tetlist)
            reftets(i,:) = index(i,[1 2 3]);
            for p=1:size(gammalperiods,2)
                if ~isempty(gammalperiods(p).sortedtetrodes)
                    if any(gammalperiods(p).sortedtetrodes(:,1) == index(i,3))
                        filtertimes = [filtertimes ; gammalperiods(p).activetimes{index(i,3)}];
                    end
                end
            end
        % if some other region (say CA3), identify CA1 tetrode that is most coherent to get events    
        else
            % check if you've previously calculated-- if so, just retrieve
            if ~isempty(mostcoherent_table)
                ind = find(mostcoherent_table(:,1)==index(i,3));
                if ind
                    mostcoherent_tet = mostcoherent_table(ind,2);
                end
                % if not previously calculated, then calculate all coherencies
            else
                gamma1 = lowgamma{index(i,1)}{index(i,2)}{index(i,3)};
                coherences = [];
                for kk = 1:length(tetlist)
                    gamma2 = lowgamma{index(i,1)}{index(i,2)}{tetlist(kk)};
                    coherences = [coherences ; tetlist(kk) meancoherence(gamma1,gamma2,freqrange)];
                end
                coherences = sortrows(coherences,-2);
                mostcoherent_tet = coherences(1,1);
                mostcoherent_table = [mostcoherent_table ; index(i,3) mostcoherent_tet];   
            end
           
            % for reporting later
            reftets(i,:) = [index(i,1) index(i,2) mostcoherent_tet];
            
            % now retrieve event times from the correct tetrode
            for p=1:size(gltimes,1)
                if any(gammalperiods(p).sortedtetrodes == mostcoherent_tet)
                    filtertimes = [filtertimes ; gammalperiods(p).activetimes{mostcoherent_tet}];
                end
            end
            
            disp(sprintf(' tet %d not a eventtet, choosing %d tet for reference based on coherence',index(i,3),mostcoherent_tet));

        end   
    elseif ~local_flag
        filtertimes = gltimes;
    end
    
    % get spikes and their phases
    spiketimes = spikes{index(i,1)}{index(i,2)}{index(i,3)}{index(i,4)}.data(:,1);
    includedspike_times = spiketimes(logical(isExcluded(spiketimes,filtertimes)));
    filteegtimes=geteegtimes(lowgamma{index(i,1)}{index(i,2)}{index(i,3)});
    includedspike_indices = lookup(includedspike_times,filteegtimes);
    
    if local_flag
        if any(index(i,3) == tetlist)
            epochphases = lowgamma{index(i,1)}{index(i,2)}{index(i,3)}.data(:,2);
        else
            epochphases = lowgamma{index(i,1)}{index(i,2)}{mostcoherent_tet}.data(:,2);
        end
    elseif ~local_flag
            % NOT IMPLEMENTED .. would have to get phases for each
            % consensus period..
    end
   
    spikephases{i} = epochphases(includedspike_indices);                            % list of phase times of spikes
    numspikes_gammal(i,:) = length(epochphases(includedspike_indices));             % # spikes firing during gammal period    
    numspikes_total(i,:) = length(spiketimes);                                          % # all spikes fired in the epoch
    cellindices(i,:) = index(i,:);

    
end


%% Lastly send to output.

out.cellindices = cellindices;
out.reftets = reftets;          % corresponds to cellindices
out.spikephases = spikephases;
out.numspikes_total = numspikes_total;
out.numspikes_gammal = numspikes_gammal;
out.gammal_durations = list_of_durations;



end