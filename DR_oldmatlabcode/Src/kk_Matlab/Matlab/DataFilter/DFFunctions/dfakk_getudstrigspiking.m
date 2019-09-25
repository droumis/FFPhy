function [out] = dfakk_getudstrigspiking(index, excludeperiods, uds, spikes, tetinfo,varargin)
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

%   'spindletetfilter' -- these tetrodes will be used to extract gamma phase
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
reftet = [];
window = [.2 .2];
numtets = 1;

for option = 1:2:length(varargin)-1
    switch varargin{option}
        case 'udstetfilter'
            udstetfilter = varargin{(option+1)};   % tetrodes used to extract the phase
        case 'minthresh'
            spindle_thresh = varargin{(option+1)};   % tetrodes used to extract the phase           
        case 'window'                                   
            window = varargin{(option+1)};   % window of the UD DU psth      
        case 'numtets'                                   
            numtets = varargin{(option+1)};   % # of tetrodes that need to report the ud/du time  (matter of whether have valid sws period at that time)    
        case 'binsize'                                   
            binsize = varargin{(option+1)};   % # of tetrodes that need to report the ud/du time  (matter of whether have valid sws period at that time)    
        otherwise
            error(['Option ''', varargin{option}, ''' not defined']);
    end
end

%% First find de novo here valid spindle periods for each spindletetfilter tetrode.
    % Then take intersection of valid lowgamma periods with includetimes.

day = index(1,1);
epoch = index(1,2);

if ~isempty(udstetfilter)
    tetlist =  evaluatefilter(tetinfo{day}{epoch},udstetfilter)'; 
else
    error('need to specify udstetfilter!')
end

if isempty(tetlist)
    disp(sprintf('no valid udstetfilter tetrodes this epoch'));
    out.cellindices = index;
    out.ctxtet = [];
    out.numspikes_total = [];
    out.numspikes_uds = [];
    return
end


ud = uds{day}{epoch}{tetlist(1)};
times = ud.timerange(1):0.001:ud.timerange(end);         % 1 ms step vector

%obtain includetimes
includetimes = ~isExcluded(times, excludeperiods);     % 1s and 0s: ones = included, zeros = excluded

% (if there are no valid periods this epoch, report no outputs) 
if sum(includetimes)==0
    disp(sprintf('no includetimes, day %d epoch %d',day,epoch));
    disp(sprintf('%d tetrodes this epoch',length(tetlist)));
    out.ctxtet = [];
    out.numspikes_total = [];
    out.numspikes_uds = [];
    return
end

% get all UD and DU times this epoch

udtimes_all = [];
dutimes_all = [];

    % retrieve all tetrodes' times
for tet = tetlist
    udtimes_all = [udtimes_all ; uds{index(1)}{index(2)}{tet}.ud_times];
    dutimes_all = [dutimes_all ; uds{index(1)}{index(2)}{tet}.du_times];
end    
    % remove events out of timefilter
udtimes_all = isExcluded(udtimes_all,includetimes);
dutimes_all = isExcluded(dutimes_all,includetimes);
    % remove events with less than consensus # of tetrodes
udtimes = unique(udtimes_all);
dutimes = unique(dutimes_all);

    for k=size(udtimes,1):-1:1
        if sum(udtimes_all == udtimes(k)) < numtets
            udtimes(k) = [];
        end
    end
    for k=size(dutimes,1):-1:1
        if sum(dutimes_all == dutimes(k)) < numtets
            dutimes(k) = [];
        end
    end    

dutimes;
udtimes;

% reportage
disp(sprintf('d%d e%d: DU%d UD%d',index(1),index(2),length(dutimes),length(udtimes)))

%%%%%%%%%%%%%%%%%%%%%%%%

% iterate over units
for i = 1:size(index,1)
           
    % get spikes and their phases
    if ~isempty(spikes{index(i,1)}{index(i,2)}{index(i,3)}{index(i,4)})
        spiketimes = spikes{index(i,1)}{index(i,2)}{index(i,3)}{index(i,4)}.data(:,1);
        % histc bin edges
        psthtime = (-window(1)-0.5*binsize):binsize:(window(2)+0.5*binsize);
        % calculate psths
            % DU
        if ~isempty(dutimes)
            dupsth = [];
            % iterate through each ripple and calculate psth
            for d=1:length(dutimes)
                onehist = histc(spiketimes , dutimes(d) + psthtime);
                dupsth = [dupsth ; onehist'];
            end
        end
            % UD
        if ~isempty(udtimes)
            udpsth = [];
            % iterate through each ripple and calculate psth
            for d=1:length(udtimes)
                onehist = histc(spiketimes , udtimes(d) + psthtime);
                udpsth = [udpsth ; onehist'];
            end
        end
        
        dupsths{i} = dupsth;        % list of phase times of spikes
        dupsths{i} = udpsth;
        
        numspikes_ud(i,:) = length(sum(sum(udpsth))); % # spikes histogrammed for UD
        numspikes_du(i,:) = length(sum(sum(dupsth)));
        numspikes_total(i,:) = length(spiketimes);                                       % # all spikes fired in the epoch
        
        cellindices(i,:) = index(i,:);
        
    else
        
        cellindices(i,:) = index(i,:);        
        %disp(sprintf('d%de%dt%dc%d, cell in iterator but empty spike data cell this epoch..',index(i,1),index(i,2),index(i,3),index(i,4)))
        continue
        
    end
    
end


%% Lastly send to output.

if exist('spikephases')
    out.cellindices = cellindices;
    out.dupsths = dupsths;
    out.udpsths = udpsths;
    out.numspikes_total = numspikes_total;
    out.numspikes_ud = numspikes_ud;
    out.numspikes_du = numspikes_du;
end
end


