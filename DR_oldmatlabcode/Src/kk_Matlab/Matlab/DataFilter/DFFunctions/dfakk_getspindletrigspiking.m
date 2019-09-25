function [out] = dfakk_getspindletrigspiking(index, excludeperiods, spindle, spindles, spikes, tetinfo,varargin)
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

for option = 1:2:length(varargin)-1
    switch varargin{option}
        case 'spindletetfilter'
            spindletetfilter = varargin{(option+1)};   % tetrodes used to extract the phase
        case 'minthresh'
            spindle_thresh = varargin{(option+1)};   % tetrodes used to extract the phase    
        case 'reftet'                                   
            reftet = varargin{(option+1)};   % set manually .. if not want maximum CV + consensus procedure (below)                
        otherwise
            error(['Option ''', varargin{option}, ''' not defined']);
    end
end

%% First find de novo here valid spindle periods for each spindletetfilter tetrode.
    % Then take intersection of valid lowgamma periods with includetimes.

day = index(1,1);
epoch = index(1,2);

if ~isempty(spindletetfilter)
    tetlist =  evaluatefilter(tetinfo{day}{epoch},spindletetfilter)'; 
else
    error('need to specify spindletetfilter!')
end

if isempty(tetlist)
    disp(sprintf('no valid spindletetfilter tetrodes this epoch'));
    out.cellindices = index;
    out.reftet = [];
    out.spikephases = [];
    out.numspikes_total = [];
    out.numspikes_spindles = [];
    out.spindle_durations = [];    
    return
end


sp = spindles{day}{epoch}{tetlist(1)};
times = sp.timerange(1):0.001:sp.timerange(end);         % 1 ms step vector

%obtain includetimes
includetimes = ~isExcluded(times, excludeperiods);     % 1s and 0s: ones = included, zeros = excluded

% (if there are no valid periods this epoch, report no outputs) 
if sum(includetimes)==0
    disp(sprintf('no includetimes, day %d epoch %d',day,epoch));
    disp(sprintf('%d tetrodes this epoch',length(tetlist)));
    out.spikephases = [];
    out.numspikes_total = [];
    out.numspikes_spindles = [];
    return
end

% get all 2 SD (or minthresh) periods for each spindletetfilter tetrode,
% take intersection w/ includetimes

spindleperiods = {};

for tet = tetlist
    
    eventperiods = [spindles{day}{epoch}{tet}.starttime spindles{day}{epoch}{tet}.endtime];
    % apply minthresh criterion
    eventperiods = eventperiods(spindles{day}{epoch}{tet}.spindle_mean_rms_zscore > spindle_thresh,:);
        eventvec = list2vec(eventperiods,times);    
    % take intersection with permissible times
    spindlevec{tet} = eventvec & includetimes;
    % report how much of epoch is 2 SD gamma
    spindleperiods{tet} = vec2list(spindlevec{tet},times);
    percentepoch = round(100*sum(spindleperiods{tet}(:,2)-spindleperiods{tet}(:,1))/(times(end)-times(1)));
        %disp(sprintf('d %d e %d t %d, %d%% spindle',day,epoch,tet,percentepoch));
    
end

%%%%%%%%


%% Then assign spikes to gamma phases, then output.


% Find the reference tet (both event and phase) to be a CTX (via
% spindletetfilter) tet w/ the largest rms CV that day.
if isempty(reftet)
    maxcv = 0;
    for tet = tetlist
        cvs = [];
        for ep=1:length(spindles{day})
            if ~isempty(spindles{day}{ep})
                cvs = [ cvs ; spindles{day}{ep}{tet}.spindle_rms_std/spindles{day}{ep}{tet}.spindle_rms_mean ];
            end
        end
        meancv = mean(cvs);
        if meancv > maxcv
            reftet = tet;
            maxcv = meancv;
        end
    end
end


% Refine the reference tetrodes valid spindle periods:
    % Specifically, look for consensus with at least ONE other ctx tetrode.

consensusvec = zeros(size(spindlevec{tet}));

for tet = tetlist
    if tet ~= reftet
        consensusvec =  ( consensusvec | spindlevec{tet} ) ;
    end
end

finalvec = spindlevec{reftet} & consensusvec ;
finalperiods = vec2list(finalvec,times);
    % reportage
    totalduration = round(sum(finalperiods(:,2)-finalperiods(:,1)));
    totalduration_reftetonly = round(sum(spindleperiods{reftet}(:,2)-spindleperiods{reftet}(:,1)));
        disp(sprintf('d %d e %d: reftet is tet %d, %d s spindle (consensus), by itself %d s',day,epoch,reftet,totalduration,totalduration_reftetonly));
    

% iterate over units
for i = 1:size(index,1)
           
    % get spikes and their phases
    if ~isempty(spikes{index(i,1)}{index(i,2)}{index(i,3)}{index(i,4)})
        spiketimes = spikes{index(i,1)}{index(i,2)}{index(i,3)}{index(i,4)}.data(:,1);
            includedspike_times = spiketimes(logical(isExcluded(spiketimes,finalperiods)));
        filteegtimes=geteegtimes(spindle{index(i,1)}{index(i,2)}{reftet});
        includedspike_indices = lookup(includedspike_times,filteegtimes);
        reftetphases = spindle{index(i,1)}{index(i,2)}{reftet}.data(:,2);
    
        spikephases{i} = reftetphases(includedspike_indices);                            % list of phase times of spikes
        numspikes_spindle(i,:) = length(includedspike_indices);                          % # spikes firing during gammal period
        numspikes_total(i,:) = length(spiketimes);                                       % # all spikes fired in the epoch
        disp(sprintf('%d %d %d %d spindle spikes: %d, ',index(i,1),index(i,2),index(i,3),index(i,4),numspikes_spindle(i,:)))
        
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
    out.reftet = reftet;
    out.spikephases = spikephases;
    out.numspikes_total = numspikes_total;
    out.numspikes_spindles = numspikes_spindle;
    out.spindle_durations = finalperiods(:,2) - finalperiods(:,1);  
end
end


