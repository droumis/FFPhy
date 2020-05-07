function [out] = kk_getconstimes(animaldir,animalprefix, epochs, eventconsname, tetfilter, varargin)
% out = kk_getconstimes(animaldir,animalprefix, epochs, eventconsname, varargin)

%%    Use in conjunction with <event>cons that have been pre-extracted and stored in the animal's day directory..
   %  following Csicsvari--Buzsaki approach to detecting fast oscillation
   %  events
%
%     animaldir and animal prefix are strings indicating the base director for
%     the animal's data and the prefix for the data files
%
%     epochs is an Nx2 list of days and epochs
%
%     tetlist is a list of tetrodes to use or an empty matrix if the
%     'cellfilter' option is used.
%
%	'minthresh', minthresh
%		     specifies a minimum threshold in stdev units for a valid
%			ripple event  (default 0)
%
% Produces a cell structure with a time field and an nripples field which
% indicates the number of electrodes with a ripple at each time point
%
% Examples:
% getriptimes('/data/name/Fre', 'fre', epochs, 1)
% getriptimes('/data/name/Fre', 'fre', epochs, [], 'cellfilter', '(isequal($area, ''CA1''))')

% default values
cellfilter = '';
minenergy = 0;
minthresh = 0;
exclusion_dur = 0;     % duration after ripple within which start of any subsequent ripple is eliminated
exclusion_nevents = 0;    % number of tetrodes that must have ripple for exclusion to kick in
exclusion2_eventconsname = '';
minvelocity = 0;
coherence_flag = 0;

optioninds = [];
for aa = 1:length(varargin)
    if ischar(varargin{aa})
        optioninds = [optioninds aa];
    end
end

for option = optioninds
    switch varargin{option}
        case 'minthresh'
            minthresh = varargin{option+1};
        case 'exclusion_dur'                 % kk implemented 8.1.13 -- exclusion of chained events + enforces minimum # of tetrodes participating
            exclusion_dur = varargin{option+1};     % duration after ripple within which start of any subsequent ripple is eliminated
        case 'consensus_numtets'         
            consensus_numtets = varargin{option+1};  
        case 'exclusion2'      % kk implemented 7.30.14, excludes events that overlap with primary event
            exclusion2_eventconsname = varargin{option+1};
            exclusion2_tetfilter = varargin{option+2};
            exclusion2_consensus_numtets = varargin{option+3};
            exclusion2_minthresh = varargin{option+4};
            exclusion2_dur = varargin{option+5};
            exclusion2_maxvelocity = varargin{option+6};
        case 'minvelocity'
            minvelocity = varargin{option+1};
        case 'maxvelocity'    % use this as a correct alternative to get2dstate .. so you don't get truncated filtered times..
            maxvelocity = varargin{option+1};
        case 'coherence'
            coherence_flag = 1;
            coh_percentile = varargin{option+1};
            coh_pairnums = varargin{option+2};
    end
end

%check to see if a cell filter is specified
if ~isempty(cellfilter)
    % this will cause us to ignore tetlist
    cellinfo = loaddatastruct(animaldir, animalprefix, 'cellinfo');
end

loaddays = unique(epochs(:,1));

tetinfo = loaddatastruct(animaldir, animalprefix, 'tetinfo');
eventcons = loaddatastruct(animaldir, animalprefix, eventconsname, loaddays);
if exist('maxvelocity','var') ||  ( minvelocity > 0 )
   pos = loaddatastruct(animaldir, animalprefix, 'pos', loaddays);
end

for i = 1:size(epochs,1)
    
        %disp(sprintf('kk_getconstimes: (START) %s %d %d',animalprefix,epochs(i,1),epochs(i,2)))
    
        % empty error checking: possibly no eventcons..(perhaps not computed because of too much noise this epoch) if so, skip
        if isempty(eventcons)
            disp(sprintf('d %d e %d skipped, eventcons is empty',epochs(i,1),epochs(i,2)))
            out{epochs(i,1)}{epochs(i,2)}.time = nan;
            out{epochs(i,1)}{epochs(i,2)}.cons = nan;  
            continue
        elseif (epochs(i,2) > length(eventcons{epochs(i,1)})) || isempty(eventcons{epochs(i,1)}{epochs(i,2)})
            disp(sprintf('d %d e %d skipped, no eventcons data',epochs(i,1),epochs(i,2)))
            out{epochs(i,1)}{epochs(i,2)}.time = nan;
            out{epochs(i,1)}{epochs(i,2)}.cons = nan;            
            continue
        end
        
        % Apply tetfilter to find which entry (grouping of tetrodes) in eventcons is the right one.
            % If a number, then it's the index.
            % If a string, then search for the TF (tetfilter) index within eventcons
        if isnumeric(tetfilter)
            TF = tetfilter;
        else
            TF = 0;
            for ttt = 1:length(eventcons{epochs(i,1)}{epochs(i,2)})
                if ~isempty(eventcons{epochs(i,1)}{epochs(i,2)}{ttt})
                    if strcmp(eventcons{epochs(i,1)}{epochs(i,2)}{ttt}.tetfilter,tetfilter)
                        TF = ttt;
                    end
                end
            end
            if ~TF
                error(sprintf('d %d e %d skipped, **could not find matching eventcons for tetfilter',epochs(i,1),epochs(i,2)))
            end
        end
           
        ec = eventcons{epochs(i,1)}{epochs(i,2)}{TF};
        
        % safety check: should have already excluded this epoch with your epochfilter
        % if you don't have minimum number of detecting consensus tetrodes this epoch
        if length(ec.tetlist) < consensus_numtets
            disp('*************eventcons data should have at least minimum number of tets')
            out{epochs(i,1)}{epochs(i,2)}.time = 0;
            out{epochs(i,1)}{epochs(i,2)}.cons = 0;                
            continue
        end        
        
        % choose your output times vector to be the same as the eegtimesvec
        % of the first participating tetrode
        times = ec.timerange(1):1/ec.samprate:ec.timerange(end);
        
        % Filter for events of threshold size
        evalid = (ec.maxthresh > minthresh);
        
        %%% If specified, filter for events based on coherence in this band %%%%%
            %  take minimum percentile coherence of event for specified
                % tetrode pair (priority given by pair number vector)
        if coherence_flag
            if ~isfield(ec,'coherence')
                disp(sprintf('d %d e %d skipped, coherence_flag ==1 but no .coherence field in eventcons',epochs(i,1),epochs(i,2)))
                out{epochs(i,1)}{epochs(i,2)}.time = nan;
                out{epochs(i,1)}{epochs(i,2)}.cons = nan;
                continue
            end
            % look for valid tetrode pairs (ordered by priority), then
            % calculate coherence_cutoff
            coh_cutoff = nan;
            pairnum = nan;
            for pp = coh_pairnums
                if (pp > length(ec.coherence))  ||  isempty(ec.coherence{pp})
                    continue
                else
                    pairnum = pp;
                    coh_cutoff = prctile(ec.coherence{pp}.eventcoherence,coh_percentile);                   
                end
            end
            % use cutoff to find qualifying events
            if ~isnan(coh_cutoff)
                goodeventinds = (ec.coherence{pairnum}.eventcoherence > coh_cutoff); 
            else               
                disp(sprintf('d %d e %d skipped, coherence_flag == 1 but no available pairs',epochs(i,1),epochs(i,2)))
                out{epochs(i,1)}{epochs(i,2)}.time = nan;
                out{epochs(i,1)}{epochs(i,2)}.cons = nan;
                continue            
            end
            
            % filter events
            evalid = evalid & goodeventinds;
                 
        end          
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        % list form of the cons "events"
        ectimes = [ec.starttime(evalid) ec.endtime(evalid)];
            
        % possibly no events.. if so skip (kk added 7.29.13)
        if isempty(ectimes)
            out{epochs(i,1)}{epochs(i,2)}.time = 0;
            out{epochs(i,1)}{epochs(i,2)}.cons = 0;                
            continue
        end
        
        % Select events that occur between min- and maxvelocity (2D velocity)
        if exist('maxvelocity','var')
            posentry = pos{epochs(i,1)}{epochs(i,2)};
            posinds = lookup(ectimes(:,1),pos{epochs(i,1)}{epochs(i,2)}.data(:,1));  % velocity at starttime
            if size(posentry.data,2) > 5
                goodinds = posentry.data(posinds,9) < maxvelocity & ...
                           posentry.data(posinds,9) >= minvelocity ;
            else
                goodinds = posentry.data(posinds,5) < maxvelocity & ...
                           posentry.data(posinds,5) >= minvelocity ;
                disp(sprintf('%s : using column 5 velocity..',animalprefix))
            end
            ectimes = ectimes(goodinds,:);
            conscount = sum(goodinds); % to report out later
        end
        
        if isempty(ectimes)
            out{epochs(i,1)}{epochs(i,2)}.time = times;
            out{epochs(i,1)}{epochs(i,2)}.cons = zeros(size(times));
            disp(sprintf('d%de%d no consensus times!',epochs(i,1),epochs(i,2)))
            clear times;
            continue
        end
        
        %check for possible last ripple that extends past times vector
        if times(end)-ectimes(end,2) < 0
            ectimes(end,:) = [];
            disp(sprintf('excluded last event, day %d epoch %d',epochs(i,1),epochs(i,2)))
        end
            
        % If exclusion2 is specified (i.e. want exclude ripplescons from gammafcons events),
            % To do this properly in conjunction with exclusion of chained
            % events, we follow this order:
                   % 1. Flag ectimes rows that have any overlap with the output of exclusion 2   (i.e. exclude gammaf events that overlap with ripples)
                   % 2. Below, enforce exclusion_dur, eliminating ectimes events that follow too soon after a previous one
                   % 3. Lastly, delete the  ectimes events flagged in step 1
        if ~isempty(exclusion2_eventconsname)
            % First obtain exclusion2 output via a recursive call
            ex2 = kk_getconstimes(animaldir,animalprefix,epochs(i,:),exclusion2_eventconsname,exclusion2_tetfilter,...
                                 'consensus_numtets',exclusion2_consensus_numtets,'minthresh',exclusion2_minthresh,...
                                 'exclusion',exclusion2_dur,'maxvelocity',exclusion2_maxvelocity);
            if isnan(ex2{epochs(i,1)}{epochs(i,2)}.time)   % ripples not available for exclusion -- shut down this epoch by setting exclude vec to 1s
                excludetimevec = times;
                excludecons = ones(size(times));
            elseif length(ex2{epochs(i,1)}{epochs(i,2)}.time) == 1 % no ripples detected this epoch..
                excludetimevec = times;
                excludecons = zeros(size(times));    
            else                % normal epoch w/ ripples
                excludecons = ex2{epochs(i,1)}{epochs(i,2)}.cons;
                excludetimevec = ex2{epochs(i,1)}{epochs(i,2)}.time;                
            end
            % Flag the ectimes entries (rows) that have any overlap with the exout.cons
            ectimes = [ectimes zeros(size(ectimes,1),1)];    % first append a column to use for flagging
            for ee = 1:size(ectimes,1)
               starttime = ectimes(ee,1);
               endtime = ectimes(ee,2);
               % convert the single consensus event to vector format to detect overlap
               singlevec = list2vec([starttime endtime],excludetimevec)';
               if any(singlevec & excludecons)
                   ectimes(ee,3) = 1;  % flag entire event
                   %disp('exclusion2 triggered .. coincident event')
               end
            end
        end
        
        % Filter out events that occur within exclusion_dur
        if exclusion_dur > 0
            for rr=fliplr(2:size(ectimes,1))
                timesinceprevevents = ectimes(rr,1) - ectimes(rr-1,2);
                if timesinceprevevents < exclusion_dur
                    ectimes(rr,:) = [ ];
                    %disp(sprintf('kk_geteventtimes: chain-excluded event: d %d e %d (%d ms)',epochs(i,1),epochs(i,2),round(timesinceprevevents*1000)))
                end
            end
        end
        
        % If exclusion2 was specified, then remove previously flagged events
        if ~isempty(exclusion2_eventconsname)
            flaggedinds = (ectimes(:,3) == 1);
            ectimes(flaggedinds,:) = [];
            ectimes(:,3) = [];
        end
        
        
        % print out number of surviving events
        %disp(sprintf('d %d e %d: %d out of %d original %s events remain',...
        %    epochs(i,1),epochs(i,2),size(ectimes,1),conscount,eventconsname))
        
        % convert from list form back to vector form
        cons = list2vec(ectimes,times)';
        
        out{epochs(i,1)}{epochs(i,2)}.time = times;
        out{epochs(i,1)}{epochs(i,2)}.cons = cons;
        clear times
        clear cons
        clear ectimes
        clear ec
        
        %disp(sprintf('kk_getconstimes: (END)',animalprefix,epochs(i,1),epochs(i,2)))
        
end



end