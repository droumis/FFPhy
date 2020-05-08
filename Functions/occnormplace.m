function out = occnormplace(allpos,eventpos,binsize,timestep,std)

% calculates occupancy & spike/event occnorm place map

    % inputs:

    % allpos : (x,y) list of all positions
    % eventpos : (x,y) list of all events / spikes
    % binsize : spatial binsize to grid the positions
    % timestep : times between successive position samples  
    % std : std of gaussian smoothing kernel
    % nonoccnormflag : set to 1 if do not want to occupancy normalize!


tmpposition = allpos;
tmpevent = eventpos;

if ~isempty(tmpposition)
    
    % initialize outputs in order of their obtaining
    
    % basic position
    out.binx = [] ;
    out.biny = [] ;
    out.occupancy = [] ;
    % events
    out.events = [] ;
    out.events_smoothed = [];
    % occupancy normalized
    out.eventrate = [] ;
    out.eventrate_smoothed = [] ;
    out.peakrate_raw = nan ;
    out.peakrate_smoothed = nan ;
    out.eventrate_smoothed_norm = [] ;
    % occupancy
    out.occupancy_smoothed  =   [] ;
    
    minx = floor(min(tmpposition(:,1)));
    maxx = ceil(max(tmpposition(:,1)));
    binx = (minx:binsize:maxx);
    miny = floor(min(tmpposition(:,2)));
    maxy = ceil(max(tmpposition(:,2)));
    biny = (miny:binsize:maxy);
    
    out.binx = binx;
    out.biny = biny;
    
    out.occupancy = fasthist2(tmpposition(:,1), tmpposition(:,2), binx, biny);

    if ~isempty(tmpevent)
        
        out.events = fasthist2(tmpevent(:,1), tmpevent(:,2), binx, biny);

        nonzero = find(out.occupancy > 0);
        out.eventrate = zeros(size(out.events));
        out.eventrate_smoothed = zeros(size(out.occupancy));
       
        if isempty(find(out.events > 0))
           out.eventrate_smoothed_norm = [];
        else
            % occupancy-normalize bins containing events
            out.eventrate(nonzero) = out.events(nonzero) ./ ( timestep * out.occupancy(nonzero) ) ;

            % smooth events and eventrate w/ Gaussian filter
            g = gaussian2(std,(6*std));         % is this the right filter?
            out.events_smoothed = filter2(g,out.events); 
            out.eventrate_smoothed = filter2(g,out.eventrate);
            
            % normalize smoothedeventrate map by dividing by peak rate
            out.eventrate_smoothed_norm = out.eventrate_smoothed / max(out.eventrate_smoothed(:)) ;
 
            % also, smooth occupancy w/ the same Gaussian filter
            out.occupancy_smoothed = filter2(g, out.occupancy);
            
        end
        
    end
    
end

end