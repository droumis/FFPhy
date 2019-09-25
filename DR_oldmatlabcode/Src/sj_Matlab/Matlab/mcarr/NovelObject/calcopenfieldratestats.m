function out = calcopenfieldratestats(f, varargin)
%
%Calculate mean firing rate, place field location, and place field correlation
% as compared to entire session firing rate map.

binsize = 1;
g = gaussian2(5,25);

for option = 1:2:length(varargin)-1
    if isstr(varargin{option})
        switch(varargin{option})
            case 'binsize'
                binsize = varargin{option+1};
            case 'filter'
                g = varargin{option+1};
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end
    else
        error('Options must be strings, followed by the variable');
    end
end

for an = 1:length(f)
    %Determine which cells go together
    epochs = 1:length(f(an).output);
    cellindex = nan(1,4*length(epochs));
    for i = epochs
        for e = 1:length(f(an).data{epochs(i)})
            for c = 1:size(f(an).data{epochs(i)}{e},1)
                tmp = [epochs(i) e f(an).data{epochs(i)}{e}(c,:)];
                if i == 1
                    cellindex(end+1,1:4) = tmp;
                else
                    matches = rowfind(tmp(:,[2 3 4]),cellindex(:,[2 3 4]));
                    if matches~=0
                        cellindex(matches,(i-1)*4+1:(i-1)*4+4) = tmp;
                    else
                        matches = rowfind(tmp(:,[2 3 4]),cellindex(:,[6 7 8]));
                        if matches~=0
                            cellindex(matches,(i-1)*4+1:(i-1)*4+4) = tmp;
                        else
                            cellindex(end+1,(i-1)*4+1:(i-1)*4+4) = tmp;
                        end
                    end
                end
            end
        end
    end
    
    if ~isempty(cellindex)
        cellindex(1,:) = []; clear e matches option varargin tmp i
    
        %Get rid of cells that don't fire in the familiar session
        cellindex(cellindex(:,1)==0,:) = [];
        
        %Get rid of cells that don't fire in either novel or supernovel
        %sessions
        invalid = find(sum(cellindex'==0)==4*(length(epochs)-1));
        cellindex(invalid,:) = [];

    
    %For each cell, compute the firing rate map, location with peak firing relative to novelty,
    %mean firing by session, location of center of mass relative to novelty, correlation between 
    %novel sessions and familiar sessions firing rate maps

    out(an).type = nan(length(cellindex),1);
    out(an).peak = nan(length(cellindex),1);
    out(an).mean = nan(length(cellindex),2);
    out(an).corr = nan(length(cellindex),2);
    out(an).location = cell(length(cellindex),1);
    
    for c = 1:size(cellindex,1)
        t = []; x = []; y = []; s = []; q = []; session = []; type = nan(length(epochs),4);
        quadcenters = [];
        for i = epochs
            if any(cellindex(c,(i-1)*4+1:(i-1)*4+4)) ~= 0
                t = [t; f(an).output{cellindex(c,(i-1)*4+1)}(cellindex(c,(i-1)*4+2)).time];
                x = [x; f(an).output{cellindex(c,(i-1)*4+1)}(cellindex(c,(i-1)*4+2)).x];
                y = [y; f(an).output{cellindex(c,(i-1)*4+1)}(cellindex(c,(i-1)*4+2)).y];
                q = [q; f(an).output{cellindex(c,(i-1)*4+1)}(cellindex(c,(i-1)*4+2)).quad];
                ind = rowfind([cellindex(c,(i-1)*4+3) cellindex(c,(i-1)*4+4)],f(an).output{cellindex(c,(i-1)*4+1)}(cellindex(c,(i-1)*4+2)).index(:,[3 4]));
                s = [s; f(an).output{cellindex(c,(i-1)*4+1)}(cellindex(c,(i-1)*4+2)).spikes{ind}];
                session = [session; i*ones(size(f(an).output{cellindex(c,(i-1)*4+1)}(cellindex(c,(i-1)*4+2)).time))];
                if isempty(quadcenters)
                    quadcenters = f(an).output{cellindex(c,(i-1)*4+1)}(cellindex(c,(i-1)*4+2)).quadcenters;
                end
                
                tmp = f(an).output{cellindex(c,(i-1)*4+1)}(cellindex(c,(i-1)*4+2)).types;
                for j = 1:4
                    %Determine what type of quadrant j is
                    if tmp(1,j) && tmp(2,j);
                        type(i,j) = 1;
                    elseif tmp(1,j) && ~tmp(2,j);
                        type(i,j) = 2;
                    elseif ~tmp(1,j) && tmp(2,j);
                        type(i,j) = 3;
                    elseif ~tmp(1,j) && ~tmp(2,j);
                        type(i,j) = 4;
                    else
                        type(i,j) = NaN;
                    end   
                end
            end
        end
        clear i j tmp
        
        %Compute meanrate by session
        spike_sessions = session(lookup(s,t));
        famrate = sum(spike_sessions==1)./(sum(session==1).*median(diff(t)));
        novrate = sum(spike_sessions==2)./(sum(session==2).*median(diff(t)));
        suprate = sum(spike_sessions==3)./(sum(session==3).*median(diff(t)));
        out(an).mean(c,:) = [novrate./famrate suprate./famrate];
        clear famrate novrate suprate
        
        %Determine bins for place field computations
        minxy = [min(x) min(y)];
        maxxy = [max(x) max(y)];
        binx = quadcenters(1)-ceil((quadcenters(1)-minxy(1))/binsize)*binsize:binsize:quadcenters(1)+ceil((-quadcenters(1)+maxxy(1))/binsize)*binsize;
        biny = quadcenters(2)-ceil((quadcenters(2)-minxy(2))/binsize)*binsize:binsize:quadcenters(2)+ceil((-quadcenters(2)+maxxy(2))/binsize)*binsize;
        edges{1} = binx; edges{2} = biny; clear minxy maxxy

        %Calculate the correlation between place field maps by session
        occupancy = hist3([x(session==1) y(session==1)],'Edges',edges);
        occupancy = occupancy.*median(diff(t));
        spikeind = lookup(s(spike_sessions==1),t);
        spikes = hist3([x(spikeind) y(spikeind)],'Edges',edges);
        ratef = conv2(spikes,g,'same')./conv2(occupancy,g,'same');
        ratef(occupancy==0) = NaN;
        clear occupancy spikeind spikes
        
        if sum(spike_sessions==2)>0
            occupancy = hist3([x(session==2) y(session==2)],'Edges',edges);
            occupancy = occupancy.*median(diff(t));
            spikeind = lookup(s(spike_sessions==2),t);
            spikes = hist3([x(spikeind) y(spikeind)],'Edges',edges);
            raten = conv2(spikes,g,'same')./conv2(occupancy,g,'same');
            raten(occupancy==0) = NaN;
            clear occupancy spikeind spikes
        else
            raten = nan(size(ratef));
        end
        
        if sum(spike_sessions==3)>0
            occupancy = hist3([x(session==3) y(session==3)],'Edges',edges);
            occupancy = occupancy.*median(diff(t));
            spikeind = lookup(s(spike_sessions==3),t);
            spikes = hist3([x(spikeind) y(spikeind)],'Edges',edges);
            rates = conv2(spikes,g,'same')./conv2(occupancy,g,'same');
            rates(occupancy==0) = NaN;
            clear occupancy spikeind spikes
        else
            rates = nan(size(ratef));
        end
        
        tmp = [reshape(raten,size(raten,1)*size(raten,2),1) ...
            reshape(ratef,size(ratef,1)*size(ratef,2),1)];
        invalid = isnan(tmp(:,1)) | isnan(tmp(:,2));
        if sum(invalid,1) ~= size(tmp,1)
            out(an).corr(c,1) = corr(tmp(~invalid,1),tmp(~invalid,2));
        end    
        tmp = [reshape(rates,size(rates,1)*size(rates,2),1) ...
            reshape(ratef,size(ratef,1)*size(ratef,2),1)];
        invalid = isnan(tmp(:,1)) | isnan(tmp(:,2));
        if sum(invalid,1) ~= size(tmp,1)
            out(an).corr(c,2) = corr(tmp(~invalid,1),tmp(~invalid,2));
        end   
        clear ratef raten rates tmp invalid
        
        %Compute global place field to determine location
        occupancy = hist3([x y],'Edges',edges);
        occupancy = occupancy.*median(diff(t));

        spikeind = lookup(s,t);
        spikes = hist3([x(spikeind) y(spikeind)],'Edges',edges);
        rate = conv2(spikes,g,'same')./conv2(occupancy,g,'same');  
        
        %Determine the location with peak firing
        [peak peakx] = max(rate);
        [peak peaky] = max(peak);
        peakx = binx(peakx(peaky));
        peaky = biny(peaky);
        out(an).peak(c) = peak;
         
        %Figure out the location of place field COM
        val = rate(~isnan(rate));
        [comx comy] = find(~isnan(rate));
        sx = sum(val.*binx(comx)')./sum(val);
        sy = sum(val.*biny(comy)')./sum(val); clear comx comy rate
        
        %Determine the peak and COM location relative to familiar object,
        %novel object, and supernovel objects
        location = zeros(4,2);
        objcenters = [median(binx(binx<quadcenters(1))) median(binx(binx>quadcenters(1))) ...
            median(biny(biny<quadcenters(1))) median(biny(biny>quadcenters(2)))];
        location(1,:) = [sqrt((peakx-objcenters(1)).^2 + (peaky-objcenters(4)).^2) sqrt((sx-objcenters(1)).^2 + (sy-objcenters(4)).^2)];
        location(2,:) = [sqrt((peakx-objcenters(2)).^2 + (peaky-objcenters(4)).^2) sqrt((sx-objcenters(2)).^2 + (sy-objcenters(4)).^2)];
        location(3,:) = [sqrt((peakx-objcenters(2)).^2 + (peaky-objcenters(3)).^2) sqrt((sx-objcenters(2)).^2 + (sy-objcenters(3)).^2)];
        location(4,:) = [sqrt((peakx-objcenters(1)).^2 + (peaky-objcenters(3)).^2) sqrt((sx-objcenters(1)).^2 + (sy-objcenters(3)).^2)];
       
        out(an).location{c} = nan(3,2);
        try
            out(an).location{c}(1,:) = location(find(any(type(2:end,:)==2)),:);
            if length(find(any(type(2:end,:)==2)))>1
                out(an).location{c}(1,:) = location(find(type(2,:)==2));
            end
        end
        if sum(spike_sessions==2)>0
            out(an).location{c}(2,:) = location(find(type(2,:)==1),:);
        end
        if sum(spike_sessions==3)>0
            out(an).location{c}(3,:) = location(find(type(3,:)==1),:);
        end
        type = type(find(location(:,1)==min(location(:,1))));
        
        %Determine the cell type
        if all(type==2)
            out(an).type(c) = 0;
        elseif any(type==1)
            out(an).type(c) = find(type==1,1);
        elseif all(type==4)
            out(an).type(c) = -1;
        end
        clear peakx peaky peaklocation quadcenters objcenters type peak q t x y s session

        out(an).index{c} = cellindex(c,:);
    end
    end        
end