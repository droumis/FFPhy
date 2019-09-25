function out = calcopenfieldresiduals(f, varargin)
%
%Calculate residuals as compared to entire session firing rate map.
% Use this to see how residuals are correlated with one another, how firing
% rate increases as a function of preference and time in the environment.

binsize = 1;
tbin = 1;
g = gaussian2(5,25);

for option = 1:2:length(varargin)-1
    if isstr(varargin{option})
        switch(varargin{option})
            case 'binsize '
                binsize = varargin{option+1};
            case 'timebin'
                tbin = varargin{option+1};
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
    cellindex(1,:) = []; clear e matches option varargin tmp i
    
    %For each cell, compute the firing rate map, quadrant with peak firing,
    %and residuals
    out(an).type = nan(length(cellindex),1);
    out(an).peak = nan(length(cellindex),1);
    
    for c = 1:size(cellindex,1)
        t = []; x = []; y = []; s = []; q = []; session = []; type = nan(length(epochs),4);
        quadcenters = []; excludetimes = [];
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
                excludetimes = [excludetimes; f(an).excludetime{i}{cellindex(c,(i-1)*4+2)}];

            end
        end
        clear i j tmp
        
        minxy = [min(x) min(y)];
        maxxy = [max(x) max(y)];
        binx = quadcenters(1)-ceil((quadcenters(1)-minxy(1))/binsize)*binsize:binsize:quadcenters(1)+ceil((-quadcenters(1)+maxxy(1))/binsize)*binsize;
        biny = quadcenters(2)-ceil((quadcenters(2)-minxy(2))/binsize)*binsize:binsize:quadcenters(2)+ceil((-quadcenters(2)+maxxy(2))/binsize)*binsize;
        edges{1} = binx; edges{2} = biny; clear minxy maxxy

        occupancy = hist3([x y],'Edges',edges);
        occupancy = occupancy.*median(diff(t));

        spikeind = lookup(s,t);
        spikes = hist3([x(spikeind) y(spikeind)],'Edges',edges);
        rate = conv2(spikes,g,'same')./conv2(occupancy,g,'same');  
        %lowocc = occupancy < 0.05;
        %rate(lowocc) = 0;% clear edges occupancy spikes spikeind
        
        %Determine the quadrant with peak firing
        [peak peakx] = max(rate);
        [peak peaky] = max(peak);
        peakx = binx(peakx(peaky));
        peaky = biny(peaky);
        
        if peakx <= quadcenters(1) && peaky >= quadcenters(2)
            peaklocation = 1;
        elseif peakx>quadcenters(1) && peaky >= quadcenters(2)
            peaklocation = 2;
        elseif peakx>quadcenters(1) && peaky < quadcenters(2)
            peaklocation = 3;
        elseif peakx<= quadcenters(1) && peaky < quadcenters(2)
            peaklocation = 4;
        end
       
        type = type(:,peaklocation); 
        %Determine the cell type
        if all(type==2)
            out(an).type(c) = 0;
        elseif any(type==1)
            out(an).type(c) = find(type==1,1);
        elseif all(type==4)
            out(an).type(c) = -1;
        end
        out(an).peak(c) = peak;
        clear peakx peaky peaklocation quadcenters type peak

        %Determine the timebins
        tbins = t(1):tbin:t(end);
        tbinind = lookup(tbins,t);
        
         
        %Get the number of spikes observed during each time
        [nspikes tmp] = flhist(s,tbins); clear tmp s
        tbins = tbins(2:end);
        tbinind = tbinind(2:end);
        
        %Compute the locations for residuals
        location = [lookup(x(tbinind),binx) lookup(y(tbinind),biny)]; clear x y i invalid
       
        %Get the expected number of spikes based on rate
        ind = sub2ind(size(rate),location(:,1),location(:,2)); clear location
        nexpected = rate(ind)*tbin; clear ind
        
        %add excludetimes
        intersession = find(diff(t)>10);
        excludetimes = [excludetimes; t(intersession) t(intersession+1)];
        excluderate = isExcluded(tbins,excludetimes);
        nexpected(logical(excluderate)) = 0; clear excludetimes excluderate
        
        %Get rid of edge effects
        invalid = [0; diff(session)];
        for i = find(invalid)'
            nexpected(tbins>t(i-1) & tbins<=t(i)) = [];
            nspikes(tbins>t(i-1) & tbins<=t(i)) = [];
            tbinind(tbins>t(i-1) & tbins<=t(i)) = [];
            tbins(tbins>t(i-1) & tbins<=t(i)) = [];
        end
        
        %Compute the residuals
        resid = nspikes(:) - nexpected(:); clear nspikes nexpected ind
        
        out(an).resid{c} = resid;
        out(an).quad{c} = q(tbinind);
        out(an).session{c} = session(tbinind);
        out(an).time{c} = tbins;
        out(an).rate{c} = rate;
        out(an).binx{c} = binx;
        out(an).biny{c} = biny;
        out(an).index{c} = cellindex(c,:);
        clear resid q tbinind session tbins rate binx biny peak type
    end
    
    
    days = min(cellindex(:,2:4:size(cellindex,2)),[],2);
    tmpout = [];
    %Compute residual correlations
    if size(cellindex,1) > 1
        pairsind = nchoosek(1:size(cellindex,1),2);
        validpairs = days(pairsind(:,1)) == days(pairsind(:,2)) & days(pairsind(:,1))~= 0; 
        pairsind = pairsind(validpairs,:); clear days validpairs

        tmpout = nan(size(pairsind,1),12);
            % 1: cell 1
            % 2: cell 2
            % 3: overlap
            % 4: correlation coefficent
            % 5: pvalue
            % 6: resid correlation in familiar session
            % 7: resid correlation in novel session
            % 8: resid correlation in super novel session
            % 9: resid correlation in familiar object quadrant novel sessions
            % 10: resid correlation in novel object quadrant novel sessions
            % 11: cell 1 peak rate location
            % 12: cell 2 peak rate location
        for i = 1:size(pairsind,1)
            % Column 1 & 2: Cell1 Cell2 (index into cindex)
            tmpout(i,[1 2]) = pairsind(i,:);

            if ~isempty(out(an).rate{pairsind(i,1)}) && ~isempty(out(an).rate{pairsind(i,2)})
                % Column3: Cell pair normalized overlap
                rate1 = out(an).rate{pairsind(i,1)};
                rate2 = out(an).rate{pairsind(i,2)};
                if all(size(rate1)==size(rate2))
                    tmpout(i,3) = calc2doverlap(rate1,rate2,'minbins',0);
                else
                    x1 = out(an).binx{pairsind(i,1)}; y1 = out(an).biny{pairsind(i,1)};
                    x2 = out(an).binx{pairsind(i,2)}; y2 = out(an).biny{pairsind(i,2)};

                    if ~isequal(x1,x2)
                        [x a b] = intersect(x1,x2);
                        rate1 = rate1(a,:); rate2 = rate2(b,:);
                    end
                    if ~isequal(y1,y2)
                        [y a b] = intersect(y1,y2);
                        rate1 = rate1(:,a); rate2 = rate2(:,b);
                    end
                    tmpout(i,3) = calc2doverlap(rate1,rate2,'minbins',0);
                    clear x y a b x1 x2 y1 y2
                end
                clear rate1 rate2

                resid1 = out(an).resid{pairsind(i,1)};
                resid2 = out(an).resid{pairsind(i,2)};
                nonzero = (resid1~=0 & resid2~=0);
                if sum(~isnan(resid1(nonzero)) & ~isnan(resid2(nonzero)))*tbin >10
                    % Column 4 & 5: correlation coefficient and pvalue
                    [tmpout(i,4) tmpout(i,5)] = corr(resid1(nonzero),resid2(nonzero), ...
                        'type','Spearman','rows','complete');
                    for s = 1:max(out(an).session{pairsind(i,1)})
                        if sum(~isnan(resid1(nonzero & out(an).session{pairsind(i,1)} == s)) &...
                                ~isnan(resid2(nonzero & out(an).session{pairsind(i,1)} == s)))*tbin > 10
                            % Column 6 - 8: residual correlation by session
                            tmpout(i,5+s) = corr(resid1(nonzero & out(an).session{pairsind(i,1)} == s),...
                                resid2(nonzero & out(an).session{pairsind(i,1)} == s),...
                                'type','Spearman','rows','complete');
                        end
                    end
                    if sum(~isnan(resid1(nonzero & out(an).quad{pairsind(i,1)}==2 & out(an).session{pairsind(i,1)} > 1)) &...
                        ~isnan(resid2(nonzero & out(an).quad{pairsind(i,1)}==2 & out(an).session{pairsind(i,1)} > 1)))*tbin > 10
                        % Column 9: residual correlation familiar object, novel sessions
                        tmpout(i,9) = corr(resid1(nonzero & out(an).quad{pairsind(i,1)}==2 & out(an).session{pairsind(i,1)} > 1),...
                            resid2(nonzero & out(an).quad{pairsind(i,1)}==2 & out(an).session{pairsind(i,1)} > 1),...
                            'type','Spearman','rows','complete');
                    end
                    if sum(~isnan(resid1(nonzero & out(an).quad{pairsind(i,1)}==1 & out(an).session{pairsind(i,1)} > 1)) & ...
                        ~isnan(resid2(nonzero & out(an).quad{pairsind(i,1)}==1 & out(an).session{pairsind(i,1)} > 1)))*tbin > 10
                        % Column 10: residual correlation novel object, novel sessions
                        tmpout(i,10) = corr(resid1(nonzero & out(an).quad{pairsind(i,1)}==1 & out(an).session{pairsind(i,1)} > 1),...
                            resid2(nonzero & out(an).quad{pairsind(i,1)}==1 & out(an).session{pairsind(i,1)} > 1),...
                            'type','Spearman','rows','complete');
                    end
                end
                tmpout(i,11) = out(an).type(pairsind(i,1));
                tmpout(i,12) = out(an).type(pairsind(i,2));
                clear resid1 resid2 nonzero s
            end
        end
    end
    
    out(an).residualcorrelation = tmpout; clear tmpout
    
end
