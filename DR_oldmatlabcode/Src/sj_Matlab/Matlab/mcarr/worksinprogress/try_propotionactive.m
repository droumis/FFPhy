%What percentage of candidate replay events have cells from both CA1 & CA3 
%and bilaterally?

load '/data13/mcarr/RipplePaper/decodefilterA.mat'
load '/data13/mcarr/RipplePaper/decodefilterB.mat'
cells = cell(length(decodefilterA),1);

for an = 1:length(decodefilterA)
    loadstring = sprintf('%s%scellinfo.mat',decodefilterA(an).animal{2},decodefilterA(an).animal{3});
    load(loadstring)
    cells{an} = cell(length(decodefilterA(an).epochs),1);
    r1 = evaluatefilter(cellinfo,'isequal($hemisphere,''right'') & isequal($area,''CA1'')');
    r3 = evaluatefilter(cellinfo,'isequal($hemisphere,''right'') & isequal($area,''CA3'')');
    l1 = evaluatefilter(cellinfo,'isequal($hemisphere,''left'') & isequal($area,''CA1'')');
    l3 = evaluatefilter(cellinfo,'isequal($hemisphere,''left'') & isequal($area,''CA3'')');
    
    for d = 1:length(decodefilterA(an).epochs)
        cells{an}{d} = cell(size(decodefilterA(an).epochs{d},1),1);
        for e = 1:length(decodefilterA(an).epochs{d})
            cells{an}{d}{e}.all = mean([rowfind(decodefilterA(an).output{d}(e).index,r1)...
                rowfind(decodefilterA(an).output{d}(e).index,r3)...
                rowfind(decodefilterA(an).output{d}(e).index,l1)...
                rowfind(decodefilterA(an).output{d}(e).index,l3)]>0);
            cells{an}{d}{e}.numcells = size(decodefilterA(an).output{d}(e).index,1);

            p = [];
            for event = 1:length(decodefilterA(an).output{d}(e).eventdata)
                ind = find([length(decodefilterA(an).output{d}(e).pvalue)>=event length(decodefilterB(an).output{d}(e).pvalue)>=event]);
                if ~isempty(ind)
                    if length(ind)==2
                        if decodefilterA(an).output{d}(e).pvalue(event)<0.05 || decodefilterB(an).output{d}(e).pvalue(event)<0.05
                            cell_index = decodefilterA(an).output{d}(e).index(unique(decodefilterA(an).output{d}(e).eventdata(event).cellindex),:);
                            p = [p; mean([rowfind(cell_index,r1) rowfind(cell_index,r3) rowfind(cell_index,l1) rowfind(cell_index,l3)]>0)];
                        end
                    elseif ind == 1
                        if decodefilterA(an).output{d}(e).pvalue(event)<0.05
                            cell_index = decodefilterA(an).output{d}(e).index(unique(decodefilterA(an).output{d}(e).eventdata(event).cellindex),:);
                            p = [p; mean([rowfind(cell_index,r1) rowfind(cell_index,r3) rowfind(cell_index,l1) rowfind(cell_index,l3)]>0)];
                        end
                    elseif ind == 2
                        if decodefilterB(an).output{d}(e).pvalue(event)<0.05
                            cell_index = decodefilterA(an).output{d}(e).index(unique(decodefilterA(an).output{d}(e).eventdata(event).cellindex),:);
                            p = [p; mean([rowfind(cell_index,r1) rowfind(cell_index,r3) rowfind(cell_index,l1) rowfind(cell_index,l3)]>0)];
                        end
                    end
                end
            end
            cells{an}{d}{e}.propactive = p;
        end
    end
end

%What % of candidate replay events have CA3 and CA1 cells
active = [];
for an = 1:length(cells)
    for d = 1:length(cells{an})
        for e = 1:length(cells{an}{d})
            if size(cells{an}{d}{e}.propactive,2)==4 && (any(cells{an}{d}{e}.all([1 3])>0) && any(cells{an}{d}{e}.all([2 4])>0))
                active = [active; [sum(cells{an}{d}{e}.propactive(:,[1 3]),2) sum(cells{an}{d}{e}.propactive(:,[2 4]),2)]];
            end
        end
    end
end
both13 = sum(all(active'>0))./size(active,1);

%What % of candidate replay events have left and right hemisphere cells
active = [];
for an = 3%:length(cells)
    for d = 1:length(cells{an})
        for e = 1:length(cells{an}{d})
            if size(cells{an}{d}{e}.propactive,2)==4 && (any(cells{an}{d}{e}.all([1 2])>0) && any(cells{an}{d}{e}.all([3 4])>0))
                active = [active; [sum(cells{an}{d}{e}.propactive(:,[1 2]),2) sum(cells{an}{d}{e}.propactive(:,[3 4]),2)]];
            end
        end
    end
end
bothlr = sum(all(active'>0))./size(active,1);

% 98% (655/667) of significant replay events occur in both CA1 and CA3,
% Bond: 99%, 465/468
% Frank: 94% 143/152
% Ten: 100% 47/47

% 89% (589/661) occur bilaterally.
% Bond: 98% 461/468
% Frank: 59% 86/146
% Ten: 89% 42/47

% candidate events were defined as having at least 5 cells active during 
% the ripple, a p-value < 0.05, and at least one cell recorded in CA1 and
% one cell in CA3 or one in each hemisphere.

% If instead of significant replay events, ask about all candidate events,
% the percentages drop to 95% for CA1 and CA3, 83% bilaterally.

