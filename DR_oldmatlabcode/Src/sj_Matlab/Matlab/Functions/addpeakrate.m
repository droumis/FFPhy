function addpeakrate(animdirect,fileprefix,index)

includeStates = 6;

load([animdirect, fileprefix,'cellinfo']);
for i = 1:size(index,1)
    if index(i,1) < 10
        dsz = '0';
    else
        dsz = '';
    end
    %Load spike and linpos structures
    load([animdirect, fileprefix,'spikes',dsz,num2str(index(i,1)),'.mat'])
    load([animdirect, fileprefix,'linpos',dsz,num2str(index(i,1)),'.mat'])
    
    %Run getbehave state
    [state, lindist] = getbehavestate(linpos, index(i,1), index(i,2), includeStates);
    
    %Determine the index of all cells recorded for this session
    cellindex = evaluatefilter(cellinfo{index(i,1)}{index(i,2)},'$numspikes > 1');
    for c = 1:size(cellindex,1)
        trajdata = calclinfields(spikes,state,lindist,linpos, [index(i,:) cellindex(c,:)]);
        peakrate = zeros(length(trajdata),1);
        for t = 1:length(trajdata)
            peakrate(t) = max(trajdata{t}(:,5));
        end
        cellinfo{index(i,1)}{index(i,2)}{cellindex(c,1)}{cellindex(c,2)}.peakrate = max(peakrate);
    end
end

save([animdirect, fileprefix,'cellinfo'], 'cellinfo'); 