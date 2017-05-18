function out = getWtracktrialstate(FFanimdir,animalID, totalepochs, varargin)

% Demetris Roumis May 2017.. get the trial state times for wtrack..
% returns vecs inbound, outbound, correct, time

if (~isempty(varargin))
    assign(varargin{:});
end

% load behavstate
load(sprintf('%s%s%s.mat',FFanimdir, animalID, 'BehaveState'));

% get the state space performance score time vec
SSallboundPerformance = getColfromMat(BehaveState.statespace, 'allbound', {'mode', 'day', 'epoch'}, 'fieldstring','pcFields');
SSallboundLastTime = getColfromMat(BehaveState.statespace, 'allepsMat', {'lasttime'}, 'fieldstring','allepsMatFields');
SSallboundCurrentTime = getColfromMat(BehaveState.statespace, 'allepsMat', {'currenttime'}, 'fieldstring','allepsMatFields');
SSinboundInd = getColfromMat(BehaveState.statespace, 'allepsMat', {'inbound'}, 'fieldstring','allepsMatFields');
SSoutboundInd = getColfromMat(BehaveState.statespace, 'allepsMat', {'outbound'}, 'fieldstring','allepsMatFields');
SSinboundPerformance = getColfromMat(BehaveState.statespace, 'inbound', {'mode', 'day', 'epoch'}, 'fieldstring','pcFields');
SSoutboundPerformance = getColfromMat(BehaveState.statespace, 'outbound', {'mode', 'day', 'epoch'}, 'fieldstring','pcFields');

SSallboundtimes = [SSallboundLastTime SSallboundCurrentTime];
% SSinboundtimes = SSallboundtimes(find(SSinboundInd),:);
% SSoutboundtimes = SSallboundtimes(find(SSoutboundInd),:);


% days = unique(totalepochs(:,1));
for dayep = 1:length(totalepochs(:,1));
    day = totalepochs(dayep,1);
    epoch = totalepochs(dayep,2);
    
    statechanges = BehaveState.statechanges{day}{epoch};
    
    %find column of interest and assign in
    [currenttime currtimecol] = getColfromMat(statechanges,'statechangeseq', {'currenttime'});
    [~, lasttimecol] = getColfromMat(statechanges, 'statechangeseq',{'lasttime'});
    correct = getColfromMat(statechanges, 'statechangeseq', {'correct'});
    inbound = getColfromMat(statechanges, 'statechangeseq',{'inbound'});
    outbound = getColfromMat(statechanges, 'statechangeseq',{'outbound'});
    
    %create time vec spanning the valid trial times for this epoch
    trialstimevec = [currenttime(1):1/1000:currenttime(end)]'; %ms accuracy.. not using 30kHz clockrate because unneccesary precision and huge vec

    correcttimes = statechanges.statechangeseq(correct==1,[lasttimecol currtimecol]); %get the correct trial start - end times
    correcttimeslogic = list2vec(correcttimes, trialstimevec);

    outboundtimes = statechanges.statechangeseq(outbound==1,[lasttimecol currtimecol]); %get the outbound trial start - end times
    outboundtimeslogic = list2vec(outboundtimes, trialstimevec);

    inboundtimes = statechanges.statechangeseq(inbound==1,[lasttimecol currtimecol]); %get the inbound trial start - end times
    inboundtimeslogic = list2vec(inboundtimes, trialstimevec);

    SSall = [SSallboundtimes SSallboundPerformance SSinboundInd SSoutboundInd];
    SSall = SSall(SSall(:,4) == day & SSall(:,5) == epoch,:);
    dayepSSallBtimes = SSall(:, [1 2]);
    dayepSSinBtimes = SSall(SSall(:,6) == 1, [1 2]);
    dayepSSoutBtimes = SSall(SSall(:,7) == 1, [1 2]);

    SSinboundtimesScore = list2vec(dayepSSinBtimes, trialstimevec, 'marker', SSinboundPerformance);
    SSoutboundtimesScore = list2vec(dayepSSoutBtimes, trialstimevec, 'marker', SSoutboundPerformance);
    SSallboundtimesScore = list2vec(dayepSSallBtimes, trialstimevec, 'marker', SSallboundPerformance);
    
    out{day}{epoch}.time = trialstimevec;
    out{day}{epoch}.inbound = inboundtimeslogic;
    out{day}{epoch}.outbound = outboundtimeslogic;
    out{day}{epoch}.correct = correcttimeslogic;
    
    out{day}{epoch}.inboundScore = SSinboundtimesScore;
    out{day}{epoch}.outboundScore = SSoutboundtimesScore;
    out{day}{epoch}.allboundScore = SSallboundtimesScore;
    

end
end

function [allcoldata, cols] = getColfromMat(statestruct, structarray, colstringcell, varargin)
if ~iscell(colstringcell)
    error('column string cell array must be a cell array of strings');
end
    fieldstring = 'fields';
if (~isempty(varargin))
    assign(varargin{:});
end
cols = [];
allcoldata = [];
for isc = 1:length(colstringcell)
    icol = eval(['find(cell2mat(cellfun(@(x) strcmp(x,colstringcell{isc}),strsplit(statestruct.', fieldstring, ', '' ''), ''UniformOutput'', false)))']);
    coldata = eval(['statestruct.' structarray '(:,icol)']);
    cols = [cols icol];
    allcoldata = [allcoldata coldata];
end
end