function [corrvec, mistvec] = getPerformanceState(BehaveState, day, epoch, LFPtimes)
%given a time vector (LFPtimes) and BehaveState struct,
%return corresponding vectors for 'correct' and 'mistake' performance state for given day, ep

%get performance state
trialIO = BehaveState.statechanges{day}{epoch}.statechangeseq;
trialIOfields = BehaveState.statechanges{day}{epoch}.fields;
corrcol = find(cell2mat(cellfun(@(x) strcmp(x,'correct'), strsplit(trialIOfields, ' '), 'UniformOutput', false)));
timeportoutcol = find(cell2mat(cellfun(@(x) strcmp(x,'timeportout'), strsplit(trialIOfields, ' '), 'UniformOutput', false)));
outBcol = find(cell2mat(cellfun(@(x) strcmp(x,'outbound'), strsplit(trialIOfields, ' '), 'UniformOutput', false)));
inBcol = find(cell2mat(cellfun(@(x) strcmp(x,'inbound'), strsplit(trialIOfields, ' '), 'UniformOutput', false)));
lastTimecol = find(cell2mat(cellfun(@(x) strcmp(x,'lasttime'), strsplit(trialIOfields, ' '), 'UniformOutput', false)));
currTimecol = find(cell2mat(cellfun(@(x) strcmp(x,'currenttime'), strsplit(trialIOfields, ' '), 'UniformOutput', false)));
outBCorr = trialIO((trialIO(:,outBcol)==1 & trialIO(:,corrcol)==1),[lastTimecol currTimecol corrcol timeportoutcol]);
outBMist = trialIO((trialIO(:,outBcol)==1 & trialIO(:,corrcol)==0),[lastTimecol currTimecol corrcol timeportoutcol]);
corrvec = list2vec(outBCorr(:,[1:2]),LFPtimes);
mistvec = list2vec(outBMist(:,[1:2]),LFPtimes);

end