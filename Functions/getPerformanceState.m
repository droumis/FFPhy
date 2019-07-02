function [corrvec, mistvec] = getPerformanceState(BehaveState, day, epoch, LFPtimes)
%given a time vector (LFPtimes) and BehaveState struct,
%return corresponding vectors for 'correct' and 'mistake' performance state for given day, ep

%get performance state
trialIO = BehaveState.statechanges{day}{epoch}.statechangeseq;
trialIOfields = BehaveState.statechanges{day}{epoch}.fields;
colfields = strsplit(trialIOfields, ' ');

corrcol = find(cellfun(@(x) strcmp(x,'correct'), colfields, 'un', 1));
portcol = find(cellfun(@(x) strcmp(x,'timeportout'), colfields, 'un', 1));
outBcol = find(cellfun(@(x) strcmp(x,'outbound'), colfields, 'un', 1));
inBcol = find(cellfun(@(x) strcmp(x,'inbound'), colfields, 'un', 1));
lastTcol = find(cellfun(@(x) strcmp(x,'lasttime'), colfields, 'un', 1));
currTcol = find(cellfun(@(x) strcmp(x,'currenttime'), colfields, 'un', 1));

outBCorr = trialIO((trialIO(:,outBcol)==1 & trialIO(:,corrcol)==1),[lastTcol currTcol corrcol portcol]);
outBMist = trialIO((trialIO(:,outBcol)==1 & trialIO(:,corrcol)==0),[lastTcol currTcol corrcol portcol]);

corrvec = list2vec(outBCorr(:,[1:2]),LFPtimes);
mistvec = list2vec(outBMist(:,[1:2]),LFPtimes);

end