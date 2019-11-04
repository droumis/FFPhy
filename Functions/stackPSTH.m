

function [PSTH, ETA, ETAShufs] = stackPSTH(signal, eventIdx, varargin)
% make a peri event/stim time hist
% signal center-stacked at eventIdx
% DR 19

runShuffle = 1;
numShufs = 200;
idxWin = -100:100;
if ~isempty(varargin)
    assign(varargin{:})
end

% stack to psth
PSTH = cell2mat(arrayfun(@(r) signal(r+idxWin), eventIdx, 'un', 0));
ETA = nanmean(PSTH);
if isnan(ETA)
    ETA = [];
end

% time shuffle psth
ETAShufs = [];
if runShuffle
for sh=1:numShufs
    PSTHShuf=[];
    for qq=1:length(eventIdx)
        shiftBy=round(rand(1)*numel(idxWin,2));
        PSTHShuf(qq,:)=circshift(PSTH(qq,:),shiftBy,2);
    end
    ETAShufs(sh,:) = nanmean(PSTHShuf);
end
end