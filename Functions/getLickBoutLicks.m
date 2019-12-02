function [LickBoutLicks] = getLickBoutLicks(an, eps, varargin)
% [out] = getWetLickBouts(andir,an,eps,varargin)
%
% args:
%
% varargs:
%
% output:
%
%{
FFPhy V0.1
@DR
%}
minILIthresh = .06; % seconds
if ~isempty(varargin)
    assign(varargin{:});
end
andef = animaldef(an);
loaddays = unique(eps(:,1));
lick = loaddatastruct(andef{2}, an, 'lick', loaddays);
LburstVec = getLickBout([], an, eps, varargin);
LickBoutLicks = [];
for e = 1:size(eps,1)
   day = eps(e,1);
   ep = eps(e,2);
   eeg = loadeegstruct(andef{2}, an, 'eeg', day, ep, 2);
   epS = eeg{day}{ep}{2}.starttime;
   epE = eeg{day}{ep}{2}.endtime;
   times = [epS:.001:epE]';% epoch ms timevec;
   LBIntervals = vec2list(LburstVec{day}{ep}.lickBout, LburstVec{day}{ep}.time);
   lickTimes = lick{day}{ep}.starttime; 
   lickTimes = lickTimes(logical(isExcluded(lickTimes, LBIntervals))); %isIncluded
   lickTimesILI = diff(lickTimes);
   % keep valid lickburst licks
   f = lickTimesILI > minILIthresh;
   LickBoutLicks = [LickBoutLicks; lickTimes(find([f(1); f]))];
end
   
