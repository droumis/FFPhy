function [intraBoutXP, boutTimes] = getLickBoutLicks(an, eps, varargin)
% [intraBoutXP, boutTimes] = getLickBoutLicks(an, eps, varargin)
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
intraBoutXP = [];
boutTimes = [];
for e = 1:size(eps,1)
   day = eps(e,1);
   ep = eps(e,2);
   eeg = loadeegstruct(andef{2}, an, 'eeg', day, ep, 2);
   epS = eeg{day}{ep}{2}.starttime;
   epE = eeg{day}{ep}{2}.endtime;
   times = [epS:.001:epE]';% epoch ms timevec;
   boutTimes{day}{ep} = vec2list(LburstVec{day}{ep}.lickBout, LburstVec{day}{ep}.time);
   lickTimes = lick{day}{ep}.starttime; 
   lickTimes = lickTimes(isIncluded(lickTimes, boutTimes{day}{ep}));
   lickTimesILI = diff(lickTimes);
   % keep valid lickburst licks
   g = 1;
   while g
       lickTimesILI = diff(lickTimes);
       if any(lickTimesILI < minILIthresh)
            lickTimes(find(lickTimesILI < minILIthresh)+1) = [];   
       else
           g = 0; 
       end
   end
   intraBoutXP = [intraBoutXP; lickTimes];
end
   
