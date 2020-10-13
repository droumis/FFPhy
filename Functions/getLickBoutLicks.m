function [intraBoutXP, boutTimes] = getLickBoutLicks(an, eps, varargin)
% [intraBoutXP, boutTimes] = getLickBoutLicks(an, eps, varargin)
% get lick bout XP events.
% calls getLickBout

% args:
%
% varargs:
%
% output:
%
%{
Contumacious Clock
FFPhy
@DR
%}
maxILIthresh = 1; % max burst ili threshold in seconds
minILIthresh = .06; % min burst ili threshold in seconds
minBoutLicks = 2; %filter out bouts with less than boutNum licks
lick = [];
if ~isempty(varargin)
    assign(varargin{:});
end
andef = animaldef(an);
if isempty(lick)
    andef = animaldef(an);
    loaddays = unique(eps(:,1));
    lick = loaddatastruct(andef{2}, an, 'lick', loaddays);
end

intraBoutXP = [];
boutTimes = [];
for e = 1:size(eps,1)
    day = eps(e,1);
    ep = eps(e,2);
    if isa(lick, 'cell')
        % get lick bouts
        try
            lickTime = lick{day}{ep}.eventtime;
        catch
            lickTime = lick{day}{ep}.starttime; % legact name
        end
    else
        lickTime = lick;
    end

%    eeg = loadeegstruct(andef{2}, an, 'eeg', day, ep, 2);
%    epS = eeg{day}{ep}{2}.starttime;
%    epE = eeg{day}{ep}{2}.endtime;
%    times = [epS:.001:epE]';% epoch ms timevec;
   LburstOut = getLickBout([], an, eps, 'lick', lickTime, 'maxILIthresh', maxILIthresh, ...
       'minILIthresh', minILIthresh, 'minBoutLicks', minBoutLicks, 'output_intervals', 1);
   boutTimes{day}{ep} = LburstOut{day}{ep}.boutTimes;
   intraBoutXP{day}{ep} = LburstOut{day}{ep}.lickBoutXP;
%    boutTimes{day}{ep} = vec2list(LburstVec{day}{ep}.lickBout, LburstVec{day}{ep}.time);   
%    XpTime = vec2list(LburstVec{day}{ep}.lickBoutXP, LburstVec{day}{ep}.time);
%    lickTime = lickTime(isIncluded(lickTime, boutTimes{day}{ep}));
end
