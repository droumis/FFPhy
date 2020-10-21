function [out] = getWetLickBouts(andir,an,eps,varargin)
% [out] = getWetLickBouts(andir,an,eps,varargin)
% -timefilter function to be called from setfilterfunction
% -returns vec of time, wetLB, dryLB
%
% args:
%
% varargs:
%
% output:
%
%
%{
FFPhy V0.1
@DR
%}

maxTimeSinceRew = 2; % seconds max since rew to classify as wet
if ~isempty(varargin)
    assign(varargin{:});
end

andef = animaldef(an);
loaddays = unique(eps(:,1));
task = loaddatastruct(andef{2}, an, 'task', loaddays);
DIO = loaddatastruct(andef{2}, an, 'DIO', loaddays);

LburstVec = getLickBout(andir, an, eps, varargin{:});

for e = 1:size(eps,1)
   day = eps(e,1);
   ep = eps(e,2);
   eeg = loadeegstruct(andef{2}, an, 'eeg', day, ep, 2);
   epS = eeg{day}{ep}{2}.starttime;
   epE = eeg{day}{ep}{2}.endtime;
   times = [epS:.001:epE]';% epoch ms timevec;
   % get reward DIO ID from taskInfo to index into DIO
   rewDIOID = task{day}{ep}.outputdio;
   isOutput = cellfun(@(x) isequal(x.input,0), DIO{day}{ep}, 'un', 1);
   dioID = cellfun(@(x) str2double(regexp(x.original_id,'\d*','Match')),...
       DIO{day}{ep}, 'un', 1);
   rewDIOIdx = find(all([ismember(dioID,rewDIOID)' isOutput'],2));
   rewPortDOut = DIO{day}{ep}(rewDIOIdx);
   
   % for each lick port Dout, gather it's times along with index
   l = cellfun(@(x) x.times, rewPortDOut,'un',0)';
   [rewTime,X] = sort(cell2mat(l));
   g = [];
   for d = 1:length(rewDIOID)
       g{d,1} = repmat(rewDIOID(d), length(l{d}),1);
   end
   i = cell2mat(g);
   id = i(X);
   LburstSE = vec2list(LburstVec{day}{ep}.lickBout, LburstVec{day}{ep}.time);
   [LBrewIdx, LBtimeSinceRew] = knnsearch(rewTime, LburstSE(:,1));   
   isWetLB = LBtimeSinceRew < maxTimeSinceRew;
    
   % convert from list form back to vector form
   isWetTimeVec = list2vec(LburstSE(isWetLB,:),times)';
   isDryTimeVec = list2vec(LburstSE(~isWetLB,:),times)';
   
%    out{day}{ep}.timeSinceRew = timeSinceRew;
   out{day}{ep}.time = times;
%    out{day}{ep}.burstID = [];
%    out{day}{ep}.wellID = [];
   out{day}{ep}.wetLB = isWetTimeVec';
   out{day}{ep}.dryLB = isDryTimeVec';
end
end