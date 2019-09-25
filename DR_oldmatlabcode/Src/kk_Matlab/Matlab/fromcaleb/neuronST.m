function [h, psth, bins] = neuronST(neur, varargin)

subset = [];
pulselength = 0.010;
psthbinsize = 0.005;
offset = 0;
tickcolor = 'k';
bounds = [-inf, inf];

[otherArgs] = procOptions(varargin);

if isempty(subset)
  spktimes = neur.timestamps;
else
  spktimes = neur.timestamps(subset);
end

N = length(spktimes);

barheight = N+1;

hr = rectangle('position',[0 (-barheight + 1.5) pulselength/10000 barheight]);
set(hr,'facecolor',[0.6 1 1],'linestyle','none');

mint = inf;
maxt = -inf;
for i = 1:N
   h{i} = spikeTrain(spktimes{i}, offset-i+1, 0.8,'linewidth',1,'color',tickcolor);
   if ~isempty(spktimes{i})
      mint = min(mint, min(spktimes{i}));
      maxt = max(maxt, max(spktimes{i}));
   end
end

mint = floor(mint/psthbinsize)*psthbinsize;
maxt = ceil(maxt/psthbinsize)*psthbinsize;
bins = linspace(mint,maxt,(maxt-mint)/psthbinsize + 1);
bins = [bins inf];

h{end+1} = hr;

for i = 1:N
   if ~isempty(spktimes{i})
      psth(:,i) = histc(spktimes{i},bins);
   else
      psth(:,i) = zeros(length(bins),1);
   end
end

psth = mean(psth,2);
