function h = plotSpikeMatrix(tVect,spkMat,varargin)

offset = 0;
tickcolor = [0 0 1];

[otherArgs] = procOptions(varargin);

M = size(spkMat,1);

holdstate = ishold;
hold on

h = [];

o = 0;
for m = 1:M
  if isnan(spkMat(m,1))
    continue;
  end
  o = o + 1;
  spktimes = tVect(find(spkMat(m,:)));
  if isempty(spktimes)
    continue;
  end
  %h{m} = spikeTrain(spktimes, offset-o+1, 0.8,'linewidth',1,'color',tickcolor,'marker','none');
  h(m) = plot(spktimes, (offset-o+1)*ones(length(spktimes),1),'.','color',tickcolor);
end

if ~holdstate
  hold off;
end
