function h = plotspikepos(spikes, pos, ind)
%function h = plotspikepos(spikes, pos, ind)
%  plots positions in grey and spikes as red dots
p = pos{ind(1)}{ind(2)}.data;
s = spikes{ind(1)}{ind(2)}{ind(3)}{ind(4)}.data;

h = plot(p(:,2), p(:,3), '.')
set(h, 'color', [.8 .8 .8]);

hold on
plot(s(:,2), s(:,3), 'r.');

