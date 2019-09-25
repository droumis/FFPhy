function h = plotlinfield(lf)
% plots the trajectory data in lf.trajdata
nplots = length(lf.trajdata)
for i = 1:nplots
    subplot(nplots, 1, i)
    plot(lf.trajdata{i}(:,1), lf.trajdata{i}(:,5));
end
if (isfield(lf, 'index'))
    subplot(nplots, 1, 1)
    title(sprintf('index %d %d %d %d', lf.index));
end
