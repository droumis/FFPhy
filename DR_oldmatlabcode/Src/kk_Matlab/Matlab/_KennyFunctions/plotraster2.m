function plotraster2(timevec,spikes,h,varargin)

%plots raster given a certain format of data :

% timevec is a histc-style vector of times
% spikes is a matrix of 0s and 1s, with each row vector corresponding to a
        % trial
% height is the height of the raster lines
% varargin are plot options, e.g. 'color', 'k', 'linewidth', 4

if isempty(varargin)
    options = {'color', 'k'};
else
    options = varargin;
end

ylevel = 0;
timevec = timevec*1000;  % convert to ms

for t=1:size(spikes,1)
        x = timevec(spikes(t,1:(end-1))==1);
        y = ylevel * ones(size(x));
        line([x ; x],[y ; y + h],options{1:end})
        ylevel = ylevel + h;
end

hold on
% plot zero line
plot(zeros(1,1+size(spikes,1)),0:size(spikes,1),'Color',[.6 .6 1],'LineWidth',1)

% axes
xlim([timevec(1) timevec(end)])
ylim([0 size(spikes,1)])
set(gca,'YDir','reverse')
set(gca,'XTick',mean([timevec(1) timevec(2)]):500:mean([timevec(end-1) timevec(end)]))
% label axes
set(gca,'FontSize',14,'FontWeight','bold')
xlabel('time (ms)','FontSize',16,'FontWeight','bold')
ylabel('trial #','FontSize',16,'FontWeight','bold')

end