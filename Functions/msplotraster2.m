function msplotraster2(timevec,spikes,plotrewardlines,h,varargin)

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

hold on

if plotrewardlines
% plot zero line
plot(zeros(1,1+size(spikes,1)),0:size(spikes,1),'-','Color',[0.5 0.5 0.5],'LineWidth',1.5)
%plot reward line
plot(ones(1,1+size(spikes,1)).*2,0:size(spikes,1),'-','Color',[1 0.6 0],'LineWidth',1.5)
end

for t=1:size(spikes,1)
        x = timevec(spikes(t,1:(end-1))==1);
        y = ylevel * ones(size(x));
        line([x ; x],[y ; y + h],options{1:end})
        ylevel = ylevel + h;
end



% axes
xlim([timevec(1) round(timevec(end))])
ylim([0 size(spikes,1)])
set(gca,'YDir','reverse')

% label axes
set(gca,'FontSize',14,'FontWeight','bold')
xlabel('time (s)','FontSize',16,'FontWeight','bold')
ylabel('trial #','FontSize',16,'FontWeight','bold')

end