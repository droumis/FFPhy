function plotraster4(timevec,outcomes,h,varargin)

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

% plot errors as long grey horizontal bars
for t=1:size(outcomes,1)
        
        y = ylevel;
        
        if outcomes(t)==0    % if error
            xvert = [timevec(1) ;  timevec(end) ; timevec(end) ; timevec(1) ];
            yvert = [ylevel     ;  ylevel       ; ylevel + h ; ylevel + h];  
            patch(xvert,yvert,[1 .9 .9],'EdgeColor','none');
        end
        
        % move to next line for next trial
        ylevel = ylevel + h;
        
end

hold on
% plot zero line
plot(zeros(1,1+size(outcomes,1)),0:size(outcomes,1),'Color',[.8 .8 .8],'LineWidth',1)

% axes
xlim([timevec(1) timevec(end)])
ylim([0 size(outcomes,1)])
set(gca,'YDir','reverse')
set(gca,'XTick',mean([timevec(1) timevec(2)]):500:mean([timevec(end-1) timevec(end)]))
% label axes
set(gca,'FontSize',14,'FontWeight','bold')
xlabel('time (ms)','FontSize',16,'FontWeight','bold')
ylabel('trial #','FontSize',16,'FontWeight','bold')

end