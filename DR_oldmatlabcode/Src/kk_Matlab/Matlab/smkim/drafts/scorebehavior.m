

% score behavior for recording animal
load('trajectories.mat');
load('smoothedpos.mat');

for d = 1:length(trajectories)
    for e = 1:3
        j = 2*e; % session index
        ncorrect(d,e) = nnz(trajectories{d}{j}.correct);
        ntotal(d,e) = nnz(~isnan(trajectories{d}{j}.correct));
        sinuosity(d,e) = pathsinuosity(smoothedpos{d}{j},10);
        speed(d,e) = mean(hypot(smoothedpos{d}{j}.data(:,4),smoothedpos{d}{j}.data(:,5)));
    end
end

graycolor = [0.7 0.7 0.7];

subject = 'S51';
for d = 5
    for e = 1:3
        j = 2*e;
        figure(e)
        line(smoothedpos{d}{j}.data(:,2),smoothedpos{d}{j}.data(:,3),'Color',graycolor);
        runtext{d,e}{1} = [num2str(ncorrect(d,e)) '/' num2str(ntotal(d,e)) ' correct (' sprintf('%0.2f',ncorrect(d,e)/ntotal(d,e)) ')'];
        runtext{d,e}{2} = ['mean speed ' sprintf('%02.1f',speed(d,e)) ' cm/s'];
        runtext{d,e}{3} = ['path sinuosity ' sprintf('%0.2f',sinuosity(d,e))];
        caption{d,e} = ['\bfrun' num2str(e)];
        if e == 2
            text(-50,-10,runtext{d,e},'FontSize',16,'Color','r', ...
              'VerticalAlignment','top','HorizontalAlignment','left');
            text(0,170,caption{d,e},'FontSize',20,'Color','r', ...
              'VerticalAlignment','bottom','HorizontalAlignment','center');
        else
            text(-50,-10,runtext{d,e},'FontSize',16,'Color','b', ...
              'VerticalAlignment','top','HorizontalAlignment','left');
            text(0,170,caption{d,e},'FontSize',20,'Color','b', ...
              'VerticalAlignment','bottom','HorizontalAlignment','center');
        end
            
        set(gca,'DataAspectRatio',[1 1 1],'Visible','off','YLim',[-30 180],'XLim',[-70 70]);

        figure(e); print('-dpng',[subject '_day' num2str(d) '_run' num2str(e) '_behavior.png']);
    end
end


