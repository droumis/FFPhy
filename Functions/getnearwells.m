

function out = getnearwells(animaldir,animalprefix,epochs,varargin)
% out = getnearwells(animaldir,animalprefix,epochs)
% Produces a cell structure with the fields:
% time, nearwell
% EPOCHS - N by 2 matrix, columns are [day epoch]
% intended to be used to seperate periods when the animal is near any well
% or away from the well

% Author Demetris Roumis June 2019

plotfigs = 0;
cmfromwell = 25;
if ~isempty(varargin)
    assign(varargin{:});
end

loaddays = unique(epochs(:,1));
linpos = loaddatastruct(animaldir, animalprefix, 'linpos', loaddays);

for i = 1:size(epochs,1)
    day = epochs(i,1);
    epoch = epochs(i,2);
    
    % mark times near any well
    lindist = linpos{day}{epoch}.statematrix.lindist;
    nearcwell = lindist < cmfromwell;
    nearouterwell = lindist > (max(lindist)-cmfromwell);
    nearwell = any([nearcwell nearouterwell],2);

    out{day}{epoch}.nearwell = nearwell;
    out{day}{epoch}.nearwell_list = nearwell;
    out{day}{epoch}.time = linpos{day}{epoch}.statematrix.time;
    
    if plotfigs
        figure(1)
        subplot(1,2,1)
        cla
        histogram(lindist)
        hold on;
        yl = ylim;
        line([cmfromwell cmfromwell], [yl(1) yl(2)], 'Color', 'black')
        line([(max(lindist)-cmfromwell) (max(lindist)-cmfromwell)], [yl(1) yl(2)], 'Color', 'black')
        title(sprintf('%s %d %d lindist', animalprefix, day, epoch))
        hold off;
        
        pos = loaddatastruct(animaldir, animalprefix, 'pos', day);
        xp = pos{day}{epoch}.data(:,6);
        yp = pos{day}{epoch}.data(:,7);
        subplot(1,2,2)
        a = scatter(xp,yp,4,nearwell, 'filled', 'MarkerFaceAlpha', .1);
        colormap(jet)
        axis tight
        title('near well 2D')
        pause
    end
end
