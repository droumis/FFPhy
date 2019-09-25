function [output] = plotsj_2dmaps(twodoccupancyrunepoch1, twodoccupancyrunepoch2, twodoccupancyrunepoch3, index)

%Plots the 2d occupancy normalized firing rate - This is organized in separate trajectories
%
%twodoccupancyrunepochn - the output from twodoccupancy for run epoch n for the day, tetrode,
%and cell you are analyzing
% twodoccupancy separates by trajectory, openfieldrate does not
%index - [day tetrode cell]
%

colormap('default');

cmap = jet(1024) ./ 1.5;
    cmap = cmap(100:920,:);
    cmap(1,:) = 1;
    colormap(cmap);

%now we need to make the subplots for each trajectory and epoch, and make
%sure that the x/y axis is the same for all subplots

    max1 = max(twodoccupancyrunepoch1(1).smoothedspikerate(:));
    max2 = max(twodoccupancyrunepoch2(1).smoothedspikerate(:));
    max3 = max(twodoccupancyrunepoch3(1).smoothedspikerate(:));
    max4 = max(twodoccupancyrunepoch1(2).smoothedspikerate(:));
    max5 = max(twodoccupancyrunepoch2(2).smoothedspikerate(:));
    max6 = max(twodoccupancyrunepoch3(2).smoothedspikerate(:));
    valuemax = max([max1 max2 max3 max4 max5 max6]);
    bounds = [-1 valuemax*.65];
    [x1] = min(twodoccupancyrunepoch1(1).xticks);
    [x2] = min(twodoccupancyrunepoch1(2).xticks);
    [x3] = min(twodoccupancyrunepoch2(1).xticks);
    [x4] = min(twodoccupancyrunepoch2(2).xticks);
    [x5] = min(twodoccupancyrunepoch3(1).xticks);
    [x6] = min(twodoccupancyrunepoch3(2).xticks);
    [x1m] = max(twodoccupancyrunepoch1(1).xticks);
    [x2m] = max(twodoccupancyrunepoch1(2).xticks);
    [x3m] = max(twodoccupancyrunepoch2(1).xticks);
    [x4m] = max(twodoccupancyrunepoch2(2).xticks);
    [x5m] = max(twodoccupancyrunepoch3(1).xticks);
    [x6m] = max(twodoccupancyrunepoch3(2).xticks);

    [y1] = min(twodoccupancyrunepoch1(1).yticks);
    [y2] = min(twodoccupancyrunepoch1(2).yticks);
    [y3] = min(twodoccupancyrunepoch2(1).yticks);
    [y4] = min(twodoccupancyrunepoch2(2).yticks);
    [y5] = min(twodoccupancyrunepoch3(1).yticks);
    [y6] = min(twodoccupancyrunepoch3(2).yticks);
    [y1m] = max(twodoccupancyrunepoch1(1).yticks);
    [y2m] = max(twodoccupancyrunepoch1(2).yticks);
    [y3m] = max(twodoccupancyrunepoch2(1).yticks);
    [y4m] = max(twodoccupancyrunepoch2(2).yticks);
    [y5m] = max(twodoccupancyrunepoch3(1).yticks);
    [y6m] = max(twodoccupancyrunepoch3(2).yticks);

    xmin = floor(min([x1 x2 x3 x4 x5 x6]));
    xmax = ceil(max([x1m x2m x3m x4m x5m x6m]));
    ymin = floor(min([y1 y2 y3 y4 y5 y6]));
    ymax = ceil(max([y1m y2m y3m y4m y5m y6m]));
    
subplot(3, 2, 1)    
h = imagesc(twodoccupancyrunepoch1(1).xticks, twodoccupancyrunepoch1(1).yticks, twodoccupancyrunepoch1(1).smoothedspikerate, bounds);
axis([xmin xmax ymin ymax]);

subplot(3, 2, 2)    
h = imagesc(twodoccupancyrunepoch1(2).xticks, twodoccupancyrunepoch1(2).yticks, twodoccupancyrunepoch1(2).smoothedspikerate, bounds);
axis([xmin xmax ymin ymax]);

subplot(3, 2, 3)    
h = imagesc(twodoccupancyrunepoch2(1).xticks, twodoccupancyrunepoch2(1).yticks, twodoccupancyrunepoch2(1).smoothedspikerate, bounds);
axis([xmin xmax ymin ymax]);

subplot(3, 2, 4)    
h = imagesc(twodoccupancyrunepoch2(2).xticks, twodoccupancyrunepoch2(2).yticks, twodoccupancyrunepoch2(2).smoothedspikerate, bounds);
axis([xmin xmax ymin ymax]);

subplot(3, 2, 5)    
h = imagesc(twodoccupancyrunepoch3(1).xticks, twodoccupancyrunepoch3(1).yticks, twodoccupancyrunepoch3(1).smoothedspikerate, bounds);
axis([xmin xmax ymin ymax]);

subplot(3, 2, 6)    
h = imagesc(twodoccupancyrunepoch3(2).xticks, twodoccupancyrunepoch3(2).yticks, twodoccupancyrunepoch3(2).smoothedspikerate, bounds);
axis([xmin xmax ymin ymax]);
colorbar('SouthOutside');