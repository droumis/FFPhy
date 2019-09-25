function out = plotvelocitydist(index, excludetimes, pos, varargin)

p = pos{index(1)}{index(2)};
clear pos

xbin = [floor(min(p.data(:,2))):2:ceil(max(p.data(:,2)))];
ybin = [floor(min(p.data(:,3))):2:ceil(max(p.data(:,3)))];

subsx = lookup(p.data(:,2),xbin);
subsy = lookup(p.data(:,3),ybin);

A = accumarray([subsx subsy],p.data(:,5),[],@(x) mean(x),NaN);

figure
imagesc(A)
colorbar

colormap('default');

cmap = jet(1024);
cmap = cmap(2:920,:);
cmap(1,:) = 1;
colormap(cmap);
set(gca,'clim',[0 32])

out = A;

end