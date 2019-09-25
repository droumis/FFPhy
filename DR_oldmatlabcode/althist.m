%# some random data
X = randn(2500,1);
Y = randn(2500,1)*2;

%# bin centers (integers)
xbins = floor(min(X)):1:ceil(max(X));
ybins = floor(min(Y)):1:ceil(max(Y));
xNumBins = numel(xbins); yNumBins = numel(ybins);

%# map X/Y values to bin indices
Xi = interp1(xbins, 1:xNumBins, X, 'nearest', 'extrap');
Yi = interp1(ybins, 1:yNumBins, Y, 'nearest', 'extrap');

%# limit indices to the range [1,numBins]
Xi = max( min(Xi,xNumBins), 1);
Yi = max( min(Yi,yNumBins), 1);

%# count number of elements in each bin
H = accumarray([Yi(:) Xi(:)], 1, [yNumBins+10 xNumBins+10]);  %%DR .....

%# plot 2D histogram
imagesc(xbins, ybins, H), axis on %# axis image
colormap hot; colorbar
hold on, plot(X, Y, 'b.', 'MarkerSize',1), hold off