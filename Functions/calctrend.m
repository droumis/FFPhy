function [out, confbounds] = calctrend(xdata,ydata,binsize,smoothwindow)
%[out, confbounds] = calctrend(xdata,ydata,binsize,smoothwindow)
%Calculates the smoothed relationship between x and y and computed 95%
%confidence bounds on the mean curve using a bootstrap measure. BINSIZE is
%the final bin step of the smoothed data.  SMOOTHWINDOW is the stddev of
%the gaussian used to smooth neighboring points (number of nearest data points to smooth
%over).

data = sortrows([xdata(:), ydata(:)]);
xdatamin = min(data(:,1));
xdatamax = max(data(:,1));
newxdata = [xdatamin:binsize:xdatamax]';

bootvalues = bootstrp(1000, @calcmeantrend, data, xdatamin, xdatamax, binsize, smoothwindow);

newydata = mean(bootvalues)';
confbounds = prctile(bootvalues,[5 95])';

out = [newxdata newydata];

function out = calcmeantrend(data, xdatamin, xdatamax, binsize, smoothwindow)

newxdata = [xdatamin:binsize:xdatamax]';
data = sortrows(data);

%smooth the data before the interpolation
data(:,2) = smoothvect(data(:,2),gaussian(smoothwindow,smoothwindow*4));

%get rid of any repeats if they exist 
[trash, ind] = unique(data(:,1));
data = data(ind,:);

newydata = interp1(data(:,1),data(:,2),newxdata);

%fill in any nan values
nanind = find(isnan(newydata));
newydata(nanind) = -1111111.11;
newydata = vectorfill(newydata,-1111111.11);
out = newydata;