function out = calcbinnedtrend(xdata,ydata,binsize)

% out = calcbinnedtrend(xdata,ydata,binsize)
% Bins the data in ydata across binned xdata (using binsize).  Also
% calculates the stderr of each bin if more than 10 data points exist.

xdata = xdata(:);
ydata = ydata(:);

xbins = [min(xdata):binsize:max(xdata)]';
binneddata = {[]};
errorvector = [];
meanvector = [];
medianvector = [];

tmpresult = [lookup(xdata,xbins) ydata];
for j = 1:size(tmpresult)
    tmpbin = tmpresult(j,1);
    if (tmpbin > 0)
        if length(binneddata) < tmpbin
            binneddata{tmpbin} = [];
        end
        if ~isnan(tmpresult(j,2))
            binneddata{tmpbin} = [binneddata{tmpbin}; tmpresult(j,2)];
        end
    end
end

for i = 1:length(binneddata)
    if isempty(binneddata{i})
        xbins(i) = nan;
        errorvector(i,1:2) = nan;
        meanvector(i,1) = nan;
        medianvector(i,1) = nan;
    else
        meanvector(i,1) = mean(binneddata{i});
        medianvector(i,1) = median(binneddata{i});
        if (length(binneddata{i}) > 5)
            tmperror = repmat(std(bootstrp(1000,@mean,binneddata{i})), [1 2]);
            errorvector(i,1:2) = tmperror;
            
            %tmperror = prctile(bootstrp(1000,@mean,binneddata{i}),[5 95]);
            %errorvector(i,1:2) = abs(meanvector(i) - tmperror);
            
        else
            errorvector(i,1:2) = nan;
        end
        
    end
end

valid = find(~isnan(xbins));
out.xbins = xbins(valid);
out.mean = meanvector(valid);
out.median = medianvector(valid);
out.error = errorvector(valid,:);
out.binneddata = binneddata;


