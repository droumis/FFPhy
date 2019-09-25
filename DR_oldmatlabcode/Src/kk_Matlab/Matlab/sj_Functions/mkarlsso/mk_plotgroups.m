function H = plotgroups(cellarray, plotopt,column)

gstats = calcgroupstats(cellarray);
if (nargin < 2)
    plotopt = '';
end
if (nargin < 3)
    column = 1;
end
gstats = calcgroupstats(cellarray,column);
for i = 1:length(gstats)
    if (length(gstats(i).mean) > 1)
        gstats(i).stderr = stderr(gstats(i).mean);
        gstats(i).mean = mean(gstats(i).mean);
    end
    if isempty(gstats(i).mean)
        gstats(i).groupnum = [];
    end
end
errorbar([gstats(:).groupnum]',[gstats(:).mean]', [gstats(:).bootstderr]', plotopt);




