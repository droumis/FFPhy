function plotgrouppoints(cellarray, xvalues)

for i = 1:length(cellarray)
    tmpdata = cellarray{i}(:);
    tmpdata(:,2) = xvalues(i);
    plot(tmpdata(:,2),tmpdata(:,1),'.');
    hold on
end
hold off