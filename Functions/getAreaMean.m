

function areaDataMean= getAreaMean(ixpc, areainds, ian, iar, datatype)
% eval(sprintf('area%soutputtmp = arrayfun(@(x) ixpc.logbased%sout{ian}{x}(:,:,:), areainds, ''un'', 0);', datatype, datatype));
% eval(sprintf('ixpc.area%smean{ian}{iar} = mean(cat(4,area%soutputtmp{:}), 4);', datatype, datatype));

eval(sprintf('areaoutputtmp = arrayfun(@(x) ixpc.logbased%sout{ian}{x}(:,:,:), areainds, ''un'', 0);', datatype));
areaDataMean = mean(cat(4,areaoutputtmp{:}), 4);

% try
% eval(sprintf('area%szmasktmp = arrayfun(@(y) arrayfun(@(x) ixpc.%szmask{ian}{x}{y}, areainds, ''un'', 0), [6 7], ''un'', 0);', datatype, datatype));
% eval(sprintf('ixpc.area%szmaskmean{ian}{iar}(6:7) = arrayfun(@(x) mean(cat(2,area%szmasktmp{x}{:}), 2), [1 2], ''un'', 0);', datatype, datatype));
% catch
% end
end