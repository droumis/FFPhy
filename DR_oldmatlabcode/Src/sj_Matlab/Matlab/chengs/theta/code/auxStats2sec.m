function [sopts, ana]= auxStats2sec(data,ana)
stept= 2;  % [sec]

auxStatsOpt

mint= data.time(1);
maxt= data.time(end);
for i=1:length(ana.a)
    ana.a{i}.phase= [0;2*pi];
    if isempty(ana.x) | (length(ana.x)<i) | ~isfield(ana.x{i},'time')
        ana.x{i}.time= mint:stept:maxt;
    end
end

