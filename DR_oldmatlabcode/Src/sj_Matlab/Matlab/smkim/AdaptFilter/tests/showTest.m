function showTest

load(['dynamic_data']);

% show AF results
load(['dynamic_data']);
load(['dynamic_result']);
load ../vis/colormap
opts.cmap= cmap;

%vis(data,genmodel,cmap);
visModel(data,result,opts);

