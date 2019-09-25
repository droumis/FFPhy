function showTest1d

load(['dynamic_data']);
load(['dynamic1d_result']);
load ../vis/colormap
opts.cmap= cmap;

%vis(data,genmodel,cmap);
visModel(data,result,opts);

