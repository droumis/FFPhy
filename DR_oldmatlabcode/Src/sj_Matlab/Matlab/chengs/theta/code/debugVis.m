function debugVis(rat, din, ein, tin, cin, ratio)
% ratio: movietime/ realtime  [only works for 2ms timesteps]
% 

d=str2num(din); e= str2num(ein); t= str2num(tin); c= str2num(cin);

load(['/home/chengs/theta/' rat '/data2/behavdata' sprintf('%.2d', d)]);
load(['/home/chengs/theta/' rat '/data2/spikedata' sprintf('%.2d', d)]);


data= behavdata{d}{e};
data.spiketimes= spikedata{d}{e}{t}{c}.time;
data.spikeindex= spikedata{d}{e}{t}{c}.index;


load(['adaptest' sprintf('%.2d', d)]);
model= adaptest{d}{e}{t}{c}.model;

global adaptest cmap
%load /home/chengs/AdaptFilter/newvis/colormap;
load /home/chengs/AdaptFilter/newvis/colormap2;
opt.playrate= str2num(ratio)*500;
opt.cmap= cmap;
%opt.cmap= jet(1024);

visModel(data, model, opt)
