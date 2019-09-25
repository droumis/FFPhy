% MATLAB startup file for smkim

% undocumented feature, may fix memory leak?
feature jitallow fast_jit_call_off;

format long g
set(0,'DefaultFigureRenderer','painters');
set(0,'DefaultLineLineWidth',1);
set(0,'DefaultLineColor','k');
set(0,'DefaultLineMarkerSize',3);
set(0,'DefaultAxesYDir','normal');
set(0,'DefaultPatchFaceColor',[0.5 0.5 0.5]);
set(0,'DefaultPatchEdgeColor','k');
set(0,'DefaultSurfaceFaceColor','interp');
set(0,'DefaultSurfaceEdgeColor','none');
set(0,'DefaultAxesTickDir','out');
set(0,'DefaultAxesColorOrder',[0 0 0]);
set(0,'DefaultAxesTickLength',[0.01 0.01]);
set(0,'DefaultAxesLayer','top');
set(0,'DefaultAxesColor','none');
set(0,'DefaultAxesYDir','normal');
set(0,'DefaultAxesXDir','normal');
set(0,'DefaultAxesLineWidth',2);
set(0,'DefaultAxesFontSize',14) 
set(0,'DefaultAxesFontName','Arial') 

addpath('/home/smkim/scratch');
addpath('/home/smkim/matlab');
addpath('/home/smkim/matlab/batch');
addpath('/home/smkim/matlab/extrema');
addpath('/data14/smkim/');
addpath('/data14/smkim/spectral_filters');
addpath(genpath('/home/smkim/matlab/chronux/spectral_analysis'));
addpath(genpath('/home/smkim/matclust'));


