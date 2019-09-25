% mutual information of the generative model vs. bin size
% make figure: mutual information of the generative model vs. bin size
%  1. from histogramm-ed spike data  2. from exact distribution

load mInfoDist 

figure(1); clf; 
subplot(2,1,1);
imagesc(nbinx,nbiny,I);
ylabel('number of bins in phase')
axis xy
colorbar
title('mutual info from generated spikes');

subplot(2,1,2);
imagesc(nbinx,nbiny,In);
xlabel('number of bins in position')
ylabel('number of bins in phase')
title('mutual info from exact distribution');
axis xy
colorbar
figopts=struct('width',3,'height',5,...
  'FontMode','fixed','FontSize',8,'LockAxes',0,'Color','rgb',...
  'BoundsCode','mcode');
exportfig(gcf,'mInfoDist',figopts);

figure(3); clf; 
subplot(3,1,1);
imagesc(nbinx,nbiny,corr)
title('lin. correll. from generated spikes');
colorbar

subplot(3,1,2);
imagesc(nbinx,nbiny,corrn)
title('lin. correll. from exact distribution');
colorbar

subplot(3,1,3);
imagesc(nbinx,nbiny,corr2)
title('lin. correll. directly from generated spikes');
colorbar

