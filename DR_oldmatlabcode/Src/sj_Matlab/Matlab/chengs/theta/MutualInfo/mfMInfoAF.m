load mInfoAF

figure(1); clf
imagesc(nbinx, nbiny, minfo');
xlabel('number of bins in position')
ylabel('number of bins in phase')
axis xy
colorbar
title('mutual info from adaptive filter estimate')
figopts=struct('width',3,'height',2.5,...
  'FontMode','fixed','FontSize',8,'LockAxes',0,'Color','rgb',...
  'BoundsCode','mcode');
exportfig(gcf,'mInfoAF',figopts);
