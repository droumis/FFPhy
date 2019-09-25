% KERNAREATEST: Test kernarea as function of h (smoothing parameter)
%
%     Usage: kernareatest(crds,alpha)
%
%         crds =    [n x 2] matrix of point coordinates.
%         alpha =   percentage for 1-alpha confidence region [default = 5].
%

function kernareatest(crds,alpha)
  d = nndist(crds);
  hmin = min(d);
  hmax = max(d);

  nh = 40;
  h = linspace(hmin,hmax,nh)';

  area0 = zeros(size(h));               % Areas without adaptive smoothing
  for i = 1:nh
    area0(i)=kernarea(crds,alpha,h(i),0);
  end; 

  area1 = zeros(size(h));               % Areas with adaptive smoothing
  for i = 1:nh
    area1(i)=kernarea(crds,alpha,h(i),1);
  end; 

  close all;
  plot(h,area0,'k',h,area1,'k--');
  legend('No adaptive smoothing','Adaptive smoothing');
  putxlab('h');
  putylab('Area');

  nplots = 6;
  h = linspace(hmin,hmax,nplots);
  for i = 1:nplots
    figure;
    kernarea(crds,alpha,h(i),1,1);
    puttitle(sprintf('h = %4.2f',h(i)));
  end;

  return;
