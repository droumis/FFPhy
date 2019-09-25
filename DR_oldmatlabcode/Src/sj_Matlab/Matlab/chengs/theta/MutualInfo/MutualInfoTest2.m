% mutual information as function of 'bin width'


n=[ceil(logspace(1,3))];
nN=length(n);
corr=zeros(nN,1);
I=zeros(nN,1);

for iN=1:nN
    z=eye(n(iN));
%    z=zeros(n(iN)); z(1,:)=ones(1,n(iN));
    [mx,my,Vx,Vy,corr(iN),I(iN)]=Hist2Stat(z,0:n(iN),0:n(iN),n(iN));
end

figure(1); clf; 

subplot(1,3,1);
plot(n,I,'.')
xlabel('n'); ylabel('mutual info');
hold on;
plot(n,log2(n),'r');
legend('numerical','theor.',4);

subplot(1,3,2);
plot(n,corr)
xlabel('n'); ylabel('r');

subplot(1,3,3);
imagesc(eye(20));
xlabel('x'); ylabel('y');
title('p(x,y)');

sen_styles
exportfig(gcf,'MutualInfo_BinWidth.eps',EPSCWideOpts);

