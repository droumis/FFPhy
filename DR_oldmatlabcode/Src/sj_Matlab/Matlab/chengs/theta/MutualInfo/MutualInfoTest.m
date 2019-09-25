% mutual information vs. linear correlation coefficient r for Gaussian PDF

scale=4;
N=10000;
bx=scale*[-4:0.1:4];
by=scale*[-4:0.1:4];

r=[-1:0.01:1];
nj=length(r);
corr=zeros(nj,1);
corr2=zeros(nj,1);
I=zeros(nj,1);
In=zeros(nj,1);

a=0.9;
for j=1:nj
    Sigma=[1 sqrt(a)*r(j); sqrt(a)*r(j) a]*scale^2;
    R=mvnrnd([0,0],Sigma,N);
    tmp=corrcoef(R);
    corr2(j)=tmp(1,2);
    C=hist2d(R,bx,by);
%    imagesc(C)
    [mx,my,Vx,Vy,corr(j),I(j)]=Hist2Stat(C',bx,by);

    if abs(r(j))~=1
	borders{1}=bx;
	borders{2}=by;
	[zn xn]=MultNormal([0 0], Sigma, borders);
	[mx,my,Vx,Vy,corrn(j),In(j)]=Hist2Stat(zn,bx,by);
    end
end

figure(1); clf; 
plot(r,I)
xlabel('r'); ylabel('mutual info from generated spikes');
title('Gaussian PDF');

figure(2); clf; plot(r,corr)
figure(3); clf; plot(corr,corr2,'.','MarkerSize',3)

figure(4); clf; 
plot(r,In)
xlabel('r'); ylabel('mutual info from exact distribution');
title('Gaussian PDF');
axis tight
sen_styles
exportfig(gcf,'MutualInfo_R.eps',EPSOpts);

