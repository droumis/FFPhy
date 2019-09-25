% calculate mutual information of the generative model vs. bin size
%  1. from histogramm-ed spike data  2. from exact distribution
%function mutualgen

nbinx=[1:2:50];
nbiny=[1:2:50];

mx= 50;
my= 3.5;
r= -0.7;
Sx= 15;
Sy= 0.7;
Sigma=[Sx sqrt(Sx*Sy)*r; sqrt(Sx*Sy)*r Sy];

minx= 0;
maxx= 100;
miny= 0;
maxy= 2*pi;

N=10000;

NX= length(nbinx);
NY= length(nbiny);
corr=zeros(NX,NY);
corr2=zeros(NX,NY);
corrn=zeros(NX,NY);
I=zeros(NX,NY);
In=zeros(NX,NY);
%keyboard
disp('start');
for ix=1:NX
    bx=[minx:(maxx-minx)/nbinx(ix):maxx];
    for iy=1:NY
        fprintf(1,'.');
        by=[miny:(maxy-miny)/nbiny(iy):maxy];

        R=mvnrnd([mx,my],Sigma,N);
        tmp=corrcoef(R);
        corr2(ix,iy)=tmp(1,2);
        [C,BX,BY]=hist2(R,bx,by);
        %    imagesc(C)
        [mx,my,Vx,Vy,corr(ix,iy),I(ix,iy)]=Hist2Stat(C',bx,by);

        borders{1}=bx;
        borders{2}=by;
        [zn xn]=MultNormal([mx my], Sigma, borders);
        [mx,my,Vx,Vy,corrn(ix,iy),In(ix,iy)]=Hist2Stat(zn,bx,by);
    end
end
disp('done');
save mInfoDist I In corr corrn corr2 nbinx nbiny

