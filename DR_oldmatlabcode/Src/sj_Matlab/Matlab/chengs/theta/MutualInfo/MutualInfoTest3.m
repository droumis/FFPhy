% mutual information for "lines"


m=100;
x=0:m;
y=0:m;

% rotating  one line
b1=1;
a=[0,0.01,0.05,.1,.2,.5,.8,.9,1:.2:4];

na=length(a);
ra=zeros(na,1);
Ia=zeros(na,1);
for ia=1:na
    z=zeros(m);
    ytmp=a(ia)*x+b1;
    im=1;
    while im <= m & ceil(ytmp(im)) <= m
	z(im,ceil(ytmp(im)))=1;
%    z(im,ceil(ytmp(im))+1)=1;
	im=im+1;
    end
    [mx,my,Vx,Vy,ra(ia),Ia(ia)]=Hist2Stat(z,x,y);
end

figure(1); clf

subplot(3,2,1)
plot(a,Ia,'k-')
hold on
plot(a,ra,'k:')
hold off
subplot(3,2,2);
imagesc(z'); axis xy;

% moving two lines apart
a=0.1;
b1=1;
b2=1:20;

nb=length(b2);
rb=zeros(nb,1);
Ib=zeros(nb,1);
ytmp=a*x+b1;
for ib=1:nb
    z=zeros(m);
    im=1;
    while im <= m & ceil(ytmp(im)) <= m
	z(im,ceil(ytmp(im)))=1;
	if ceil(ytmp(im))+b2(ib) <= m
	    z(im,ceil(ytmp(im))+b2(ib))=1;
	end
	im=im+1;
    end
    [mx,my,Vx,Vy,rb(ib),Ib(ib)]=Hist2Stat(z,x,y);
end


subplot(3,2,3)
plot(b2,Ib,'k-')
hold on
plot(b2,rb,'k:')
hold off
subplot(3,2,4);
imagesc(z'); axis xy;


% thickening a line
a=0.1;
b1=1;
b2=0:20;

nb=length(b2);
rt=zeros(nb,1);
It=zeros(nb,1);
ytmp=a*x+b1;
z=zeros(m);
for ib=1:nb
    im=1;
    while im <= m & ceil(ytmp(im)) <= m
	z(im,ceil(ytmp(im)))=1;
	if ceil(ytmp(im))+b2(ib) <= m
	    z(im,ceil(ytmp(im))+b2(ib))=1;
	end
	im=im+1;
    end
    [mx,my,Vx,Vy,rt(ib),It(ib)]=Hist2Stat(z,x,y);
end


subplot(3,2,5)
plot(b2,It,'k-')
hold on
plot(b2,rt,'k:')
hold off
subplot(3,2,6);
imagesc(z'); axis xy;

sen_styles
exportfig(gcf,'MutualInfo_lines.eps',EPSCLargeOpts);


% non-linearity, power
p=.5:.1:4;
np=length(p);
rp=zeros(np,1);
Ip=zeros(np,1);
for ip=1:np
    z=zeros(m);
    ytmp=(x/m).^p(ip);
    im=1;
    for im=1:m
	z(im,1+floor(m*ytmp(im)))=1;
    end
    [mx,my,Vx,Vy,rp(ip),Ip(ip)]=Hist2Stat(z,x,y);
end

figure(2); clf
subplot(3,2,1)
plot(p,Ip,'k-')
hold on
plot(p,rp,'k:')
hold off
subplot(3,2,2);
imagesc(z'); axis xy;


% non-linearity, circle
f=0:.2:2*pi+.2;
R=40;
nf=length(f)-1;
rf=zeros(nf,1);
If=zeros(nf,1);
z=zeros(m);
for iF=1:nf
    fStart=f(iF); fEnd=f(iF+1);
    for k=0:20
	angle=(fEnd-fStart)/20*k+fStart;
	z(ceil(R*cos(angle)+m/2), ceil(R*sin(angle)+m/2))=1;
    end
    [mx,my,Vx,Vy,rf(iF),If(iF)]=Hist2Stat(z,x,y);
end

subplot(3,2,3)
plot(f(1:nf),If,'k-')
hold on
plot(f(1:nf),rf,'k:')
hold off
subplot(3,2,4);
imagesc(z'); axis xy;

exportfig(gcf,'MutualInfo_nonlinear.eps',EPSCLargeOpts);

