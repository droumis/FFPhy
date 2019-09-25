clear all;
randn('state', 0);

x = 0:0.05:10; x = x';
yf = (10 * (sin(2.*pi.*0.1.*x - pi./6) + 1.1 + 5.*kaiser(length(x),150)));
y = zeros(size(yf));
sig = 6;

iknt = 4;
pr = 'Uniform';
prp = [1 50];
brn = 300;
sm = 1000;
tu = 60;
cn = 0.4;
ft = 2.*length(y);
cnf = 95;
bn = length(x);
vb = logical(0);

for ndx = 1:length(y)
    y(ndx) = normrnd(yf(ndx), sig);
end
tic;
nbars(x,y,iknt,pr,prp,brn,sm,tu,cn,ft,cnf,bn,vb);
toc
load -ascii samp_mugrid;
load -ascii samp_mu;

%xx = linspace(min(x),max(x),ft);

confint = zeros(2, size(samp_mu,2));
for ndx = 1:size(samp_mu,2)
    confint(1,ndx) = prctile(samp_mu(:,ndx),1);
    confint(2,ndx) = prctile(samp_mu(:,ndx),99);
end
fmu = mean(samp_mu,1);
plot(x,y,'.');
hyf = line(x,yf);
set(hyf,'Color','b');
set(hyf,'LineWidth',1);
% plot(x,[(yf - 2.*sig), (yf + 2.*sig)],'--r');
hfmu = line(x,fmu);
set(hfmu,'Color','g');
set(hfmu,'LineWidth',2);

hconfl = line(x, confint(1,:));
set(hconfl, 'Color', 'g');
set(hconfl, 'LineStyle', '--');
hconfu = line(x, confint(2,:));
set(hconfu, 'Color', 'g');
set(hconfu, 'LineStyle', '--');
%plot(xx,samp_mugrid(end,:),'g');hold off;
    
