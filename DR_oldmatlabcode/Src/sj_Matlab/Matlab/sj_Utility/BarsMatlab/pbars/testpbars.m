clear all;
randn('state', 0);

% load -ascii example.data;
% x = example(2:end,1);
% y = int32(example(2:end,2));
% yf = example(2:end,2);
x = linspace(0,10,200); x = x';
yf = (10 * (sin(2.*pi.*0.1.*x - pi./6) + 2 + 5.*kaiser(length(x),150)));
y = zeros(size(yf));
sig = 6;
for ndx = 1:length(y)
    y(ndx) = poissrnd(yf(ndx));
end
y = int32(y);
y(find(y < 0)) = 0;

%     save y y;


lgsp = true;
iknt = 5;
pr = 'Uniform';
prp = [1 15];
brn = 100;
sm = 1000;
tu = 50;
cn = 0.4;
ft = 2.*length(y);
cnf = 95;
trl = 60;
bn = length(x);
vb = true;


tic;
pbars(x,y,lgsp,iknt,pr,prp,brn,sm,tu,cn,ft,cnf,trl,bn,vb);
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
clear pbars.dll
