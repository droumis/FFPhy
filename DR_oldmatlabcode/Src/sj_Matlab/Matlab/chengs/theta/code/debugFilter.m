%function debugFilter(rat, din, ein, tin, cin, niter, Qxp, Qt)
function debugFilter(mode, rat, din, ein, tin, cin, niter, Qxp, Qt, Q1)

model.name= 'PosPhase_Isi';
%model.name= 'Pos_Isi';

d=str2num(din); e= str2num(ein); t= str2num(tin); c= str2num(cin);

load(['/home/chengs/theta/' rat '/data2/behavdata' sprintf('%.2d', d)]);
load(['/home/chengs/theta/' rat '/data2/spikedata' sprintf('%.2d', d)]);


data= behavdata{d}{e};
data.spiketimes= spikedata{d}{e}{t}{c}.time;
data.spikeindex= spikedata{d}{e}{t}{c}.index;


%%% Tuning parameters


ncpx= 20; % spacing n=20, ca. 7.5cm (n=30, ~5cm)
ncpy= 12; % separation 30deg

ncpt= 20; % number of control points for isi- distribution


%%%

small=1e-4;

minpos=min(data.linpos)-small; % around 6cm
maxpos=max(data.linpos)+small; % around 160cm


model.xp.name='CSplines2d';
%model.xp.rectify= 1;
model.xp.periodic= 'y';
for k=1:4
    model.xp.cpx{k}= [minpos-small, linspace(minpos,maxpos,ncpx), maxpos+small]; 
    model.xp.cpy{k}= 2*pi*[0:ncpy]/ncpy;
end
model.xp.conv= 0.1;
model.xp.convper= .1;
model.xp.outputCompress= 1;
model.xp.outputInterval= 1000; % timesteps!

model.x.name='CSplines';
%model.x.rectify= 1;
for k=1:4
    model.x.cpx{k}= [minpos-small, linspace(minpos,maxpos,ncpx), maxpos+small]; 
end
odel.x.conv= 0.1;
model.x.convper= .1;
model.x.outputCompress= 1;
model.x.outputInterval= 1000; % timesteps!

model.isi.name='CSplines'; 
%model.isi.name='CSplinesNorm'; 
%model.isi.rectify= 1;
model.isi.conv= .1;
model.isi.convper= .1;
model.isi.outputCompress= 1;
model.isi.outputInterval= 1000;  % timesteps!
% model.isi.name='Const'; 
% model.isi.a=1;

totalncp= 0;
switch model.name
case 'Pos_Isi';
    for k=1:4; totalncp= totalncp+length(model.x.cpx{k}); end
    model.max_isi= 0.2;           % in sec!
    model.isi.cpx=[-small 0 1:4:25 30:20:200 200+small 201]/1000;
    fopts.eps= [str2num(Qxp)* ones(1,4), str2num(Qt) * ones(1,4)];
case 'PosPhase_Isi';
    for k=1:4; totalncp= totalncp+length(model.xp.cpx{k})*length(model.xp.cpy{k}); end
    model.max_isi= 0.08;           % in sec!
    model.isi.cpx=[-small 0 1:2:17 17:5:37 47:10:67 80+small 81]/1000;
%    fopts.eps= [2* ones(1,16), 1 * ones(1,4)];
%    fopts.eps= [0.15* ones(1,16), 0.15 * ones(1,4)];
    fopts.eps= [str2num(Qxp)* ones(1,16), str2num(Qt) * ones(1,4)];
end
nisi= length(model.isi.cpx);

%totalncp+ nisi
%totalncp
%nisi

switch mode
case 'AF'
    fopts.name= 'AscentFilter';
    fopts.alternatePass= 1;
case 'AS'
    fopts.name= 'AFSmooth';
    fopts.alternatePass= 1;
case 'EM'
    fopts.name= 'EM';
case 'KS'
    fopts.name= 'KalmanSmooth';
case 'KB'
    fopts.name= 'KSBackFirst';
end
fopts.maxGradient=50;
%fopts.eps= [5* ones(1,16), 0.5 * ones(1,4)];
fopts.niter= str2num(niter);

switch mode
case {'AF', 'AS'}
    model.xp.operator.name= 'Rectify';
    model.x.operator.name= 'Rectify';
    model.isi.operator.name= 'Rectify';
otherwise
%    model.xp.operator.name= 'Rectify';
%    model.x.operator.name= 'Rectify';
%    model.isi.operator.name= 'Rectify';
    model.xp.operator.name= 'Exp';
    model.x.operator.name= 'Exp';
    model.isi.operator.name= 'Exp';
    fopts.Q1= str2num(Q1);
    fopts.Qt= str2num(Qxp)*ones(1,totalncp);
    fopts.Qt(1,totalncp+1:totalncp+nisi)= str2num(Qt)*ones(1,nisi);
end


%[model, filter]= adaptFilter(data, model, fopts);
model= adaptFilter(data, model, fopts);

%keyboard
%load(['adaptest' sprintf('%.2d', d)]);
%model= adaptest{d}{e}{t}{c}.model;

%load /bach/theta/AdaptFilter/newvis/colormap;
%load /bach/theta/AdaptFilter/newvis/colormap2;
%opt.cmap= cmap;
%opt.playrate= 0.0003; % nearly realtime
%opt.playrate= 0.003;
%opt.cmap= jet(1024);

%keyboard

%visModel(data, model, opt)
%disp('Done.');

