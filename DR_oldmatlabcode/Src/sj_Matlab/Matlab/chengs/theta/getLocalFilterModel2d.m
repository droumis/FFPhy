function [model, fopts]=getLocalFilterModel(data)

%%% Tuning parameters

fopts.alternatePass= 1;
%fopts.eps= [.15* ones(1,16), 0.05 * ones(1,16)];
%fopts.alternatePass= 0;
fopts.niter= 20;

ncpx= 30; % spacing n=20, ca. 7.5cm (n=30, ~5cm)
ncpy= 12; % separation 30deg

max_isi= 0.080;           % in sec!


%%% relatively fixed parameters

small=1e-4;

fopts.maxGradient=1e4;

minpos=min(data.linpos)-small; % around 6cm
maxpos=max(data.linpos)+small; % around 160cm

model.name= 'PosPhase_Isi';
model.max_isi= max_isi;           % in sec!

model.spatial.name='CSplines2d';
model.spatial.periodic= 'y';
model.spatial.operator.name= 'Rectify';
for k=1:4
    model.spatial.cpx{k}= [minpos-small, linspace(minpos,maxpos,ncpx), maxpos+small]; 
    model.spatial.cpy{k}= 2*pi*[0:ncpy]/ncpy;
end
model.spatial.conv= 0.005;
model.spatial.convper= .1;
model.spatial.outputCompress= 1;
model.spatial.outputInterval= 1000; % timesteps!

model.isi.name='CSplines'; 
model.isi.operator.name= 'Rectify';
%model.isi.name='CSplinesNorm'; 
%model.isi.cpx=[-small 0 1:4:25 30:10:70 80+small 81]/1000;
model.isi.cpx=[-small 1 3 5 7 11 15:10:65 80+small 81]/1000;
model.isi.conv= .005;
model.isi.convper= .1;
model.isi.outputCompress= 1;
model.isi.outputInterval= 1000; % timesteps!

fopts.name= 'AscentFilter';
