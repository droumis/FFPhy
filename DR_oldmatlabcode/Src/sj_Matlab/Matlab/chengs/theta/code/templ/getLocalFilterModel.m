function [model, fopts]=getLocalFilterModel(data)

minpos=min(data.linpos); % around 6cm
maxpos=max(data.linpos); % around 160cm

model.name= 'PosPhase_Isi';
model.max_isi= 0.08;
model.xp.name='CSplines2d';
for k=1:4
% spacing ca. 7.5cm (5cm was pretty good as well)
    model.xp.cpx{k}= [minpos, linspace(minpos,maxpos,20), maxpos]; 
% spacing 30 deg
    model.xp.cpy{k}= [0, linspace(0,2*pi,12), 2*pi];
end
model.xp.conv= 1.5;
model.xp.convper= .1;

model.isi.name='CSplines'; 
model.isi.cpx=[0,1:4:25,30:10:100,100]/1000;
model.isi.conv= .3;
model.isi.convper= .1;
% model.isi.name='Const'; 
% model.isi.a=1;

fopts.name= 'AscentFilter';
%fopts.eps= [1 1];
%fopts.eps= [1, 0.5 * ones(1,4)];
%fopts.eps= [5 * ones(1,16), 0];
fopts.eps= [3 * ones(1,16), .5 * ones(1,4)];
fopts.niter= 20;
fopts.maxGradient=50;
