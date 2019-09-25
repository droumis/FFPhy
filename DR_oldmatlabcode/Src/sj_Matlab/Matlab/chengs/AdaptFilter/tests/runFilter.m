function runFilter

load(['dynamic_data']);


%%%%%%%% 

fopts.name= 'AscentFilter';
fopts.maxGradient=50;
fopts.eps= [.15* ones(1,16), 0.05 * ones(1,4)];
fopts.alternatePass= 1;
%fopts.eps= [.15* ones(1,16), 0.05 * ones(1,16)];
%fopts.alternatePass= 0;
fopts.niter= 20;
%fopts.niter= 10;

ncpx= 20; % spacing n=20, ca. 7.5cm (n=30, ~5cm)
ncpy= 12; % separation 30deg

max_isi= 0.080;           % in sec!

small=1e-4;


minpos=min(data.linpos)-small; % around 6cm
maxpos=max(data.linpos)+small; % around 160cm

model.name= 'PosPhase_Isi';
model.max_isi= max_isi;           % in sec!

model.spatial.name='CSplines2d';
model.spatial.periodic= 'y';
model.spatial.cpx= [minpos-small, linspace(minpos,maxpos,ncpx), maxpos+small]; 
model.spatial.cpy= 2*pi*[0:ncpy]/ncpy;
model.spatial.conv= 0.3;
model.spatial.convper= .1;
model.spatial.outputCompress= 1;
model.spatial.outputInterval= 1000; % timesteps!

model.isi.name='CSplines'; 
%model.isi.name='CSplinesNorm'; 
%model.isi.cpx=[-small 0 1:2:25 30:5:45  50:10:70 80+small 81]/1000;
model.isi.cpx=[-small 0 1:4:25 30:10:70 80+small 81]/1000;
model.isi.conv= .1;
model.isi.convper= .1;
model.isi.outputCompress= 1;
model.isi.outputInterval= 1000; % timesteps!

%%%%%%%%

 
save(['dynamic_init'], 'model');

result= adaptFilter(data,model, fopts);
result.max_isi= model.max_isi;
%result= adaptFilterIni(data, genmodel,model, fopts);

save(['dynamic_result'], 'result');

