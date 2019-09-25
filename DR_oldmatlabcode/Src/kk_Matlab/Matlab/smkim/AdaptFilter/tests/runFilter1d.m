function runFilter1d

load(['dynamic_data']);

model.name= 'Pos_Isi';
model.max_isi= 0.08;
model.spatial.name='CSplines';
model.spatial.cpx= [-5, linspace(0,150,30), 155 160];
model.spatial.conv= 0.1;
model.spatial.convper= .1;
% model.spatial.name='Const';
% model.spatial.a=2;

model.isi.name='CSplines'; 
model.isi.cpx=[-0.1, linspace(0,100,10), 100.1]/1000;
model.isi.conv= .3;
model.isi.convper= .1;
%model.isi.name='Const'; 
%model.isi.a=1;
 
save(['dynamic1d_init'], 'model');

fopts.name= 'AscentFilter';
fopts.alternatePass= 1;
fopts.eps= [2 * ones(1,8)];
%fopts.eps= [2 * ones(1,4), 0.1*ones(1,4)];
fopts.niter= 20;
fopts.maxGradient=50;

result= adaptFilter(data, model, fopts);

%visData(data,model.spatial, vis_opts);
 
save(['dynamic1d_result'], 'result');

