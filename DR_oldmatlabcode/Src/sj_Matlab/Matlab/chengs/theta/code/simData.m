%function simData
% Generate artificial datasets with spikes based on the AF estimate of a real
% cell and the real behavior of the animal.

nsets= 100;
load genmodel-longISI_x1_t003
%load genmodel-longISI_x003_t03
%load genmodel-longISI_x003_t003
%load genmodel-longISI_x18_t0003
%load genmodel-longISI_x003_t1
%load genmodel-longISI_x3e-4_t003
%load genmodel-longISI_x0001_t10

load behavdata01
gen.name='RescalingGen';
gen.rand_seed= round(datenum(clock)/100+cputime);
gen.timestep= 1e-3;

select= [];
select.cellnum= [];
spikedata= {};

genmodel.isi.parameters.convertFac= genmodel.isi.parameters.convertFac*440/760;

%t= behavdata{1}{4}.time;
%dt= mean(diff(t))/2;
%nt= 2*length(t);
%data.time= [0:nt-1]'*dt + t(1);
%data.linpos= interp1(t, behavdata{1}{4}.linpos, data.time, 'spine');
%data.traj= round(interp1(t, behavdata{1}{4}.traj, data.time, 'linear'));

data= behavdata{1}{4};

nspikes= zeros(nsets, 1);
for i=1:nsets
    select.cellnum(i,:)= [1 4 27 i];
    sp= genSpikes(data, genmodel, gen);
    gen.rand_seed= gen.rand_seed+ round(rand*100);
    spikedata{1}{4}{27}{i}.time= sp.spiketimes;
    spikedata{1}{4}{27}{i}.index= sp.ispikes;
    nspikes(i)= length(sp.spiketimes);
end

save spikedata01 spikedata
save select-gen select 

fprintf(1, 'average number of spikes= %d\n', round(mean(nspikes)));
