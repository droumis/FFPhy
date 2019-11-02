%%%%run decoder_prep_bond04.m first

for i=1:length(vecL0)
    if vecL0(i)>=3
        vecL(i)=6-vecL0(i);
    elseif vecL0(i)<=-3
        vecL(i)=-6-vecL0(i);
    else
        vecL(i)=vecL0(i);
    end
end
vecLF0=[vecT vecL'];
[d1,d2]=sort(vecLF0(:,1));
vecLF=vecLF0(d2,:); %column 1: ascending time; column 2: linearized position


A=pos{ex}{ep}.data(:,1); %time stamps for animal's trajectory
ti=round(A(1)*1000):1:round(A(end)*1000); %binning time stamps at 1 ms

stateV=linspace(-3,3,61);

%%
[sn,state_bin]=histc(vecLF(:,2),stateV);
state_disM=[state_bin(1:end-1) state_bin(2:end)];
n=size(stateV,2);
%by column=departuring
for s=1:n
    sp0=state_disM(find(state_disM(:,1)==s),2); %by departure x_k-1 (by column); sp0 is the departuring x_(k-1);
    if isempty(sp0)==0
        stateM(:,s)=histc(sp0,linspace(1,n,n))./size(sp0,1);
    elseif isempty(sp0)==1
        stateM(:,s)=zeros(1,n);
    end
end
K=inline('exp(-(x.^2+y.^2)/2/sig^2)'); %gaussian
[dx,dy]=meshgrid([-1:1]);
sig=0.5;
weight=K(sig,dx,dy)/sum(sum(K(sig,dx,dy))); %normalizing weights
stateM_gaus=conv2(stateM,weight,'same'); %gaussian smoothed
stateM_gausnorm=stateM_gaus*diag(1./sum(stateM_gaus,1)); %normalized to confine probability to 1




%%
segT_I1=[2520 2565;2647 2705;2748 2777;2835 2856;2880 2909;2959 2975;3023.5 3032;3078 3097;3139 3179;3241.5 3255.5;3316.4 3330.6];
segT_I0=[2585 2620;2705 2730;2795 2820;2858 2880;2909 2936;2991 3023.5;3037 3062;3113.5 3129;3190.5 3227.5;3274.5 3301;3340 3387.5];

ind_I1=[];
for i=1:size(segT_I1,1)
    ind_I1a=find(vecLF(:,1)>=segT_I1(i,1)&vecLF(:,1)<=segT_I1(i,2));
    ind_I1=[ind_I1;ind_I1a];
end
ind_I0=[];
for i=1:size(segT_I0,1)
    ind_I0a=find(vecLF(:,1)>=segT_I0(i,1)&vecLF(:,1)<=segT_I0(i,2));
    ind_I0=[ind_I0;ind_I0a];
end

%%
n=length(stateV);
stateM_I1=zeros(n,n);
for i=1:size(segT_I1,1)
    vecLF_seg=vecLF(find(vecLF(:,1)>=segT_I1(i,1)&vecLF(:,1)<=segT_I1(i,2)),:);
    [sn,state_bin]=histc(vecLF_seg(:,2),stateV);
    state_disM=[state_bin(1:end-1) state_bin(2:end)];
    stateM_seg=zeros(n,n);
    for s=1:n
        sp0=state_disM(find(state_disM(:,1)==s),2); %by departure x_k-1 (by column); sp0 is the departuring x_(k-1);
        if isempty(sp0)==0
            stateM_seg(:,s)=histc(sp0,linspace(1,n,n))./size(sp0,1);
        elseif isempty(sp0)==1
            stateM_seg(:,s)=zeros(1,n);
        end
    end
    stateM_I1=stateM_I1+stateM_seg;
end
stateM_I1=stateM_I1*diag(1./sum(stateM_I1,1));
K=inline('exp(-(x.^2+y.^2)/2/sig^2)'); %gaussian
[dx,dy]=meshgrid([-1:1]);
sig=0.5;
weight=K(sig,dx,dy)/sum(sum(K(sig,dx,dy))); %normalizing weights
stateM_gaus=conv2(stateM_I1,weight,'same'); %gaussian smoothed
stateM_I1_gausnorm=stateM_gaus*diag(1./sum(stateM_gaus,1)); %normalized to confine probability to 1

stateM_I0=zeros(n,n);
for i=1:size(segT_I0,1)
    vecLF_seg=vecLF(find(vecLF(:,1)>=segT_I0(i,1)&vecLF(:,1)<=segT_I0(i,2)),:);
    [sn,state_bin]=histc(vecLF_seg(:,2),stateV);
    state_disM=[state_bin(1:end-1) state_bin(2:end)];
    stateM_seg=zeros(n,n);
    for s=1:n
        sp0=state_disM(find(state_disM(:,1)==s),2); %by departure x_k-1 (by column); sp0 is the departuring x_(k-1);
        if isempty(sp0)==0
            stateM_seg(:,s)=histc(sp0,linspace(1,n,n))./size(sp0,1);
        elseif isempty(sp0)==1
            stateM_seg(:,s)=zeros(1,n);
        end
    end
    stateM_I0=stateM_I0+stateM_seg;
end
stateM_I0=stateM_I0*diag(1./sum(stateM_I0,1));
K=inline('exp(-(x.^2+y.^2)/2/sig^2)'); %gaussian
[dx,dy]=meshgrid([-1:1]);
sig=0.5;
weight=K(sig,dx,dy)/sum(sum(K(sig,dx,dy))); %normalizing weights
stateM_gaus=conv2(stateM_I0,weight,'same'); %gaussian smoothed
stateM_I0_gausnorm=stateM_gaus*diag(1./sum(stateM_gaus,1)); %normalized to confine probability to 1





%%
ripplescons = loaddatastruct('/opt/data13/kkay/Bon/','bon','ripplescons',ex);
startT=ripplescons{ex}{ep}{1}.starttime;
endT=ripplescons{ex}{ep}{1}.endtime;
ripple_seg=[round(startT*1000)-ti(1)-1 round(endT*1000)-ti(1)-1]; %index for ripple segments


for kk=1:size(a,2)
    j=a(kk);i=b(kk);
    B=spikes{ex}{ep}{j}{i}.data(:,1); %spiking times for tetrode j, cell i
    xi=round(B*1000); %binning spiking times at 1 ms
    [sptrain2,~]=ismember(ti,xi); %sptrain2: spike train binned at 1 ms instead of 33.4ms (sptrain0)
    sptrain2_list{kk}=sptrain2;
end

for k=1:size(ripple_seg,1)
    spike_r=[];
    for kk=1:size(a,2)
        sptrain2=sptrain2_list{kk};
        spike_r=[spike_r;sptrain2(ripple_seg(k,1):ripple_seg(k,2))];
    end
    spike_r_all{k}=spike_r;
end


for k=1:size(ripple_seg,1)
    spike_r=spike_r_all{k};
    sumR(k)=sum(spike_r(:));
end
rippleI=find(sumR>10);
%length(kk)



%%
%%%%choose ripple index
%%%%set stdev value and number of bins for numeric integral

rIndV=60;
spike_r=spike_r_all{rippleI(rIndV)};

%%





%%
linpos = loaddatastruct('/opt/data13/kkay/Bon/','bon','linpos',ex);
tlin=linpos{4}{2}.statematrix.time;

poslin=vecLF(:,2);xdel=0.1; 
xs=min(poslin):xdel:max(poslin);
dt=tlin(2)-tlin(1);
xtrain=poslin';


%%
%%%%encode per tetrode
load('bond_data\bond04\01-149\bond04-01_params.mat');
%ind_t1=find(filedata.params(:,1)/10000>=tlin(1)&filedata.params(:,1)/10000<=tlin(end)&filedata.params(:,8)>10&(filedata.params(:,2)>100|filedata.params(:,3)>100|filedata.params(:,4)>100|filedata.params(:,5)>100));
%threshold at 100uV
ind_t1=find(filedata.params(:,1)/10000>=tlin(1)&filedata.params(:,1)/10000<=tlin(end));
time_t1=filedata.params(ind_t1,1);
mark0_t1=[filedata.params(ind_t1,2) filedata.params(ind_t1,3) filedata.params(ind_t1,4) filedata.params(ind_t1,5)];
time2_t1=time_t1/10000;
spikeT0_t1=time2_t1;
[procInd0_t1,procInd1_t1]=histc(spikeT0_t1,tlin);
procInd_t1=find(procInd0_t1);
spikeT_t1=tlin(procInd_t1);
spike_t1=procInd0_t1';
[~,rawInd0_t1]=histc(spikeT0_t1,time2_t1);
markAll_t1(:,1)=procInd1_t1;markAll_t1(:,2:5)=mark0_t1(rawInd0_t1(rawInd0_t1~=0),:);
mdel=20; ms=min(min(markAll_t1(:,2:5))):mdel:max(max(markAll_t1(:,2:5))); 
sxker=2*xdel; smker = mdel; T=size(tlin,1);
occ=normpdf(xs'*ones(1,T),ones(length(xs),1)*xtrain,sxker)*ones(T,length(ms));
%occ: columns are identical; occupancy based on position; denominator
Xnum_t1=normpdf(xs'*ones(1,length(spikeT0_t1)),ones(length(xs),1)*xtrain(procInd1_t1),sxker);
%Xnum: Gaussian kernel estimators for position
Lint_t1=sum(Xnum_t1,2)./occ(:,1)./dt; %integral
Lint_t1=Lint_t1./sum(Lint_t1);

load('bond_data\bond04\02-151\bond04-02_params.mat');
%ind_t2=find(filedata.params(:,1)/10000>=tlin(1)&filedata.params(:,1)/10000<=tlin(end)&filedata.params(:,8)>10&(filedata.params(:,2)>100|filedata.params(:,3)>100|filedata.params(:,4)>100|filedata.params(:,5)>100));
ind_t2=find(filedata.params(:,1)/10000>=tlin(1)&filedata.params(:,1)/10000<=tlin(end));
time_t2=filedata.params(ind_t2,1);
mark0_t2=[filedata.params(ind_t2,2) filedata.params(ind_t2,3) filedata.params(ind_t2,4) filedata.params(ind_t2,5)];
time2_t2=time_t2/10000;
spikeT0_t2=time2_t2;
[procInd0_t2,procInd1_t2]=histc(spikeT0_t2,tlin);
procInd_t2=find(procInd0_t2);
spikeT_t2=tlin(procInd_t2);
spike_t2=procInd0_t2';
[~,rawInd0_t2]=histc(spikeT0_t2,time2_t2);
markAll_t2(:,1)=procInd1_t2;markAll_t2(:,2:5)=mark0_t2(rawInd0_t2(rawInd0_t2~=0),:);
mdel=20; ms=min(min(markAll_t2(:,2:5))):mdel:max(max(markAll_t2(:,2:5))); 
sxker=2*xdel; smker = mdel; T=size(tlin,1);
occ=normpdf(xs'*ones(1,T),ones(length(xs),1)*xtrain,sxker)*ones(T,length(ms));
%occ: columns are identical; occupancy based on position; denominator
Xnum_t2=normpdf(xs'*ones(1,length(spikeT0_t2)),ones(length(xs),1)*xtrain(procInd1_t2),sxker);
%Xnum: Gaussian kernel estimators for position
Lint_t2=sum(Xnum_t2,2)./occ(:,1)./dt; %integral
Lint_t2=Lint_t2./sum(Lint_t2);

load('bond_data\bond04\04-107\bond04-04_params.mat');
%ind_t4=find(filedata.params(:,1)/10000>=tlin(1)&filedata.params(:,1)/10000<=tlin(end)&filedata.params(:,8)>10&(filedata.params(:,2)>100|filedata.params(:,3)>100|filedata.params(:,4)>100|filedata.params(:,5)>100));
ind_t4=find(filedata.params(:,1)/10000>=tlin(1)&filedata.params(:,1)/10000<=tlin(end));
time_t4=filedata.params(ind_t4,1);
mark0_t4=[filedata.params(ind_t4,2) filedata.params(ind_t4,3) filedata.params(ind_t4,4) filedata.params(ind_t4,5)];
time2_t4=time_t4/10000;
spikeT0_t4=time2_t4;
[procInd0_t4,procInd1_t4]=histc(spikeT0_t4,tlin);
procInd_t4=find(procInd0_t4);
spikeT_t4=tlin(procInd_t4);
spike_t4=procInd0_t4';
[~,rawInd0_t4]=histc(spikeT0_t4,time2_t4);
markAll_t4(:,1)=procInd1_t4;markAll_t4(:,2:5)=mark0_t4(rawInd0_t4(rawInd0_t4~=0),:);
mdel=20; ms=min(min(markAll_t4(:,2:5))):mdel:max(max(markAll_t4(:,2:5))); 
sxker=2*xdel; smker = mdel; T=size(tlin,1);
occ=normpdf(xs'*ones(1,T),ones(length(xs),1)*xtrain,sxker)*ones(T,length(ms));
%occ: columns are identical; occupancy based on position; denominator
Xnum_t4=normpdf(xs'*ones(1,length(spikeT0_t4)),ones(length(xs),1)*xtrain(procInd1_t4),sxker);
%Xnum: Gaussian kernel estimators for position
Lint_t4=sum(Xnum_t4,2)./occ(:,1)./dt; %integral
Lint_t4=Lint_t4./sum(Lint_t4);

load('bond_data\bond04\05-105\bond04-05_params.mat');
%ind_t5=find(filedata.params(:,1)/10000>=tlin(1)&filedata.params(:,1)/10000<=tlin(end)&filedata.params(:,8)>10&(filedata.params(:,2)>100|filedata.params(:,3)>100|filedata.params(:,4)>100|filedata.params(:,5)>100));
ind_t5=find(filedata.params(:,1)/10000>=tlin(1)&filedata.params(:,1)/10000<=tlin(end));
time_t5=filedata.params(ind_t5,1);
mark0_t5=[filedata.params(ind_t5,2) filedata.params(ind_t5,3) filedata.params(ind_t5,4) filedata.params(ind_t5,5)];
time2_t5=time_t5/10000;
spikeT0_t5=time2_t5;
[procInd0_t5,procInd1_t5]=histc(spikeT0_t5,tlin);
procInd_t5=find(procInd0_t5);
spikeT_t5=tlin(procInd_t5);
spike_t5=procInd0_t5';
[~,rawInd0_t5]=histc(spikeT0_t5,time2_t5);
markAll_t5(:,1)=procInd1_t5;markAll_t5(:,2:5)=mark0_t5(rawInd0_t5(rawInd0_t5~=0),:);
mdel=20; ms=min(min(markAll_t5(:,2:5))):mdel:max(max(markAll_t5(:,2:5))); 
sxker=2*xdel; smker = mdel; T=size(tlin,1);
occ=normpdf(xs'*ones(1,T),ones(length(xs),1)*xtrain,sxker)*ones(T,length(ms));
%occ: columns are identical; occupancy based on position; denominator
Xnum_t5=normpdf(xs'*ones(1,length(spikeT0_t5)),ones(length(xs),1)*xtrain(procInd1_t5),sxker);
%Xnum: Gaussian kernel estimators for position
Lint_t5=sum(Xnum_t5,2)./occ(:,1)./dt; %integral
Lint_t5=Lint_t5./sum(Lint_t5);

load('bond_data\bond04\07-162\bond04-07_params.mat');
%ind_t7=find(filedata.params(:,1)/10000>=tlin(1)&filedata.params(:,1)/10000<=tlin(end)&filedata.params(:,8)>10&(filedata.params(:,2)>100|filedata.params(:,3)>100|filedata.params(:,4)>100|filedata.params(:,5)>100));
ind_t7=find(filedata.params(:,1)/10000>=tlin(1)&filedata.params(:,1)/10000<=tlin(end));
time_t7=filedata.params(ind_t7,1);
mark0_t7=[filedata.params(ind_t7,2) filedata.params(ind_t7,3) filedata.params(ind_t7,4) filedata.params(ind_t7,5)];
time2_t7=time_t7/10000;
spikeT0_t7=time2_t7;
[procInd0_t7,procInd1_t7]=histc(spikeT0_t7,tlin);
procInd_t7=find(procInd0_t7);
spikeT_t7=tlin(procInd_t7);
spike_t7=procInd0_t7';
[~,rawInd0_t7]=histc(spikeT0_t7,time2_t7);
markAll_t7(:,1)=procInd1_t7;markAll_t7(:,2:5)=mark0_t7(rawInd0_t7(rawInd0_t7~=0),:);
mdel=20; ms=min(min(markAll_t7(:,2:5))):mdel:max(max(markAll_t7(:,2:5))); 
sxker=2*xdel; smker = mdel; T=size(tlin,1);
occ=normpdf(xs'*ones(1,T),ones(length(xs),1)*xtrain,sxker)*ones(T,length(ms));
%occ: columns are identical; occupancy based on position; denominator
Xnum_t7=normpdf(xs'*ones(1,length(spikeT0_t7)),ones(length(xs),1)*xtrain(procInd1_t7),sxker);
%Xnum: Gaussian kernel estimators for position
Lint_t7=sum(Xnum_t7,2)./occ(:,1)./dt; %integral
Lint_t7=Lint_t7./sum(Lint_t7);

load('bond_data\bond04\10-146\bond04-10_params.mat');
%ind_t10=find(filedata.params(:,1)/10000>=tlin(1)&filedata.params(:,1)/10000<=tlin(end)&filedata.params(:,8)>10&(filedata.params(:,2)>100|filedata.params(:,3)>100|filedata.params(:,4)>100|filedata.params(:,5)>100));
ind_t10=find(filedata.params(:,1)/10000>=tlin(1)&filedata.params(:,1)/10000<=tlin(end));
time_t10=filedata.params(ind_t10,1);
mark0_t10=[filedata.params(ind_t10,2) filedata.params(ind_t10,3) filedata.params(ind_t10,4) filedata.params(ind_t10,5)];
time2_t10=time_t10/10000;
spikeT0_t10=time2_t10;
[procInd0_t10,procInd1_t10]=histc(spikeT0_t10,tlin);
procInd_t10=find(procInd0_t10);
spikeT_t10=tlin(procInd_t10);
spike_t10=procInd0_t10';
[~,rawInd0_t10]=histc(spikeT0_t10,time2_t10);
markAll_t10(:,1)=procInd1_t10;markAll_t10(:,2:5)=mark0_t10(rawInd0_t10(rawInd0_t10~=0),:);
mdel=20; ms=min(min(markAll_t10(:,2:5))):mdel:max(max(markAll_t10(:,2:5))); 
sxker=2*xdel; smker = mdel; T=size(tlin,1);
occ=normpdf(xs'*ones(1,T),ones(length(xs),1)*xtrain,sxker)*ones(T,length(ms));
%occ: columns are identical; occupancy based on position; denominator
Xnum_t10=normpdf(xs'*ones(1,length(spikeT0_t10)),ones(length(xs),1)*xtrain(procInd1_t10),sxker);
%Xnum: Gaussian kernel estimators for position
Lint_t10=sum(Xnum_t10,2)./occ(:,1)./dt; %integral
Lint_t10=Lint_t10./sum(Lint_t10);

load('bond_data\bond04\11-096\bond04-11_params.mat');
%ind_t11=find(filedata.params(:,1)/10000>=tlin(1)&filedata.params(:,1)/10000<=tlin(end)&filedata.params(:,8)>10&(filedata.params(:,2)>100|filedata.params(:,3)>100|filedata.params(:,4)>100|filedata.params(:,5)>100));
ind_t11=find(filedata.params(:,1)/10000>=tlin(1)&filedata.params(:,1)/10000<=tlin(end));
time_t11=filedata.params(ind_t11,1);
mark0_t11=[filedata.params(ind_t11,2) filedata.params(ind_t11,3) filedata.params(ind_t11,4) filedata.params(ind_t11,5)];
time2_t11=time_t11/10000;
spikeT0_t11=time2_t11;
[procInd0_t11,procInd1_t11]=histc(spikeT0_t11,tlin);
procInd_t11=find(procInd0_t11);
spikeT_t11=tlin(procInd_t11);
spike_t11=procInd0_t11';
[~,rawInd0_t11]=histc(spikeT0_t11,time2_t11);
markAll_t11(:,1)=procInd1_t11;markAll_t11(:,2:5)=mark0_t11(rawInd0_t11(rawInd0_t11~=0),:);
mdel=20; ms=min(min(markAll_t11(:,2:5))):mdel:max(max(markAll_t11(:,2:5)));
sxker=2*xdel; smker = mdel; T=size(tlin,1);
occ=normpdf(xs'*ones(1,T),ones(length(xs),1)*xtrain,sxker)*ones(T,length(ms));
%occ: columns are identical; occupancy based on position; denominator
Xnum_t11=normpdf(xs'*ones(1,length(spikeT0_t11)),ones(length(xs),1)*xtrain(procInd1_t11),sxker);
%Xnum: Gaussian kernel estimators for position
Lint_t11=sum(Xnum_t11,2)./occ(:,1)./dt; %integral
Lint_t11=Lint_t11./sum(Lint_t11);

load('bond_data\bond04\12-136\bond04-12_params.mat');
%ind_t12=find(filedata.params(:,1)/10000>=tlin(1)&filedata.params(:,1)/10000<=tlin(end)&filedata.params(:,8)>10&(filedata.params(:,2)>100|filedata.params(:,3)>100|filedata.params(:,4)>100|filedata.params(:,5)>100));
ind_t12=find(filedata.params(:,1)/10000>=tlin(1)&filedata.params(:,1)/10000<=tlin(end));
time_t12=filedata.params(ind_t12,1);
mark0_t12=[filedata.params(ind_t12,2) filedata.params(ind_t12,3) filedata.params(ind_t12,4) filedata.params(ind_t12,5)];
time2_t12=time_t12/10000;
spikeT0_t12=time2_t12;
[procInd0_t12,procInd1_t12]=histc(spikeT0_t12,tlin);
procInd_t12=find(procInd0_t12);
spikeT_t12=tlin(procInd_t12);
spike_t12=procInd0_t12';
[~,rawInd0_t12]=histc(spikeT0_t12,time2_t12);
markAll_t12(:,1)=procInd1_t12;markAll_t12(:,2:5)=mark0_t12(rawInd0_t12(rawInd0_t12~=0),:);
mdel=20; ms=min(min(markAll_t12(:,2:5))):mdel:max(max(markAll_t12(:,2:5)));
sxker=2*xdel; smker = mdel; T=size(tlin,1);
occ=normpdf(xs'*ones(1,T),ones(length(xs),1)*xtrain,sxker)*ones(T,length(ms));
%occ: columns are identical; occupancy based on position; denominator
Xnum_t12=normpdf(xs'*ones(1,length(spikeT0_t12)),ones(length(xs),1)*xtrain(procInd1_t12),sxker);
%Xnum: Gaussian kernel estimators for position
Lint_t12=sum(Xnum_t12,2)./occ(:,1)./dt; %integral
Lint_t12=Lint_t12./sum(Lint_t12);

load('bond_data\bond04\13-094\bond04-13_params.mat');
%ind_t13=find(filedata.params(:,1)/10000>=tlin(1)&filedata.params(:,1)/10000<=tlin(end)&filedata.params(:,8)>10&(filedata.params(:,2)>100|filedata.params(:,3)>100|filedata.params(:,4)>100|filedata.params(:,5)>100));
ind_t13=find(filedata.params(:,1)/10000>=tlin(1)&filedata.params(:,1)/10000<=tlin(end));
time_t13=filedata.params(ind_t13,1);
mark0_t13=[filedata.params(ind_t13,2) filedata.params(ind_t13,3) filedata.params(ind_t13,4) filedata.params(ind_t13,5)];
time2_t13=time_t13/10000;
spikeT0_t13=time2_t13;
[procInd0_t13,procInd1_t13]=histc(spikeT0_t13,tlin);
procInd_t13=find(procInd0_t13);
spikeT_t13=tlin(procInd_t13);
spike_t13=procInd0_t13';
[~,rawInd0_t13]=histc(spikeT0_t13,time2_t13);
markAll_t13(:,1)=procInd1_t13;markAll_t13(:,2:5)=mark0_t13(rawInd0_t13(rawInd0_t13~=0),:);
mdel=20; ms=min(min(markAll_t13(:,2:5))):mdel:max(max(markAll_t13(:,2:5)));
sxker=2*xdel; smker = mdel; T=size(tlin,1);
occ=normpdf(xs'*ones(1,T),ones(length(xs),1)*xtrain,sxker)*ones(T,length(ms));
%occ: columns are identical; occupancy based on position; denominator
Xnum_t13=normpdf(xs'*ones(1,length(spikeT0_t13)),ones(length(xs),1)*xtrain(procInd1_t13),sxker);
%Xnum: Gaussian kernel estimators for position
Lint_t13=sum(Xnum_t13,2)./occ(:,1)./dt; %integral
Lint_t13=Lint_t13./sum(Lint_t13);

load('bond_data\bond04\14-096\bond04-14_params.mat');
%ind_t14=find(filedata.params(:,1)/10000>=tlin(1)&filedata.params(:,1)/10000<=tlin(end)&filedata.params(:,8)>10&(filedata.params(:,2)>100|filedata.params(:,3)>100|filedata.params(:,4)>100|filedata.params(:,5)>100));
ind_t14=find(filedata.params(:,1)/10000>=tlin(1)&filedata.params(:,1)/10000<=tlin(end));
time_t14=filedata.params(ind_t14,1);
mark0_t14=[filedata.params(ind_t14,2) filedata.params(ind_t14,3) filedata.params(ind_t14,4) filedata.params(ind_t14,5)];
time2_t14=time_t14/10000;
spikeT0_t14=time2_t14;
[procInd0_t14,procInd1_t14]=histc(spikeT0_t14,tlin);
procInd_t14=find(procInd0_t14);
spikeT_t14=tlin(procInd_t14);
spike_t14=procInd0_t14';
[~,rawInd0_t14]=histc(spikeT0_t14,time2_t14);
markAll_t14(:,1)=procInd1_t14;markAll_t14(:,2:5)=mark0_t14(rawInd0_t14(rawInd0_t14~=0),:);
mdel=20; ms=min(min(markAll_t14(:,2:5))):mdel:max(max(markAll_t14(:,2:5)));
sxker=2*xdel; smker = mdel; T=size(tlin,1);
occ=normpdf(xs'*ones(1,T),ones(length(xs),1)*xtrain,sxker)*ones(T,length(ms));
%occ: columns are identical; occupancy based on position; denominator
Xnum_t14=normpdf(xs'*ones(1,length(spikeT0_t14)),ones(length(xs),1)*xtrain(procInd1_t14),sxker);
%Xnum: Gaussian kernel estimators for position
Lint_t14=sum(Xnum_t14,2)./occ(:,1)./dt; %integral
Lint_t14=Lint_t14./sum(Lint_t14);

load('bond_data\bond04\17-108\bond04-17_params.mat');
%ind_t17=find(filedata.params(:,1)/10000>=tlin(1)&filedata.params(:,1)/10000<=tlin(end)&filedata.params(:,8)>10&(filedata.params(:,2)>100|filedata.params(:,3)>100|filedata.params(:,4)>100|filedata.params(:,5)>100));
ind_t17=find(filedata.params(:,1)/10000>=tlin(1)&filedata.params(:,1)/10000<=tlin(end));
time_t17=filedata.params(ind_t17,1);
mark0_t17=[filedata.params(ind_t17,2) filedata.params(ind_t17,3) filedata.params(ind_t17,4) filedata.params(ind_t17,5)];
time2_t17=time_t17/10000;
spikeT0_t17=time2_t17;
[procInd0_t17,procInd1_t17]=histc(spikeT0_t17,tlin);
procInd_t17=find(procInd0_t17);
spikeT_t17=tlin(procInd_t17);
spike_t17=procInd0_t17';
[~,rawInd0_t17]=histc(spikeT0_t17,time2_t17);
markAll_t17(:,1)=procInd1_t17;markAll_t17(:,2:5)=mark0_t17(rawInd0_t17(rawInd0_t17~=0),:);
mdel=20; ms=min(min(markAll_t17(:,2:5))):mdel:max(max(markAll_t17(:,2:5)));
sxker=2*xdel; smker = mdel; T=size(tlin,1);
occ=normpdf(xs'*ones(1,T),ones(length(xs),1)*xtrain,sxker)*ones(T,length(ms));
%occ: columns are identical; occupancy based on position; denominator
Xnum_t17=normpdf(xs'*ones(1,length(spikeT0_t17)),ones(length(xs),1)*xtrain(procInd1_t17),sxker);
%Xnum: Gaussian kernel estimators for position
Lint_t17=sum(Xnum_t17,2)./occ(:,1)./dt; %integral
Lint_t17=Lint_t17./sum(Lint_t17);

load('bond_data\bond04\18-138\bond04-18_params.mat');
%ind_t18=find(filedata.params(:,1)/10000>=tlin(1)&filedata.params(:,1)/10000<=tlin(end)&filedata.params(:,8)>10&(filedata.params(:,2)>100|filedata.params(:,3)>100|filedata.params(:,4)>100|filedata.params(:,5)>100));
ind_t18=find(filedata.params(:,1)/10000>=tlin(1)&filedata.params(:,1)/10000<=tlin(end));
time_t18=filedata.params(ind_t18,1);
mark0_t18=[filedata.params(ind_t18,2) filedata.params(ind_t18,3) filedata.params(ind_t18,4) filedata.params(ind_t18,5)];
time2_t18=time_t18/10000;
spikeT0_t18=time2_t18;
[procInd0_t18,procInd1_t18]=histc(spikeT0_t18,tlin);
procInd_t18=find(procInd0_t18);
spikeT_t18=tlin(procInd_t18);
spike_t18=procInd0_t18';
[~,rawInd0_t18]=histc(spikeT0_t18,time2_t18);
markAll_t18(:,1)=procInd1_t18;markAll_t18(:,2:5)=mark0_t18(rawInd0_t18(rawInd0_t18~=0),:);
mdel=20; ms=min(min(markAll_t18(:,2:5))):mdel:max(max(markAll_t18(:,2:5)));
sxker=2*xdel; smker = mdel; T=size(tlin,1);
occ=normpdf(xs'*ones(1,T),ones(length(xs),1)*xtrain,sxker)*ones(T,length(ms));
%occ: columns are identical; occupancy based on position; denominator
Xnum_t18=normpdf(xs'*ones(1,length(spikeT0_t18)),ones(length(xs),1)*xtrain(procInd1_t18),sxker);
%Xnum: Gaussian kernel estimators for position
Lint_t18=sum(Xnum_t18,2)./occ(:,1)./dt; %integral
Lint_t18=Lint_t18./sum(Lint_t18);

load('bond_data\bond04\19-130\bond04-19_params.mat');
%ind_t19=find(filedata.params(:,1)/10000>=tlin(1)&filedata.params(:,1)/10000<=tlin(end)&filedata.params(:,8)>10&(filedata.params(:,2)>100|filedata.params(:,3)>100|filedata.params(:,4)>100|filedata.params(:,5)>100));
ind_t19=find(filedata.params(:,1)/10000>=tlin(1)&filedata.params(:,1)/10000<=tlin(end));
time_t19=filedata.params(ind_t19,1);
mark0_t19=[filedata.params(ind_t19,2) filedata.params(ind_t19,3) filedata.params(ind_t19,4) filedata.params(ind_t19,5)];
time2_t19=time_t19/10000;
spikeT0_t19=time2_t19;
[procInd0_t19,procInd1_t19]=histc(spikeT0_t19,tlin);
procInd_t19=find(procInd0_t19);
spikeT_t19=tlin(procInd_t19);
spike_t19=procInd0_t19';
[~,rawInd0_t19]=histc(spikeT0_t19,time2_t19);
markAll_t19(:,1)=procInd1_t19;markAll_t19(:,2:5)=mark0_t19(rawInd0_t19(rawInd0_t19~=0),:);
mdel=20; ms=min(min(markAll_t19(:,2:5))):mdel:max(max(markAll_t19(:,2:5)));
sxker=2*xdel; smker = mdel; T=size(tlin,1);
occ=normpdf(xs'*ones(1,T),ones(length(xs),1)*xtrain,sxker)*ones(T,length(ms));
%occ: columns are identical; occupancy based on position; denominator
Xnum_t19=normpdf(xs'*ones(1,length(spikeT0_t19)),ones(length(xs),1)*xtrain(procInd1_t19),sxker);
%Xnum: Gaussian kernel estimators for position
Lint_t19=sum(Xnum_t19,2)./occ(:,1)./dt; %integral
Lint_t19=Lint_t19./sum(Lint_t19);

load('bond_data\bond04\20-132\bond04-20_params.mat');
%ind_t20=find(filedata.params(:,1)/10000>=tlin(1)&filedata.params(:,1)/10000<=tlin(end)&filedata.params(:,8)>10&(filedata.params(:,2)>100|filedata.params(:,3)>100|filedata.params(:,4)>100|filedata.params(:,5)>100));
ind_t20=find(filedata.params(:,1)/10000>=tlin(1)&filedata.params(:,1)/10000<=tlin(end));
time_t20=filedata.params(ind_t20,1);
mark0_t20=[filedata.params(ind_t20,2) filedata.params(ind_t20,3) filedata.params(ind_t20,4) filedata.params(ind_t20,5)];
time2_t20=time_t20/10000;
spikeT0_t20=time2_t20;
[procInd0_t20,procInd1_t20]=histc(spikeT0_t20,tlin);
procInd_t20=find(procInd0_t20);
spikeT_t20=tlin(procInd_t20);
spike_t20=procInd0_t20';
[~,rawInd0_t20]=histc(spikeT0_t20,time2_t20);
markAll_t20(:,1)=procInd1_t20;markAll_t20(:,2:5)=mark0_t20(rawInd0_t20(rawInd0_t20~=0),:);
mdel=20; ms=min(min(markAll_t20(:,2:5))):mdel:max(max(markAll_t20(:,2:5)));
sxker=2*xdel; smker = mdel; T=size(tlin,1);
occ=normpdf(xs'*ones(1,T),ones(length(xs),1)*xtrain,sxker)*ones(T,length(ms));
%occ: columns are identical; occupancy based on position; denominator
Xnum_t20=normpdf(xs'*ones(1,length(spikeT0_t20)),ones(length(xs),1)*xtrain(procInd1_t20),sxker);
%Xnum: Gaussian kernel estimators for position
Lint_t20=sum(Xnum_t20,2)./occ(:,1)./dt; %integral
Lint_t20=Lint_t20./sum(Lint_t20);

load('bond_data\bond04\22-171\bond04-22_params.mat');
%ind_t22=find(filedata.params(:,1)/10000>=tlin(1)&filedata.params(:,1)/10000<=tlin(end)&filedata.params(:,8)>10&(filedata.params(:,2)>100|filedata.params(:,3)>100|filedata.params(:,4)>100|filedata.params(:,5)>100));
ind_t22=find(filedata.params(:,1)/10000>=tlin(1)&filedata.params(:,1)/10000<=tlin(end));
time_t22=filedata.params(ind_t22,1);
mark0_t22=[filedata.params(ind_t22,2) filedata.params(ind_t22,3) filedata.params(ind_t22,4) filedata.params(ind_t22,5)];
time2_t22=time_t22/10000;
spikeT0_t22=time2_t22;
[procInd0_t22,procInd1_t22]=histc(spikeT0_t22,tlin);
procInd_t22=find(procInd0_t22);
spikeT_t22=tlin(procInd_t22);
spike_t22=procInd0_t22';
[~,rawInd0_t22]=histc(spikeT0_t22,time2_t22);
markAll_t22(:,1)=procInd1_t22;markAll_t22(:,2:5)=mark0_t22(rawInd0_t22(rawInd0_t22~=0),:);
mdel=20; ms=min(min(markAll_t22(:,2:5))):mdel:max(max(markAll_t22(:,2:5)));
sxker=2*xdel; smker = mdel; T=size(tlin,1);
occ=normpdf(xs'*ones(1,T),ones(length(xs),1)*xtrain,sxker)*ones(T,length(ms));
%occ: columns are identical; occupancy based on position; denominator
Xnum_t22=normpdf(xs'*ones(1,length(spikeT0_t22)),ones(length(xs),1)*xtrain(procInd1_t22),sxker);
%Xnum: Gaussian kernel estimators for position
Lint_t22=sum(Xnum_t22,2)./occ(:,1)./dt; %integral
Lint_t22=Lint_t22./sum(Lint_t22);

load('bond_data\bond04\23-152\bond04-23_params.mat');
%ind_t23=find(filedata.params(:,1)/10000>=tlin(1)&filedata.params(:,1)/10000<=tlin(end)&filedata.params(:,8)>10&(filedata.params(:,2)>100|filedata.params(:,3)>100|filedata.params(:,4)>100|filedata.params(:,5)>100));
ind_t23=find(filedata.params(:,1)/10000>=tlin(1)&filedata.params(:,1)/10000<=tlin(end));
time_t23=filedata.params(ind_t23,1);
mark0_t23=[filedata.params(ind_t23,2) filedata.params(ind_t23,3) filedata.params(ind_t23,4) filedata.params(ind_t23,5)];
time2_t23=time_t23/10000;
spikeT0_t23=time2_t23;
[procInd0_t23,procInd1_t23]=histc(spikeT0_t23,tlin);
procInd_t23=find(procInd0_t23);
spikeT_t23=tlin(procInd_t23);
spike_t23=procInd0_t23';
[~,rawInd0_t23]=histc(spikeT0_t23,time2_t23);
markAll_t23(:,1)=procInd1_t23;markAll_t23(:,2:5)=mark0_t23(rawInd0_t23(rawInd0_t23~=0),:);
mdel=20; ms=min(min(markAll_t23(:,2:5))):mdel:max(max(markAll_t23(:,2:5)));
sxker=2*xdel; smker = mdel; T=size(tlin,1);
occ=normpdf(xs'*ones(1,T),ones(length(xs),1)*xtrain,sxker)*ones(T,length(ms));
%occ: columns are identical; occupancy based on position; denominator
Xnum_t23=normpdf(xs'*ones(1,length(spikeT0_t23)),ones(length(xs),1)*xtrain(procInd1_t23),sxker);
%Xnum: Gaussian kernel estimators for position
Lint_t23=sum(Xnum_t23,2)./occ(:,1)./dt; %integral
Lint_t23=Lint_t23./sum(Lint_t23);

load('bond_data\bond04\27-158\bond04-27_params.mat');
%ind_t27=find(filedata.params(:,1)/10000>=tlin(1)&filedata.params(:,1)/10000<=tlin(end)&filedata.params(:,8)>10&(filedata.params(:,2)>100|filedata.params(:,3)>100|filedata.params(:,4)>100|filedata.params(:,5)>100));
ind_t27=find(filedata.params(:,1)/10000>=tlin(1)&filedata.params(:,1)/10000<=tlin(end));
time_t27=filedata.params(ind_t27,1);
mark0_t27=[filedata.params(ind_t27,2) filedata.params(ind_t27,3) filedata.params(ind_t27,4) filedata.params(ind_t27,5)];
time2_t27=time_t27/10000;
spikeT0_t27=time2_t27;
[procInd0_t27,procInd1_t27]=histc(spikeT0_t27,tlin);
procInd_t27=find(procInd0_t27);
spikeT_t27=tlin(procInd_t27);
spike_t27=procInd0_t27';
[~,rawInd0_t27]=histc(spikeT0_t27,time2_t27);
markAll_t27(:,1)=procInd1_t27;markAll_t27(:,2:5)=mark0_t27(rawInd0_t27(rawInd0_t27~=0),:);
mdel=20; ms=min(min(markAll_t27(:,2:5))):mdel:max(max(markAll_t27(:,2:5)));
sxker=2*xdel; smker = mdel; T=size(tlin,1);
occ=normpdf(xs'*ones(1,T),ones(length(xs),1)*xtrain,sxker)*ones(T,length(ms));
%occ: columns are identical; occupancy based on position; denominator
Xnum_t27=normpdf(xs'*ones(1,length(spikeT0_t27)),ones(length(xs),1)*xtrain(procInd1_t27),sxker);
%Xnum: Gaussian kernel estimators for position
Lint_t27=sum(Xnum_t27,2)./occ(:,1)./dt; %integral
Lint_t27=Lint_t27./sum(Lint_t27);

load('bond_data\bond04\29-118\bond04-29_params.mat');
%ind_t29=find(filedata.params(:,1)/10000>=tlin(1)&filedata.params(:,1)/10000<=tlin(end)&filedata.params(:,8)>10&(filedata.params(:,2)>100|filedata.params(:,3)>100|filedata.params(:,4)>100|filedata.params(:,5)>100));
ind_t29=find(filedata.params(:,1)/10000>=tlin(1)&filedata.params(:,1)/10000<=tlin(end));
time_t29=filedata.params(ind_t29,1);
mark0_t29=[filedata.params(ind_t29,2) filedata.params(ind_t29,3) filedata.params(ind_t29,4) filedata.params(ind_t29,5)];
time2_t29=time_t29/10000;
spikeT0_t29=time2_t29;
[procInd0_t29,procInd1_t29]=histc(spikeT0_t29,tlin);
procInd_t29=find(procInd0_t29);
spikeT_t29=tlin(procInd_t29);
spike_t29=procInd0_t29';
[~,rawInd0_t29]=histc(spikeT0_t29,time2_t29);
markAll_t29(:,1)=procInd1_t29;markAll_t29(:,2:5)=mark0_t29(rawInd0_t29(rawInd0_t29~=0),:);
mdel=20; ms=min(min(markAll_t29(:,2:5))):mdel:max(max(markAll_t29(:,2:5)));
sxker=2*xdel; smker = mdel; T=size(tlin,1);
occ=normpdf(xs'*ones(1,T),ones(length(xs),1)*xtrain,sxker)*ones(T,length(ms));
%occ: columns are identical; occupancy based on position; denominator
Xnum_t29=normpdf(xs'*ones(1,length(spikeT0_t29)),ones(length(xs),1)*xtrain(procInd1_t29),sxker);
%Xnum: Gaussian kernel estimators for position
Lint_t29=sum(Xnum_t29,2)./occ(:,1)./dt; %integral
Lint_t29=Lint_t29./sum(Lint_t29);



%%
time0=[time_t1;time_t2;time_t4;time_t5;time_t7;time_t10;time_t11;time_t12;time_t13;time_t14;time_t17;time_t18;time_t19;time_t20;time_t22;time_t23;time_t27;time_t29];
[time,timeInd]=sort(time0);
mark0=[mark0_t1;mark0_t2;mark0_t4;mark0_t5;mark0_t7;mark0_t10;mark0_t11;mark0_t12;mark0_t13;mark0_t14;mark0_t17;mark0_t18;mark0_t19;mark0_t20;mark0_t22;mark0_t23;mark0_t27;mark0_t29];
mark0=mark0(timeInd,:);
procInd1=[procInd1_t1;procInd1_t2;procInd1_t4;procInd1_t5;procInd1_t7;procInd1_t10;procInd1_t11;procInd1_t12;procInd1_t13;procInd1_t14;procInd1_t17;procInd1_t18;procInd1_t19;procInd1_t20;procInd1_t22;procInd1_t23;procInd1_t27;procInd1_t29];
procInd1=procInd1(timeInd,:);

len_t1=length(time_t1);len_t2=length(time_t2);len_t4=length(time_t4);len_t5=length(time_t5);
len_t7=length(time_t7);len_t11=length(time_t11);len_t10=length(time_t10);len_t12=length(time_t12);
len_t13=length(time_t13);len_t14=length(time_t14);len_t17=length(time_t17);len_t18=length(time_t18);
len_t19=length(time_t19);len_t20=length(time_t20);len_t22=length(time_t22);len_t23=length(time_t23);
len_t27=length(time_t27);len_t29=length(time_t29);

tet_ind=zeros(length(time),5);
%indicator matrix
%row: time point; column: which tetrode spikes
for i=1:length(time0)
    if timeInd(i)>=1 && timeInd(i)<=len_t1
        tet_ind(i,1)=1; %tet1
    elseif timeInd(i)>=len_t1+1 && timeInd(i)<=len_t1+len_t2
        tet_ind(i,2)=1; %tet2
    elseif timeInd(i)>=len_t1+len_t2+1 && timeInd(i)<=len_t1+len_t2+len_t4
        tet_ind(i,3)=1; %tet4
    elseif timeInd(i)>=len_t1+len_t2+len_t4+1 && timeInd(i)<=len_t1+len_t2+len_t4+len_t5
        tet_ind(i,4)=1; %tet5
    elseif timeInd(i)>=len_t1+len_t2+len_t4+len_t5+1 && timeInd(i)<=len_t1+len_t2+len_t4+len_t5+len_t7
        tet_ind(i,5)=1; %tet7
    elseif timeInd(i)>=len_t1+len_t2+len_t4+len_t5+len_t7+1 && timeInd(i)<=len_t1+len_t2+len_t4+len_t5+len_t7+len_t10
        tet_ind(i,6)=1; %tet10
    elseif timeInd(i)>=len_t1+len_t2+len_t4+len_t5+len_t7+len_t10+1 && timeInd(i)<=len_t1+len_t2+len_t4+len_t5+len_t7+len_t10+len_t11
        tet_ind(i,7)=1; %tet11
    elseif timeInd(i)>=len_t1+len_t2+len_t4+len_t5+len_t7+len_t10+len_t11+1 && timeInd(i)<=len_t1+len_t2+len_t4+len_t5+len_t7+len_t10+len_t11+len_t12
        tet_ind(i,8)=1; %tet12
    elseif timeInd(i)>=len_t1+len_t2+len_t4+len_t5+len_t7+len_t10+len_t11+len_t12+1 && timeInd(i)<=len_t1+len_t2+len_t4+len_t5+len_t7+len_t10+len_t11+len_t12+len_t13
        tet_ind(i,9)=1; %tet13
    elseif timeInd(i)>=len_t1+len_t2+len_t4+len_t5+len_t7+len_t10+len_t11+len_t12+len_t13+1 && timeInd(i)<=len_t1+len_t2+len_t4+len_t5+len_t7+len_t10+len_t11+len_t12+len_t13+len_t14
        tet_ind(i,10)=1; %tet14
    elseif timeInd(i)>=len_t1+len_t2+len_t4+len_t5+len_t7+len_t10+len_t11+len_t12+len_t13+len_t14+1 && timeInd(i)<=len_t1+len_t2+len_t4+len_t5+len_t7+len_t10+len_t11+len_t12+len_t13+len_t14+len_t17
        tet_ind(i,11)=1; %tet17
    elseif timeInd(i)>=len_t1+len_t2+len_t4+len_t5+len_t7+len_t10+len_t11+len_t12+len_t13+len_t14+len_t17+1 && timeInd(i)<=len_t1+len_t2+len_t4+len_t5+len_t7+len_t10+len_t11+len_t12+len_t13+len_t14+len_t17+len_t18
        tet_ind(i,12)=1; %tet18
    elseif timeInd(i)>=len_t1+len_t2+len_t4+len_t5+len_t7+len_t10+len_t11+len_t12+len_t13+len_t14+len_t17+len_t18+1 && timeInd(i)<=len_t1+len_t2+len_t4+len_t5+len_t7+len_t10+len_t11+len_t12+len_t13+len_t14+len_t17+len_t18+len_t19
        tet_ind(i,13)=1; %tet19
    elseif timeInd(i)>=len_t1+len_t2+len_t4+len_t5+len_t7+len_t10+len_t11+len_t12+len_t13+len_t14+len_t17+len_t18+len_t19+1 && timeInd(i)<=len_t1+len_t2+len_t4+len_t5+len_t7+len_t10+len_t11+len_t12+len_t13+len_t14+len_t17+len_t18+len_t19+len_t20
        tet_ind(i,14)=1; %tet20
    elseif timeInd(i)>=len_t1+len_t2+len_t4+len_t5+len_t7+len_t10+len_t11+len_t12+len_t13+len_t14+len_t17+len_t18+len_t19+len_t20+1 && timeInd(i)<=len_t1+len_t2+len_t4+len_t5+len_t7+len_t10+len_t11+len_t12+len_t13+len_t14+len_t17+len_t18+len_t19+len_t20+len_t22
        tet_ind(i,15)=1; %tet22
    elseif timeInd(i)>=len_t1+len_t2+len_t4+len_t5+len_t7+len_t10+len_t11+len_t12+len_t13+len_t14+len_t17+len_t18+len_t19+len_t20+len_t22+1 && timeInd(i)<=len_t1+len_t2+len_t4+len_t5+len_t7+len_t10+len_t11+len_t12+len_t13+len_t14+len_t17+len_t18+len_t19+len_t20+len_t22+len_t23
        tet_ind(i,16)=1; %tet23
    elseif timeInd(i)>=len_t1+len_t2+len_t4+len_t5+len_t7+len_t10+len_t11+len_t12+len_t13+len_t14+len_t17+len_t18+len_t19+len_t20+len_t22+len_t23+1 && timeInd(i)<=len_t1+len_t2+len_t4+len_t5+len_t7+len_t10+len_t11+len_t12+len_t13+len_t14+len_t17+len_t18+len_t19+len_t20+len_t22+len_t23+len_t27
        tet_ind(i,17)=1; %tet27
    elseif timeInd(i)>=len_t1+len_t2+len_t4+len_t5+len_t7+len_t10+len_t11+len_t12+len_t13+len_t14+len_t17+len_t18+len_t19+len_t20+len_t22+len_t23+len_t27+1 && timeInd(i)<=len_t1+len_t2+len_t4+len_t5+len_t7+len_t10+len_t11+len_t12+len_t13+len_t14+len_t17+len_t18+len_t19+len_t20+len_t22+len_t23+len_t27+len_t29
        tet_ind(i,18)=1; %tet29
    end
end

tet_sum=tet_ind.*cumsum(tet_ind,1); %row: time point; column: index of spike per tetrode


%%
mdel=20; ms=100:mdel:max(mark0(:)); 
sxker=2*xdel; smker = mdel; T=size(tlin,1);
occ=normpdf(xs'*ones(1,T),ones(length(xs),1)*xtrain,sxker)*ones(T,length(ms));
%occ: columns are identical; occupancy based on position; denominator
Xnum=normpdf(xs'*ones(1,length(time0)),ones(length(xs),1)*xtrain(procInd1),sxker);
%Xnum: Gaussian kernel estimators for position
Lint=sum(Xnum,2)./occ(:,1)./dt; %integral
Lint=Lint./sum(Lint);
%Lint: conditional intensity function for the unmarked case


%%




%%
%%%%binary
procInd1_I1=procInd1(ismember(procInd1,ind_I1));
Xnum_I1=normpdf(xs'*ones(1,length(xtrain(procInd1_I1))),ones(length(xs),1)*xtrain(procInd1_I1),sxker);
%Xnum: Gaussian kernel estimators for position
Lint_I1=sum(Xnum_I1,2)./occ(:,1)./dt; %integral
Lint_I1=Lint_I1./sum(Lint_I1);


procInd1_I0=procInd1(ismember(procInd1,ind_I0));
Xnum_I0=normpdf(xs'*ones(1,length(xtrain(procInd1_I0))),ones(length(xs),1)*xtrain(procInd1_I0),sxker);
%Xnum: Gaussian kernel estimators for position
Lint_I0=sum(Xnum_I0,2)./occ(:,1)./dt; %integral
Lint_I0=Lint_I0./sum(Lint_I0);




%%
%%%%encode per tetrode, binary
load('bond_data\bond04\01-149\bond04-01_params.mat');
%ind_t1=find(filedata.params(:,1)/10000>=tlin(1)&filedata.params(:,1)/10000<=tlin(end)&filedata.params(:,8)>10&(filedata.params(:,2)>100|filedata.params(:,3)>100|filedata.params(:,4)>100|filedata.params(:,5)>100));
%threshold at 100uV
ind_t1=find(filedata.params(:,1)/10000>=tlin(1)&filedata.params(:,1)/10000<=tlin(end));
time_t1=filedata.params(ind_t1,1);
mark0_t1=[filedata.params(ind_t1,2) filedata.params(ind_t1,3) filedata.params(ind_t1,4) filedata.params(ind_t1,5)];
time2_t1=time_t1/10000;
spikeT0_t1=time2_t1;
[procInd0_t1,procInd1_t1]=histc(spikeT0_t1,tlin);
procInd_t1=find(procInd0_t1);
spikeT_t1=tlin(procInd_t1);
spike_t1=procInd0_t1';
[~,rawInd0_t1]=histc(spikeT0_t1,time2_t1);
markAll_t1(:,1)=procInd1_t1;markAll_t1(:,2:5)=mark0_t1(rawInd0_t1(rawInd0_t1~=0),:);
mdel=20; ms=min(min(markAll_t1(:,2:5))):mdel:max(max(markAll_t1(:,2:5))); 
sxker=2*xdel; smker = mdel; T=size(tlin,1);
occ=normpdf(xs'*ones(1,T),ones(length(xs),1)*xtrain,sxker)*ones(T,length(ms));
%occ: columns are identical; occupancy based on position; denominator
procInd1_t1_I1a=procInd1_t1(ismember(procInd1_t1,ind_I1));
procInd1_t1_I1=find(ismember(procInd1_t1,ind_I1));
Xnum_t1_I1=normpdf(xs'*ones(1,length(procInd1_t1_I1a)),ones(length(xs),1)*xtrain(procInd1_t1_I1a),sxker);
%Xnum: Gaussian kernel estimators for position
Lint_t1_I1=sum(Xnum_t1_I1,2)./occ(:,1)./dt; %integral
Lint_t1_I1=Lint_t1_I1./sum(Lint_t1_I1);
procInd1_t1_I0a=procInd1_t1(ismember(procInd1_t1,ind_I0));
procInd1_t1_I0=find(ismember(procInd1_t1,ind_I0));
Xnum_t1_I0=normpdf(xs'*ones(1,length(procInd1_t1_I0a)),ones(length(xs),1)*xtrain(procInd1_t1_I0a),sxker);
%Xnum: Gaussian kernel estimators for position
Lint_t1_I0=sum(Xnum_t1_I0,2)./occ(:,1)./dt; %integral
Lint_t1_I0=Lint_t1_I0./sum(Lint_t1_I0);

load('bond_data\bond04\02-151\bond04-02_params.mat');
%ind_t2=find(filedata.params(:,1)/10000>=tlin(1)&filedata.params(:,1)/10000<=tlin(end)&filedata.params(:,8)>10&(filedata.params(:,2)>100|filedata.params(:,3)>100|filedata.params(:,4)>100|filedata.params(:,5)>100));
ind_t2=find(filedata.params(:,1)/10000>=tlin(1)&filedata.params(:,1)/10000<=tlin(end));
time_t2=filedata.params(ind_t2,1);
mark0_t2=[filedata.params(ind_t2,2) filedata.params(ind_t2,3) filedata.params(ind_t2,4) filedata.params(ind_t2,5)];
time2_t2=time_t2/10000;
spikeT0_t2=time2_t2;
[procInd0_t2,procInd1_t2]=histc(spikeT0_t2,tlin);
procInd_t2=find(procInd0_t2);
spikeT_t2=tlin(procInd_t2);
spike_t2=procInd0_t2';
[~,rawInd0_t2]=histc(spikeT0_t2,time2_t2);
markAll_t2(:,1)=procInd1_t2;markAll_t2(:,2:5)=mark0_t2(rawInd0_t2(rawInd0_t2~=0),:);
mdel=20; ms=min(min(markAll_t2(:,2:5))):mdel:max(max(markAll_t2(:,2:5))); 
sxker=2*xdel; smker = mdel; T=size(tlin,1);
occ=normpdf(xs'*ones(1,T),ones(length(xs),1)*xtrain,sxker)*ones(T,length(ms));
%occ: columns are identical; occupancy based on position; denominator
procInd1_t2_I1a=procInd1_t2(ismember(procInd1_t2,ind_I1));
procInd1_t2_I1=find(ismember(procInd1_t2,ind_I1));
Xnum_t2_I1=normpdf(xs'*ones(1,length(procInd1_t2_I1a)),ones(length(xs),1)*xtrain(procInd1_t2_I1a),sxker);
%Xnum: Gaussian kernel estimators for position
Lint_t2_I1=sum(Xnum_t2_I1,2)./occ(:,1)./dt; %integral
Lint_t2_I1=Lint_t2_I1./sum(Lint_t2_I1);
procInd1_t2_I0a=procInd1_t2(ismember(procInd1_t2,ind_I0));
procInd1_t2_I0=find(ismember(procInd1_t2,ind_I0));
Xnum_t2_I0=normpdf(xs'*ones(1,length(procInd1_t2_I0a)),ones(length(xs),1)*xtrain(procInd1_t2_I0a),sxker);
%Xnum: Gaussian kernel estimators for position
Lint_t2_I0=sum(Xnum_t2_I0,2)./occ(:,1)./dt; %integral
Lint_t2_I0=Lint_t2_I0./sum(Lint_t2_I0);

load('bond_data\bond04\04-107\bond04-04_params.mat');
%ind_t4=find(filedata.params(:,1)/10000>=tlin(1)&filedata.params(:,1)/10000<=tlin(end)&filedata.params(:,8)>10&(filedata.params(:,2)>100|filedata.params(:,3)>100|filedata.params(:,4)>100|filedata.params(:,5)>100));
ind_t4=find(filedata.params(:,1)/10000>=tlin(1)&filedata.params(:,1)/10000<=tlin(end));
time_t4=filedata.params(ind_t4,1);
mark0_t4=[filedata.params(ind_t4,2) filedata.params(ind_t4,3) filedata.params(ind_t4,4) filedata.params(ind_t4,5)];
time2_t4=time_t4/10000;
spikeT0_t4=time2_t4;
[procInd0_t4,procInd1_t4]=histc(spikeT0_t4,tlin);
procInd_t4=find(procInd0_t4);
spikeT_t4=tlin(procInd_t4);
spike_t4=procInd0_t4';
[~,rawInd0_t4]=histc(spikeT0_t4,time2_t4);
markAll_t4(:,1)=procInd1_t4;markAll_t4(:,2:5)=mark0_t4(rawInd0_t4(rawInd0_t4~=0),:);
mdel=20; ms=min(min(markAll_t4(:,2:5))):mdel:max(max(markAll_t4(:,2:5))); 
sxker=2*xdel; smker = mdel; T=size(tlin,1);
occ=normpdf(xs'*ones(1,T),ones(length(xs),1)*xtrain,sxker)*ones(T,length(ms));
%occ: columns are identical; occupancy based on position; denominator
procInd1_t4_I1a=procInd1_t4(ismember(procInd1_t4,ind_I1));
procInd1_t4_I1=find(ismember(procInd1_t4,ind_I1));
Xnum_t4_I1=normpdf(xs'*ones(1,length(procInd1_t4_I1a)),ones(length(xs),1)*xtrain(procInd1_t4_I1a),sxker);
%Xnum: Gaussian kernel estimators for position
Lint_t4_I1=sum(Xnum_t4_I1,2)./occ(:,1)./dt; %integral
Lint_t4_I1=Lint_t4_I1./sum(Lint_t4_I1);
procInd1_t4_I0a=procInd1_t4(ismember(procInd1_t4,ind_I0));
procInd1_t4_I0=find(ismember(procInd1_t4,ind_I0));
Xnum_t4_I0=normpdf(xs'*ones(1,length(procInd1_t4_I0a)),ones(length(xs),1)*xtrain(procInd1_t4_I0a),sxker);
%Xnum: Gaussian kernel estimators for position
Lint_t4_I0=sum(Xnum_t4_I0,2)./occ(:,1)./dt; %integral
Lint_t4_I0=Lint_t4_I0./sum(Lint_t4_I0);

load('bond_data\bond04\05-105\bond04-05_params.mat');
%ind_t5=find(filedata.params(:,1)/10000>=tlin(1)&filedata.params(:,1)/10000<=tlin(end)&filedata.params(:,8)>10&(filedata.params(:,2)>100|filedata.params(:,3)>100|filedata.params(:,4)>100|filedata.params(:,5)>100));
ind_t5=find(filedata.params(:,1)/10000>=tlin(1)&filedata.params(:,1)/10000<=tlin(end));
time_t5=filedata.params(ind_t5,1);
mark0_t5=[filedata.params(ind_t5,2) filedata.params(ind_t5,3) filedata.params(ind_t5,4) filedata.params(ind_t5,5)];
time2_t5=time_t5/10000;
spikeT0_t5=time2_t5;
[procInd0_t5,procInd1_t5]=histc(spikeT0_t5,tlin);
procInd_t5=find(procInd0_t5);
spikeT_t5=tlin(procInd_t5);
spike_t5=procInd0_t5';
[~,rawInd0_t5]=histc(spikeT0_t5,time2_t5);
markAll_t5(:,1)=procInd1_t5;markAll_t5(:,2:5)=mark0_t5(rawInd0_t5(rawInd0_t5~=0),:);
mdel=20; ms=min(min(markAll_t5(:,2:5))):mdel:max(max(markAll_t5(:,2:5))); 
sxker=2*xdel; smker = mdel; T=size(tlin,1);
occ=normpdf(xs'*ones(1,T),ones(length(xs),1)*xtrain,sxker)*ones(T,length(ms));
%occ: columns are identical; occupancy based on position; denominator
procInd1_t5_I1a=procInd1_t5(ismember(procInd1_t5,ind_I1));
procInd1_t5_I1=find(ismember(procInd1_t5,ind_I1));
Xnum_t5_I1=normpdf(xs'*ones(1,length(procInd1_t5_I1a)),ones(length(xs),1)*xtrain(procInd1_t5_I1a),sxker);
%Xnum: Gaussian kernel estimators for position
Lint_t5_I1=sum(Xnum_t5_I1,2)./occ(:,1)./dt; %integral
Lint_t5_I1=Lint_t5_I1./sum(Lint_t5_I1);
procInd1_t5_I0a=procInd1_t5(ismember(procInd1_t5,ind_I0));
procInd1_t5_I0=find(ismember(procInd1_t5,ind_I0));
Xnum_t5_I0=normpdf(xs'*ones(1,length(procInd1_t5_I0a)),ones(length(xs),1)*xtrain(procInd1_t5_I0a),sxker);
%Xnum: Gaussian kernel estimators for position
Lint_t5_I0=sum(Xnum_t5_I0,2)./occ(:,1)./dt; %integral
Lint_t5_I0=Lint_t5_I0./sum(Lint_t5_I0);

load('bond_data\bond04\07-162\bond04-07_params.mat');
%ind_t7=find(filedata.params(:,1)/10000>=tlin(1)&filedata.params(:,1)/10000<=tlin(end)&filedata.params(:,8)>10&(filedata.params(:,2)>100|filedata.params(:,3)>100|filedata.params(:,4)>100|filedata.params(:,5)>100));
ind_t7=find(filedata.params(:,1)/10000>=tlin(1)&filedata.params(:,1)/10000<=tlin(end));
time_t7=filedata.params(ind_t7,1);
mark0_t7=[filedata.params(ind_t7,2) filedata.params(ind_t7,3) filedata.params(ind_t7,4) filedata.params(ind_t7,5)];
time2_t7=time_t7/10000;
spikeT0_t7=time2_t7;
[procInd0_t7,procInd1_t7]=histc(spikeT0_t7,tlin);
procInd_t7=find(procInd0_t7);
spikeT_t7=tlin(procInd_t7);
spike_t7=procInd0_t7';
[~,rawInd0_t7]=histc(spikeT0_t7,time2_t7);
markAll_t7(:,1)=procInd1_t7;markAll_t7(:,2:5)=mark0_t7(rawInd0_t7(rawInd0_t7~=0),:);
mdel=20; ms=min(min(markAll_t7(:,2:5))):mdel:max(max(markAll_t7(:,2:5))); 
sxker=2*xdel; smker = mdel; T=size(tlin,1);
occ=normpdf(xs'*ones(1,T),ones(length(xs),1)*xtrain,sxker)*ones(T,length(ms));
%occ: columns are identical; occupancy based on position; denominator
procInd1_t7_I1a=procInd1_t7(ismember(procInd1_t7,ind_I1));
procInd1_t7_I1=find(ismember(procInd1_t7,ind_I1));
Xnum_t7_I1=normpdf(xs'*ones(1,length(procInd1_t7_I1a)),ones(length(xs),1)*xtrain(procInd1_t7_I1a),sxker);
%Xnum: Gaussian kernel estimators for position
Lint_t7_I1=sum(Xnum_t7_I1,2)./occ(:,1)./dt; %integral
Lint_t7_I1=Lint_t7_I1./sum(Lint_t7_I1);
procInd1_t7_I0a=procInd1_t7(ismember(procInd1_t7,ind_I0));
procInd1_t7_I0=find(ismember(procInd1_t7,ind_I0));
Xnum_t7_I0=normpdf(xs'*ones(1,length(procInd1_t7_I0a)),ones(length(xs),1)*xtrain(procInd1_t7_I0a),sxker);
%Xnum: Gaussian kernel estimators for position
Lint_t7_I0=sum(Xnum_t7_I0,2)./occ(:,1)./dt; %integral
Lint_t7_I0=Lint_t7_I0./sum(Lint_t7_I0);

load('bond_data\bond04\10-146\bond04-10_params.mat');
%ind_t10=find(filedata.params(:,1)/10000>=tlin(1)&filedata.params(:,1)/10000<=tlin(end)&filedata.params(:,8)>10&(filedata.params(:,2)>100|filedata.params(:,3)>100|filedata.params(:,4)>100|filedata.params(:,5)>100));
ind_t10=find(filedata.params(:,1)/10000>=tlin(1)&filedata.params(:,1)/10000<=tlin(end));
time_t10=filedata.params(ind_t10,1);
mark0_t10=[filedata.params(ind_t10,2) filedata.params(ind_t10,3) filedata.params(ind_t10,4) filedata.params(ind_t10,5)];
time2_t10=time_t10/10000;
spikeT0_t10=time2_t10;
[procInd0_t10,procInd1_t10]=histc(spikeT0_t10,tlin);
procInd_t10=find(procInd0_t10);
spikeT_t10=tlin(procInd_t10);
spike_t10=procInd0_t10';
[~,rawInd0_t10]=histc(spikeT0_t10,time2_t10);
markAll_t10(:,1)=procInd1_t10;markAll_t10(:,2:5)=mark0_t10(rawInd0_t10(rawInd0_t10~=0),:);
mdel=20; ms=min(min(markAll_t10(:,2:5))):mdel:max(max(markAll_t10(:,2:5))); 
sxker=2*xdel; smker = mdel; T=size(tlin,1);
occ=normpdf(xs'*ones(1,T),ones(length(xs),1)*xtrain,sxker)*ones(T,length(ms));
%occ: columns are identical; occupancy based on position; denominator
procInd1_t10_I1a=procInd1_t10(ismember(procInd1_t10,ind_I1));
procInd1_t10_I1=find(ismember(procInd1_t10,ind_I1));
Xnum_t10_I1=normpdf(xs'*ones(1,length(procInd1_t10_I1a)),ones(length(xs),1)*xtrain(procInd1_t10_I1a),sxker);
%Xnum: Gaussian kernel estimators for position
Lint_t10_I1=sum(Xnum_t10_I1,2)./occ(:,1)./dt; %integral
Lint_t10_I1=Lint_t10_I1./sum(Lint_t10_I1);
procInd1_t10_I0a=procInd1_t10(ismember(procInd1_t10,ind_I0));
procInd1_t10_I0=find(ismember(procInd1_t10,ind_I0));
Xnum_t10_I0=normpdf(xs'*ones(1,length(procInd1_t10_I0a)),ones(length(xs),1)*xtrain(procInd1_t10_I0a),sxker);
%Xnum: Gaussian kernel estimators for position
Lint_t10_I0=sum(Xnum_t10_I0,2)./occ(:,1)./dt; %integral
Lint_t10_I0=Lint_t10_I0./sum(Lint_t10_I0);

load('bond_data\bond04\11-096\bond04-11_params.mat');
%ind_t11=find(filedata.params(:,1)/10000>=tlin(1)&filedata.params(:,1)/10000<=tlin(end)&filedata.params(:,8)>10&(filedata.params(:,2)>100|filedata.params(:,3)>100|filedata.params(:,4)>100|filedata.params(:,5)>100));
ind_t11=find(filedata.params(:,1)/10000>=tlin(1)&filedata.params(:,1)/10000<=tlin(end));
time_t11=filedata.params(ind_t11,1);
mark0_t11=[filedata.params(ind_t11,2) filedata.params(ind_t11,3) filedata.params(ind_t11,4) filedata.params(ind_t11,5)];
time2_t11=time_t11/10000;
spikeT0_t11=time2_t11;
[procInd0_t11,procInd1_t11]=histc(spikeT0_t11,tlin);
procInd_t11=find(procInd0_t11);
spikeT_t11=tlin(procInd_t11);
spike_t11=procInd0_t11';
[~,rawInd0_t11]=histc(spikeT0_t11,time2_t11);
markAll_t11(:,1)=procInd1_t11;markAll_t11(:,2:5)=mark0_t11(rawInd0_t11(rawInd0_t11~=0),:);
mdel=20; ms=min(min(markAll_t11(:,2:5))):mdel:max(max(markAll_t11(:,2:5)));
sxker=2*xdel; smker = mdel; T=size(tlin,1);
occ=normpdf(xs'*ones(1,T),ones(length(xs),1)*xtrain,sxker)*ones(T,length(ms));
%occ: columns are identical; occupancy based on position; denominator
procInd1_t11_I1a=procInd1_t11(ismember(procInd1_t11,ind_I1));
procInd1_t11_I1=find(ismember(procInd1_t11,ind_I1));
Xnum_t11_I1=normpdf(xs'*ones(1,length(procInd1_t11_I1a)),ones(length(xs),1)*xtrain(procInd1_t11_I1a),sxker);
%Xnum: Gaussian kernel estimators for position
Lint_t11_I1=sum(Xnum_t11_I1,2)./occ(:,1)./dt; %integral
Lint_t11_I1=Lint_t11_I1./sum(Lint_t11_I1);
procInd1_t11_I0a=procInd1_t11(ismember(procInd1_t11,ind_I0));
procInd1_t11_I0=find(ismember(procInd1_t11,ind_I0));
Xnum_t11_I0=normpdf(xs'*ones(1,length(procInd1_t11_I0a)),ones(length(xs),1)*xtrain(procInd1_t11_I0a),sxker);
%Xnum: Gaussian kernel estimators for position
Lint_t11_I0=sum(Xnum_t11_I0,2)./occ(:,1)./dt; %integral
Lint_t11_I0=Lint_t11_I0./sum(Lint_t11_I0);

load('bond_data\bond04\12-136\bond04-12_params.mat');
%ind_t12=find(filedata.params(:,1)/10000>=tlin(1)&filedata.params(:,1)/10000<=tlin(end)&filedata.params(:,8)>10&(filedata.params(:,2)>100|filedata.params(:,3)>100|filedata.params(:,4)>100|filedata.params(:,5)>100));
ind_t12=find(filedata.params(:,1)/10000>=tlin(1)&filedata.params(:,1)/10000<=tlin(end));
time_t12=filedata.params(ind_t12,1);
mark0_t12=[filedata.params(ind_t12,2) filedata.params(ind_t12,3) filedata.params(ind_t12,4) filedata.params(ind_t12,5)];
time2_t12=time_t12/10000;
spikeT0_t12=time2_t12;
[procInd0_t12,procInd1_t12]=histc(spikeT0_t12,tlin);
procInd_t12=find(procInd0_t12);
spikeT_t12=tlin(procInd_t12);
spike_t12=procInd0_t12';
[~,rawInd0_t12]=histc(spikeT0_t12,time2_t12);
markAll_t12(:,1)=procInd1_t12;markAll_t12(:,2:5)=mark0_t12(rawInd0_t12(rawInd0_t12~=0),:);
mdel=20; ms=min(min(markAll_t12(:,2:5))):mdel:max(max(markAll_t12(:,2:5)));
sxker=2*xdel; smker = mdel; T=size(tlin,1);
occ=normpdf(xs'*ones(1,T),ones(length(xs),1)*xtrain,sxker)*ones(T,length(ms));
%occ: columns are identical; occupancy based on position; denominator
procInd1_t12_I1a=procInd1_t12(ismember(procInd1_t12,ind_I1));
procInd1_t12_I1=find(ismember(procInd1_t12,ind_I1));
Xnum_t12_I1=normpdf(xs'*ones(1,length(procInd1_t12_I1a)),ones(length(xs),1)*xtrain(procInd1_t12_I1a),sxker);
%Xnum: Gaussian kernel estimators for position
Lint_t12_I1=sum(Xnum_t12_I1,2)./occ(:,1)./dt; %integral
Lint_t12_I1=Lint_t12_I1./sum(Lint_t12_I1);
procInd1_t12_I0a=procInd1_t12(ismember(procInd1_t12,ind_I0));
procInd1_t12_I0=find(ismember(procInd1_t12,ind_I0));
Xnum_t12_I0=normpdf(xs'*ones(1,length(procInd1_t12_I0a)),ones(length(xs),1)*xtrain(procInd1_t12_I0a),sxker);
%Xnum: Gaussian kernel estimators for position
Lint_t12_I0=sum(Xnum_t12_I0,2)./occ(:,1)./dt; %integral
Lint_t12_I0=Lint_t12_I0./sum(Lint_t12_I0);

load('bond_data\bond04\13-094\bond04-13_params.mat');
%ind_t13=find(filedata.params(:,1)/10000>=tlin(1)&filedata.params(:,1)/10000<=tlin(end)&filedata.params(:,8)>10&(filedata.params(:,2)>100|filedata.params(:,3)>100|filedata.params(:,4)>100|filedata.params(:,5)>100));
ind_t13=find(filedata.params(:,1)/10000>=tlin(1)&filedata.params(:,1)/10000<=tlin(end));
time_t13=filedata.params(ind_t13,1);
mark0_t13=[filedata.params(ind_t13,2) filedata.params(ind_t13,3) filedata.params(ind_t13,4) filedata.params(ind_t13,5)];
time2_t13=time_t13/10000;
spikeT0_t13=time2_t13;
[procInd0_t13,procInd1_t13]=histc(spikeT0_t13,tlin);
procInd_t13=find(procInd0_t13);
spikeT_t13=tlin(procInd_t13);
spike_t13=procInd0_t13';
[~,rawInd0_t13]=histc(spikeT0_t13,time2_t13);
markAll_t13(:,1)=procInd1_t13;markAll_t13(:,2:5)=mark0_t13(rawInd0_t13(rawInd0_t13~=0),:);
mdel=20; ms=min(min(markAll_t13(:,2:5))):mdel:max(max(markAll_t13(:,2:5)));
sxker=2*xdel; smker = mdel; T=size(tlin,1);
occ=normpdf(xs'*ones(1,T),ones(length(xs),1)*xtrain,sxker)*ones(T,length(ms));
%occ: columns are identical; occupancy based on position; denominator
procInd1_t13_I1a=procInd1_t13(ismember(procInd1_t13,ind_I1));
procInd1_t13_I1=find(ismember(procInd1_t13,ind_I1));
Xnum_t13_I1=normpdf(xs'*ones(1,length(procInd1_t13_I1a)),ones(length(xs),1)*xtrain(procInd1_t13_I1a),sxker);
%Xnum: Gaussian kernel estimators for position
Lint_t13_I1=sum(Xnum_t13_I1,2)./occ(:,1)./dt; %integral
Lint_t13_I1=Lint_t13_I1./sum(Lint_t13_I1);
procInd1_t13_I0a=procInd1_t13(ismember(procInd1_t13,ind_I0));
procInd1_t13_I0=find(ismember(procInd1_t13,ind_I0));
Xnum_t13_I0=normpdf(xs'*ones(1,length(procInd1_t13_I0a)),ones(length(xs),1)*xtrain(procInd1_t13_I0a),sxker);
%Xnum: Gaussian kernel estimators for position
Lint_t13_I0=sum(Xnum_t13_I0,2)./occ(:,1)./dt; %integral
Lint_t13_I0=Lint_t13_I0./sum(Lint_t13_I0);

load('bond_data\bond04\14-096\bond04-14_params.mat');
%ind_t14=find(filedata.params(:,1)/10000>=tlin(1)&filedata.params(:,1)/10000<=tlin(end)&filedata.params(:,8)>10&(filedata.params(:,2)>100|filedata.params(:,3)>100|filedata.params(:,4)>100|filedata.params(:,5)>100));
ind_t14=find(filedata.params(:,1)/10000>=tlin(1)&filedata.params(:,1)/10000<=tlin(end));
time_t14=filedata.params(ind_t14,1);
mark0_t14=[filedata.params(ind_t14,2) filedata.params(ind_t14,3) filedata.params(ind_t14,4) filedata.params(ind_t14,5)];
time2_t14=time_t14/10000;
spikeT0_t14=time2_t14;
[procInd0_t14,procInd1_t14]=histc(spikeT0_t14,tlin);
procInd_t14=find(procInd0_t14);
spikeT_t14=tlin(procInd_t14);
spike_t14=procInd0_t14';
[~,rawInd0_t14]=histc(spikeT0_t14,time2_t14);
markAll_t14(:,1)=procInd1_t14;markAll_t14(:,2:5)=mark0_t14(rawInd0_t14(rawInd0_t14~=0),:);
mdel=20; ms=min(min(markAll_t14(:,2:5))):mdel:max(max(markAll_t14(:,2:5)));
sxker=2*xdel; smker = mdel; T=size(tlin,1);
occ=normpdf(xs'*ones(1,T),ones(length(xs),1)*xtrain,sxker)*ones(T,length(ms));
%occ: columns are identical; occupancy based on position; denominator
procInd1_t14_I1a=procInd1_t14(ismember(procInd1_t14,ind_I1));
procInd1_t14_I1=find(ismember(procInd1_t14,ind_I1));
Xnum_t14_I1=normpdf(xs'*ones(1,length(procInd1_t14_I1a)),ones(length(xs),1)*xtrain(procInd1_t14_I1a),sxker);
%Xnum: Gaussian kernel estimators for position
Lint_t14_I1=sum(Xnum_t14_I1,2)./occ(:,1)./dt; %integral
Lint_t14_I1=Lint_t14_I1./sum(Lint_t14_I1);
procInd1_t14_I0a=procInd1_t14(ismember(procInd1_t14,ind_I0));
procInd1_t14_I0=find(ismember(procInd1_t14,ind_I0));
Xnum_t14_I0=normpdf(xs'*ones(1,length(procInd1_t14_I0a)),ones(length(xs),1)*xtrain(procInd1_t14_I0a),sxker);
%Xnum: Gaussian kernel estimators for position
Lint_t14_I0=sum(Xnum_t14_I0,2)./occ(:,1)./dt; %integral
Lint_t14_I0=Lint_t14_I0./sum(Lint_t14_I0);

load('bond_data\bond04\17-108\bond04-17_params.mat');
%ind_t17=find(filedata.params(:,1)/10000>=tlin(1)&filedata.params(:,1)/10000<=tlin(end)&filedata.params(:,8)>10&(filedata.params(:,2)>100|filedata.params(:,3)>100|filedata.params(:,4)>100|filedata.params(:,5)>100));
ind_t17=find(filedata.params(:,1)/10000>=tlin(1)&filedata.params(:,1)/10000<=tlin(end));
time_t17=filedata.params(ind_t17,1);
mark0_t17=[filedata.params(ind_t17,2) filedata.params(ind_t17,3) filedata.params(ind_t17,4) filedata.params(ind_t17,5)];
time2_t17=time_t17/10000;
spikeT0_t17=time2_t17;
[procInd0_t17,procInd1_t17]=histc(spikeT0_t17,tlin);
procInd_t17=find(procInd0_t17);
spikeT_t17=tlin(procInd_t17);
spike_t17=procInd0_t17';
[~,rawInd0_t17]=histc(spikeT0_t17,time2_t17);
markAll_t17(:,1)=procInd1_t17;markAll_t17(:,2:5)=mark0_t17(rawInd0_t17(rawInd0_t17~=0),:);
mdel=20; ms=min(min(markAll_t17(:,2:5))):mdel:max(max(markAll_t17(:,2:5)));
sxker=2*xdel; smker = mdel; T=size(tlin,1);
occ=normpdf(xs'*ones(1,T),ones(length(xs),1)*xtrain,sxker)*ones(T,length(ms));
%occ: columns are identical; occupancy based on position; denominator
procInd1_t17_I1a=procInd1_t17(ismember(procInd1_t17,ind_I1));
procInd1_t17_I1=find(ismember(procInd1_t17,ind_I1));
Xnum_t17_I1=normpdf(xs'*ones(1,length(procInd1_t17_I1a)),ones(length(xs),1)*xtrain(procInd1_t17_I1a),sxker);
%Xnum: Gaussian kernel estimators for position
Lint_t17_I1=sum(Xnum_t17_I1,2)./occ(:,1)./dt; %integral
Lint_t17_I1=Lint_t17_I1./sum(Lint_t17_I1);
procInd1_t17_I0a=procInd1_t17(ismember(procInd1_t17,ind_I0));
procInd1_t17_I0=find(ismember(procInd1_t17,ind_I0));
Xnum_t17_I0=normpdf(xs'*ones(1,length(procInd1_t17_I0a)),ones(length(xs),1)*xtrain(procInd1_t17_I0a),sxker);
%Xnum: Gaussian kernel estimators for position
Lint_t17_I0=sum(Xnum_t17_I0,2)./occ(:,1)./dt; %integral
Lint_t17_I0=Lint_t17_I0./sum(Lint_t17_I0);

load('bond_data\bond04\18-138\bond04-18_params.mat');
%ind_t18=find(filedata.params(:,1)/10000>=tlin(1)&filedata.params(:,1)/10000<=tlin(end)&filedata.params(:,8)>10&(filedata.params(:,2)>100|filedata.params(:,3)>100|filedata.params(:,4)>100|filedata.params(:,5)>100));
ind_t18=find(filedata.params(:,1)/10000>=tlin(1)&filedata.params(:,1)/10000<=tlin(end));
time_t18=filedata.params(ind_t18,1);
mark0_t18=[filedata.params(ind_t18,2) filedata.params(ind_t18,3) filedata.params(ind_t18,4) filedata.params(ind_t18,5)];
time2_t18=time_t18/10000;
spikeT0_t18=time2_t18;
[procInd0_t18,procInd1_t18]=histc(spikeT0_t18,tlin);
procInd_t18=find(procInd0_t18);
spikeT_t18=tlin(procInd_t18);
spike_t18=procInd0_t18';
[~,rawInd0_t18]=histc(spikeT0_t18,time2_t18);
markAll_t18(:,1)=procInd1_t18;markAll_t18(:,2:5)=mark0_t18(rawInd0_t18(rawInd0_t18~=0),:);
mdel=20; ms=min(min(markAll_t18(:,2:5))):mdel:max(max(markAll_t18(:,2:5)));
sxker=2*xdel; smker = mdel; T=size(tlin,1);
occ=normpdf(xs'*ones(1,T),ones(length(xs),1)*xtrain,sxker)*ones(T,length(ms));
%occ: columns are identical; occupancy based on position; denominator
procInd1_t18_I1a=procInd1_t18(ismember(procInd1_t18,ind_I1));
procInd1_t18_I1=find(ismember(procInd1_t18,ind_I1));
Xnum_t18_I1=normpdf(xs'*ones(1,length(procInd1_t18_I1a)),ones(length(xs),1)*xtrain(procInd1_t18_I1a),sxker);
%Xnum: Gaussian kernel estimators for position
Lint_t18_I1=sum(Xnum_t18_I1,2)./occ(:,1)./dt; %integral
Lint_t18_I1=Lint_t18_I1./sum(Lint_t18_I1);
procInd1_t18_I0a=procInd1_t18(ismember(procInd1_t18,ind_I0));
procInd1_t18_I0=find(ismember(procInd1_t18,ind_I0));
Xnum_t18_I0=normpdf(xs'*ones(1,length(procInd1_t18_I0a)),ones(length(xs),1)*xtrain(procInd1_t18_I0a),sxker);
%Xnum: Gaussian kernel estimators for position
Lint_t18_I0=sum(Xnum_t18_I0,2)./occ(:,1)./dt; %integral
Lint_t18_I0=Lint_t18_I0./sum(Lint_t18_I0);

load('bond_data\bond04\19-130\bond04-19_params.mat');
%ind_t19=find(filedata.params(:,1)/10000>=tlin(1)&filedata.params(:,1)/10000<=tlin(end)&filedata.params(:,8)>10&(filedata.params(:,2)>100|filedata.params(:,3)>100|filedata.params(:,4)>100|filedata.params(:,5)>100));
ind_t19=find(filedata.params(:,1)/10000>=tlin(1)&filedata.params(:,1)/10000<=tlin(end));
time_t19=filedata.params(ind_t19,1);
mark0_t19=[filedata.params(ind_t19,2) filedata.params(ind_t19,3) filedata.params(ind_t19,4) filedata.params(ind_t19,5)];
time2_t19=time_t19/10000;
spikeT0_t19=time2_t19;
[procInd0_t19,procInd1_t19]=histc(spikeT0_t19,tlin);
procInd_t19=find(procInd0_t19);
spikeT_t19=tlin(procInd_t19);
spike_t19=procInd0_t19';
[~,rawInd0_t19]=histc(spikeT0_t19,time2_t19);
markAll_t19(:,1)=procInd1_t19;markAll_t19(:,2:5)=mark0_t19(rawInd0_t19(rawInd0_t19~=0),:);
mdel=20; ms=min(min(markAll_t19(:,2:5))):mdel:max(max(markAll_t19(:,2:5)));
sxker=2*xdel; smker = mdel; T=size(tlin,1);
occ=normpdf(xs'*ones(1,T),ones(length(xs),1)*xtrain,sxker)*ones(T,length(ms));
%occ: columns are identical; occupancy based on position; denominator
procInd1_t19_I1a=procInd1_t19(ismember(procInd1_t19,ind_I1));
procInd1_t19_I1=find(ismember(procInd1_t19,ind_I1));
Xnum_t19_I1=normpdf(xs'*ones(1,length(procInd1_t19_I1a)),ones(length(xs),1)*xtrain(procInd1_t19_I1a),sxker);
%Xnum: Gaussian kernel estimators for position
Lint_t19_I1=sum(Xnum_t19_I1,2)./occ(:,1)./dt; %integral
Lint_t19_I1=Lint_t19_I1./sum(Lint_t19_I1);
procInd1_t19_I0a=procInd1_t19(ismember(procInd1_t19,ind_I1));
procInd1_t19_I0=find(ismember(procInd1_t19,ind_I1));
Xnum_t19_I0=normpdf(xs'*ones(1,length(procInd1_t19_I0a)),ones(length(xs),1)*xtrain(procInd1_t19_I0a),sxker);
%Xnum: Gaussian kernel estimators for position
Lint_t19_I0=sum(Xnum_t19_I0,2)./occ(:,1)./dt; %integral
Lint_t19_I0=Lint_t19_I0./sum(Lint_t19_I0);

load('bond_data\bond04\20-132\bond04-20_params.mat');
%ind_t20=find(filedata.params(:,1)/10000>=tlin(1)&filedata.params(:,1)/10000<=tlin(end)&filedata.params(:,8)>10&(filedata.params(:,2)>100|filedata.params(:,3)>100|filedata.params(:,4)>100|filedata.params(:,5)>100));
ind_t20=find(filedata.params(:,1)/10000>=tlin(1)&filedata.params(:,1)/10000<=tlin(end));
time_t20=filedata.params(ind_t20,1);
mark0_t20=[filedata.params(ind_t20,2) filedata.params(ind_t20,3) filedata.params(ind_t20,4) filedata.params(ind_t20,5)];
time2_t20=time_t20/10000;
spikeT0_t20=time2_t20;
[procInd0_t20,procInd1_t20]=histc(spikeT0_t20,tlin);
procInd_t20=find(procInd0_t20);
spikeT_t20=tlin(procInd_t20);
spike_t20=procInd0_t20';
[~,rawInd0_t20]=histc(spikeT0_t20,time2_t20);
markAll_t20(:,1)=procInd1_t20;markAll_t20(:,2:5)=mark0_t20(rawInd0_t20(rawInd0_t20~=0),:);
mdel=20; ms=min(min(markAll_t20(:,2:5))):mdel:max(max(markAll_t20(:,2:5)));
sxker=2*xdel; smker = mdel; T=size(tlin,1);
occ=normpdf(xs'*ones(1,T),ones(length(xs),1)*xtrain,sxker)*ones(T,length(ms));
%occ: columns are identical; occupancy based on position; denominator
procInd1_t20_I1a=procInd1_t20(ismember(procInd1_t20,ind_I1));
procInd1_t20_I1=find(ismember(procInd1_t20,ind_I1));
Xnum_t20_I1=normpdf(xs'*ones(1,length(procInd1_t20_I1a)),ones(length(xs),1)*xtrain(procInd1_t20_I1a),sxker);
%Xnum: Gaussian kernel estimators for position
Lint_t20_I1=sum(Xnum_t20_I1,2)./occ(:,1)./dt; %integral
Lint_t20_I1=Lint_t20_I1./sum(Lint_t20_I1);
procInd1_t20_I0a=procInd1_t20(ismember(procInd1_t20,ind_I0));
procInd1_t20_I0=find(ismember(procInd1_t20,ind_I0));
Xnum_t20_I0=normpdf(xs'*ones(1,length(procInd1_t20_I0a)),ones(length(xs),1)*xtrain(procInd1_t20_I0a),sxker);
%Xnum: Gaussian kernel estimators for position
Lint_t20_I0=sum(Xnum_t20_I0,2)./occ(:,1)./dt; %integral
Lint_t20_I0=Lint_t20_I0./sum(Lint_t20_I0);

load('bond_data\bond04\22-171\bond04-22_params.mat');
%ind_t22=find(filedata.params(:,1)/10000>=tlin(1)&filedata.params(:,1)/10000<=tlin(end)&filedata.params(:,8)>10&(filedata.params(:,2)>100|filedata.params(:,3)>100|filedata.params(:,4)>100|filedata.params(:,5)>100));
ind_t22=find(filedata.params(:,1)/10000>=tlin(1)&filedata.params(:,1)/10000<=tlin(end));
time_t22=filedata.params(ind_t22,1);
mark0_t22=[filedata.params(ind_t22,2) filedata.params(ind_t22,3) filedata.params(ind_t22,4) filedata.params(ind_t22,5)];
time2_t22=time_t22/10000;
spikeT0_t22=time2_t22;
[procInd0_t22,procInd1_t22]=histc(spikeT0_t22,tlin);
procInd_t22=find(procInd0_t22);
spikeT_t22=tlin(procInd_t22);
spike_t22=procInd0_t22';
[~,rawInd0_t22]=histc(spikeT0_t22,time2_t22);
markAll_t22(:,1)=procInd1_t22;markAll_t22(:,2:5)=mark0_t22(rawInd0_t22(rawInd0_t22~=0),:);
mdel=20; ms=min(min(markAll_t22(:,2:5))):mdel:max(max(markAll_t22(:,2:5)));
sxker=2*xdel; smker = mdel; T=size(tlin,1);
occ=normpdf(xs'*ones(1,T),ones(length(xs),1)*xtrain,sxker)*ones(T,length(ms));
%occ: columns are identical; occupancy based on position; denominator
procInd1_t22_I1a=procInd1_t22(ismember(procInd1_t22,ind_I1));
procInd1_t22_I1=find(ismember(procInd1_t22,ind_I1));
Xnum_t22_I1=normpdf(xs'*ones(1,length(procInd1_t22_I1a)),ones(length(xs),1)*xtrain(procInd1_t22_I1a),sxker);
%Xnum: Gaussian kernel estimators for position
Lint_t22_I1=sum(Xnum_t22_I1,2)./occ(:,1)./dt; %integral
Lint_t22_I1=Lint_t22_I1./sum(Lint_t22_I1);
procInd1_t22_I0a=procInd1_t22(ismember(procInd1_t22,ind_I1));
procInd1_t22_I0=find(ismember(procInd1_t22,ind_I1));
Xnum_t22_I0=normpdf(xs'*ones(1,length(procInd1_t22_I0a)),ones(length(xs),1)*xtrain(procInd1_t22_I0a),sxker);
%Xnum: Gaussian kernel estimators for position
Lint_t22_I0=sum(Xnum_t22_I0,2)./occ(:,1)./dt; %integral
Lint_t22_I0=Lint_t22_I0./sum(Lint_t22_I0);

load('bond_data\bond04\23-152\bond04-23_params.mat');
%ind_t23=find(filedata.params(:,1)/10000>=tlin(1)&filedata.params(:,1)/10000<=tlin(end)&filedata.params(:,8)>10&(filedata.params(:,2)>100|filedata.params(:,3)>100|filedata.params(:,4)>100|filedata.params(:,5)>100));
ind_t23=find(filedata.params(:,1)/10000>=tlin(1)&filedata.params(:,1)/10000<=tlin(end));
time_t23=filedata.params(ind_t23,1);
mark0_t23=[filedata.params(ind_t23,2) filedata.params(ind_t23,3) filedata.params(ind_t23,4) filedata.params(ind_t23,5)];
time2_t23=time_t23/10000;
spikeT0_t23=time2_t23;
[procInd0_t23,procInd1_t23]=histc(spikeT0_t23,tlin);
procInd_t23=find(procInd0_t23);
spikeT_t23=tlin(procInd_t23);
spike_t23=procInd0_t23';
[~,rawInd0_t23]=histc(spikeT0_t23,time2_t23);
markAll_t23(:,1)=procInd1_t23;markAll_t23(:,2:5)=mark0_t23(rawInd0_t23(rawInd0_t23~=0),:);
mdel=20; ms=min(min(markAll_t23(:,2:5))):mdel:max(max(markAll_t23(:,2:5)));
sxker=2*xdel; smker = mdel; T=size(tlin,1);
occ=normpdf(xs'*ones(1,T),ones(length(xs),1)*xtrain,sxker)*ones(T,length(ms));
%occ: columns are identical; occupancy based on position; denominator
procInd1_t23_I1a=procInd1_t23(ismember(procInd1_t23,ind_I1));
procInd1_t23_I1=find(ismember(procInd1_t23,ind_I1));
Xnum_t23_I1=normpdf(xs'*ones(1,length(procInd1_t23_I1a)),ones(length(xs),1)*xtrain(procInd1_t23_I1a),sxker);
%Xnum: Gaussian kernel estimators for position
Lint_t23_I1=sum(Xnum_t23_I1,2)./occ(:,1)./dt; %integral
Lint_t23_I1=Lint_t23_I1./sum(Lint_t23_I1);
procInd1_t23_I0a=procInd1_t23(ismember(procInd1_t23,ind_I0));
procInd1_t23_I0=find(ismember(procInd1_t23,ind_I0));
Xnum_t23_I0=normpdf(xs'*ones(1,length(procInd1_t23_I0a)),ones(length(xs),1)*xtrain(procInd1_t23_I0a),sxker);
%Xnum: Gaussian kernel estimators for position
Lint_t23_I0=sum(Xnum_t23_I0,2)./occ(:,1)./dt; %integral
Lint_t23_I0=Lint_t23_I0./sum(Lint_t23_I0);

load('bond_data\bond04\27-158\bond04-27_params.mat');
%ind_t27=find(filedata.params(:,1)/10000>=tlin(1)&filedata.params(:,1)/10000<=tlin(end)&filedata.params(:,8)>10&(filedata.params(:,2)>100|filedata.params(:,3)>100|filedata.params(:,4)>100|filedata.params(:,5)>100));
ind_t27=find(filedata.params(:,1)/10000>=tlin(1)&filedata.params(:,1)/10000<=tlin(end));
time_t27=filedata.params(ind_t27,1);
mark0_t27=[filedata.params(ind_t27,2) filedata.params(ind_t27,3) filedata.params(ind_t27,4) filedata.params(ind_t27,5)];
time2_t27=time_t27/10000;
spikeT0_t27=time2_t27;
[procInd0_t27,procInd1_t27]=histc(spikeT0_t27,tlin);
procInd_t27=find(procInd0_t27);
spikeT_t27=tlin(procInd_t27);
spike_t27=procInd0_t27';
[~,rawInd0_t27]=histc(spikeT0_t27,time2_t27);
markAll_t27(:,1)=procInd1_t27;markAll_t27(:,2:5)=mark0_t27(rawInd0_t27(rawInd0_t27~=0),:);
mdel=20; ms=min(min(markAll_t27(:,2:5))):mdel:max(max(markAll_t27(:,2:5)));
sxker=2*xdel; smker = mdel; T=size(tlin,1);
occ=normpdf(xs'*ones(1,T),ones(length(xs),1)*xtrain,sxker)*ones(T,length(ms));
%occ: columns are identical; occupancy based on position; denominator
procInd1_t27_I1a=procInd1_t27(ismember(procInd1_t27,ind_I1));
procInd1_t27_I1=find(ismember(procInd1_t27,ind_I1));
Xnum_t27_I1=normpdf(xs'*ones(1,length(procInd1_t27_I1a)),ones(length(xs),1)*xtrain(procInd1_t27_I1a),sxker);
%Xnum: Gaussian kernel estimators for position
Lint_t27_I1=sum(Xnum_t27_I1,2)./occ(:,1)./dt; %integral
Lint_t27_I1=Lint_t27_I1./sum(Lint_t27_I1);
procInd1_t27_I0a=procInd1_t27(ismember(procInd1_t27,ind_I0));
procInd1_t27_I0=find(ismember(procInd1_t27,ind_I0));
Xnum_t27_I0=normpdf(xs'*ones(1,length(procInd1_t27_I0a)),ones(length(xs),1)*xtrain(procInd1_t27_I0a),sxker);
%Xnum: Gaussian kernel estimators for position
Lint_t27_I0=sum(Xnum_t27_I0,2)./occ(:,1)./dt; %integral
Lint_t27_I0=Lint_t27_I0./sum(Lint_t27_I0);

load('bond_data\bond04\29-118\bond04-29_params.mat');
%ind_t29=find(filedata.params(:,1)/10000>=tlin(1)&filedata.params(:,1)/10000<=tlin(end)&filedata.params(:,8)>10&(filedata.params(:,2)>100|filedata.params(:,3)>100|filedata.params(:,4)>100|filedata.params(:,5)>100));
ind_t29=find(filedata.params(:,1)/10000>=tlin(1)&filedata.params(:,1)/10000<=tlin(end));
time_t29=filedata.params(ind_t29,1);
mark0_t29=[filedata.params(ind_t29,2) filedata.params(ind_t29,3) filedata.params(ind_t29,4) filedata.params(ind_t29,5)];
time2_t29=time_t29/10000;
spikeT0_t29=time2_t29;
[procInd0_t29,procInd1_t29]=histc(spikeT0_t29,tlin);
procInd_t29=find(procInd0_t29);
spikeT_t29=tlin(procInd_t29);
spike_t29=procInd0_t29';
[~,rawInd0_t29]=histc(spikeT0_t29,time2_t29);
markAll_t29(:,1)=procInd1_t29;markAll_t29(:,2:5)=mark0_t29(rawInd0_t29(rawInd0_t29~=0),:);
mdel=20; ms=min(min(markAll_t29(:,2:5))):mdel:max(max(markAll_t29(:,2:5)));
sxker=2*xdel; smker = mdel; T=size(tlin,1);
occ=normpdf(xs'*ones(1,T),ones(length(xs),1)*xtrain,sxker)*ones(T,length(ms));
%occ: columns are identical; occupancy based on position; denominator
procInd1_t29_I1a=procInd1_t29(ismember(procInd1_t29,ind_I1));
procInd1_t29_I1=find(ismember(procInd1_t29,ind_I1));
Xnum_t29_I1=normpdf(xs'*ones(1,length(procInd1_t29_I1a)),ones(length(xs),1)*xtrain(procInd1_t29_I1a),sxker);
%Xnum: Gaussian kernel estimators for position
Lint_t29_I1=sum(Xnum_t29_I1,2)./occ(:,1)./dt; %integral
Lint_t29_I1=Lint_t29_I1./sum(Lint_t29_I1);
procInd1_t29_I0a=procInd1_t29(ismember(procInd1_t29,ind_I0));
procInd1_t29_I0=find(ismember(procInd1_t29,ind_I0));
Xnum_t29_I0=normpdf(xs'*ones(1,length(procInd1_t29_I0a)),ones(length(xs),1)*xtrain(procInd1_t29_I0a),sxker);
%Xnum: Gaussian kernel estimators for position
Lint_t29_I0=sum(Xnum_t29_I0,2)./occ(:,1)./dt; %integral
Lint_t29_I0=Lint_t29_I0./sum(Lint_t29_I0);