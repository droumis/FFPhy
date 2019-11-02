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

%%
X=[];
s=0.5;
c_pt_lin=linspace(-3, 3.0001, 15); %13 bins for frank; 15 for bond
c_pt_lin_all=[0.2 c_pt_lin -0.2];
num_c_pts = length(c_pt_lin_all);
vecInd=find(vecLF(:,1)~=0);
for k0=1:length(vecInd)
    k=vecInd(k0);
    nearest_c_pt_index = find(c_pt_lin<=vecLF(k,2), 1, 'last')+1;
    nearest_c_pt_lin = c_pt_lin_all(nearest_c_pt_index);
    next_c_pt_lin = c_pt_lin_all(nearest_c_pt_index+1);
    u = (vecLF(k,2)-nearest_c_pt_lin)/(next_c_pt_lin-nearest_c_pt_lin);
    p=[u^3 u^2 u 1]*[-s 2-s s-2 s;2*s s-3 3-2*s -s;-s 0 s 0;0 1 0 0];
    X(k0,:) = [zeros(1,nearest_c_pt_index-2) p zeros(1,num_c_pts-4-(nearest_c_pt_index-2))];
end


newX=[];
stateV=linspace(-3,3,121);
num_c_pts = length(c_pt_lin_all);
for k=1:size(stateV,2)
    nearest_c_pt_index = find(c_pt_lin<=stateV(k), 1, 'last')+1;
    nearest_c_pt_lin = c_pt_lin_all(nearest_c_pt_index);
    next_c_pt_lin = c_pt_lin_all(nearest_c_pt_index+1);
    u = (stateV(k)-nearest_c_pt_lin)/(next_c_pt_lin-nearest_c_pt_lin);
    p=[u^3 u^2 u 1]*[-s 2-s s-2 s;2*s s-3 3-2*s -s;-s 0 s 0;0 1 0 0];
    newX(k,:) = [zeros(1,nearest_c_pt_index-2) p zeros(1,num_c_pts-4-(nearest_c_pt_index-2))];
end


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
lambda_non=[]; spike=[];
spRate=zeros(size(a,2),1);
tic;
for kk=1:size(a,2)
    j=a(kk);i=b(kk);
    B=spikes{ex}{ep}{j}{i}.data(:,1); %spiking times for tetrode j, cell i
    [sptrain0,bin]=histc(B,A);

    spRate(kk)=sum(sptrain0)/size(sptrain0,1);
    
    %X=X(find(vecLF(:,1)~=0),:);
    sptrain0=sptrain0(find(vecLF(:,1)~=0),:);
    spike=[spike;sptrain0'];

    betahat=glmfit(X,sptrain0,'poisson','constant','off');
    
    lambda_pause=glmval(betahat,newX,'log','constant','off');
    %lambda_pause: place field for non-ripples for individual cell
    lambda_non=[lambda_non lambda_pause]; %stacking place field by cells
    clear sptrain0 lambda_pause;
end
toc


%%
segT_I1=[2520 2565;2647 2705;2748 2777;2835 2856;2880 2909;2959 2975;3023.5 3032;3078 3097;3139 3179;3241.5 3255.5;3316.4 3330.6];
segT_I0=[2585 2620;2705 2730;2795 2820;2858 2880;2909 2936;2991 3023.5;3037 3062;3113.5 3209;3190.5 3227.5;3274.5 3301;3340 3387.5];

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

ind_I=zeros(size(vecLF,1),1);
ind_I(ind_I1)=1;

lambda_non_I1=zeros(length(stateV),size(a,2));
lambda_non_I0=lambda_non_I1;
tic;
for kk=1:size(a,2)
    j=a(kk);i=b(kk);
    B=spikes{ex}{ep}{j}{i}.data(:,1); %spiking times for tetrode j, cell i
    [sptrain0,bin]=histc(B,A);
    
    block1Time=repmat(ind_I,1,size(X,2));

    %%%%
    [betahat,~,stats]=glmfit([X X.*block1Time],sptrain0,'poisson','constant','off');
    
    [lambdaI1,dyloI1,dyhiI1]=glmval(betahat,[newX newX],'log',stats,'constant','off');
    [lambdaI0,dyloI0,dyhiI0]=glmval(betahat,[newX zeros(size(newX,1),size(newX,2))],'log',stats,'constant','off');

    lambda_non_I1(:,kk)=lambdaI1;
    lambda_non_I0(:,kk)=lambdaI0;

    clear betahat lambda dylo dyhi
end
toc

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



%%
%%%%ripples
%%%%%%new ripple times by Kenny; April 15, 2015
%%%%%%Create indicator function of ripples for 1ms binned spike trains
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
rippeI=find(sumR>10);
%length(kk)



%%
%%%%choose ripple index
%%%%set stdev value and number of bins for numeric integral

rIndV=10;
spike_r=spike_r_all{rippeI(rIndV)};



%%
%%%%1. initial different, state and obs identical

clear onestep_I1 onestep_I0 postx_I1 postx_I0 postxM_r;
C=size(spike_r,1);
dt=1/33.4;
n=length(stateV);numSteps=size(spike_r,2);
%P(x0|I)
Px_I1=exp(-stateV.^2./(2*0.2^2));Px_I1=Px_I1./sum(Px_I1);
Px_I0=max(Px_I1)*ones(1,n)-Px_I1; Px_I0=Px_I0./sum(Px_I0);
%P(x0)=P(x0|I)P(I);
postx_I1=0.5*Px_I1';postx_I0=0.5*Px_I0';
pI1_vec=zeros(numSteps,1);pI0_vec=zeros(numSteps,1);
for t=1:numSteps
    onestep_I1=stateM_gausnorm*postx_I1;
    onestep_I0=stateM_gausnorm*postx_I0;
    %figure;plot(stateV,onestep_I1,'r',stateV,onestep_I0,'b');
    L=ones(n,1);
    for c=1:C
        L=L.*((lambda_non(:,c).*dt).^spike_r(c,t).*exp(-lambda_non(:,c).*dt));
    end
    totnorm=sum(onestep_I1.*L)+sum(onestep_I0.*L);
    postx_I1=onestep_I1.*L./totnorm;
    postx_I0=onestep_I0.*L./totnorm;
    pI1_vec(t)=sum(postx_I1);
    pI0_vec(t)=sum(postx_I0);
    clear onestep_I1 onestep_I2
end

figure;
subplot(3,3,1:2);
imagesc(spike_r);colormap(flipud(gray));freezeColors
set(gca,'YTick',1:5:30,'YTickLabel',1:5:30,'FontSize',14);
title(['ripple = ',num2str(rippeI(rIndV))]);
box off
subplot(3,3,4);
plot(1:numSteps,pI1_vec,'r.',1:numSteps,pI0_vec,'b.');
set(gca,'FontSize',14); title('Initial');
legend('Outbound','Inbound','Location','SouthEast');
xlim([0 numSteps]);ylim([0 1]);



%%
%%%%2. state different, initial and obs identical
%bond04, transition|I=1(center arm)
clear onestep_I1 onestep_I0 postx_I1 postx_I0 postxM_r;
C=size(spike_r,1);
%C=size(lambda_non,2);
dt=1/33.4;
n=length(stateV);numSteps=size(spike_r,2);
postx_I1=0.5*ones(n,1)./n;postx_I0=0.5*ones(n,1)./n; %random prior
pI1_vec=zeros(numSteps,1);pI0_vec=zeros(numSteps,1);
for t=1:numSteps
    onestep_I1=stateM_I1_gausnorm*postx_I1;
    onestep_I0=stateM_I0_gausnorm*postx_I0;
    L=ones(n,1);
    for c=1:C
        L=L.*((lambda_non(:,c).*dt).^spike_r(c,t).*exp(-lambda_non(:,c).*dt));
    end
    totnorm=sum(onestep_I1.*L)+sum(onestep_I0.*L);
    postx_I1=onestep_I1.*L./totnorm;
    postx_I0=onestep_I0.*L./totnorm;
    pI1_vec(t)=sum(postx_I1);
    pI0_vec(t)=sum(postx_I0);
    clear onestep_I1 onestep_I2
end

subplot(3,3,5);
plot(1:numSteps,pI1_vec,'r.',1:numSteps,pI0_vec,'b.');
set(gca,'FontSize',14);title('State transition');
%legend('Outbound','Inbound','Location','SouthEast');
xlim([0 numSteps]);ylim([0 1]);


%%
%%%%2.5 initial+state different; obs same
clear onestep_I1 onestep_I0 postx_I1 postx_I0 postxM_r;
C=size(spike_r,1);
dt=1/33.4;
n=length(stateV);numSteps=size(spike_r,2);
%P(x0|I)
Px_I1=exp(-stateV.^2./(2*0.2^2));Px_I1=Px_I1./sum(Px_I1);
Px_I0=max(Px_I1)*ones(1,n)-Px_I1; Px_I0=Px_I0./sum(Px_I0);
%P(x0)=P(x0|I)P(I);
postx_I1=0.5*Px_I1';postx_I0=0.5*Px_I0';
pI1_vec=zeros(numSteps,1);pI0_vec=zeros(numSteps,1);
for t=1:numSteps
    onestep_I1=stateM_I1_gausnorm*postx_I1;
    onestep_I0=stateM_I0_gausnorm*postx_I0;
    L=ones(n,1);
    L=ones(n,1);
    for c=1:C
        L=L.*((lambda_non(:,c).*dt).^spike_r(c,t).*exp(-lambda_non(:,c).*dt));
    end
    totnorm=sum(onestep_I1.*L)+sum(onestep_I0.*L);
    postx_I1=onestep_I1.*L./totnorm;
    postx_I0=onestep_I0.*L./totnorm;
    pI1_vec(t)=sum(postx_I1);
    pI0_vec(t)=sum(postx_I0);
    clear onestep_I1 onestep_I2
end

subplot(3,3,7)
plot(1:numSteps,pI1_vec,'r.',1:numSteps,pI0_vec,'b.');
set(gca,'FontSize',14); title('Initial & State transition');
%legend('Outbound','Inbound','Location','SouthEast');
xlim([0 numSteps]);ylim([0 1]);




%%
%%%%3. obs different, initial and state identical
%bond04, transition|I=1(center arm)

clear onestep_I1 onestep_I0 postx_I1 postx_I0 postxM_r;
C=size(spike_r,1);
dt=1/33.4;
n=length(stateV);numSteps=size(spike_r,2);
postx_I1=0.5*ones(n,1)./n;postx_I0=0.5*ones(n,1)./n; %random prior
pI1_vec=zeros(numSteps,1);pI0_vec=zeros(numSteps,1);
for t=1:numSteps
    onestep_I1=stateM_gausnorm*postx_I1;
    onestep_I0=stateM_gausnorm*postx_I0;
    L_I1=ones(n,1);L_I0=ones(n,1);
    for c=1:C
        L_I1=L_I1.*((lambda_non_I1(:,c).*dt).^spike_r(c,t).*exp(-lambda_non_I1(:,c).*dt));
        L_I0=L_I0.*((lambda_non_I0(:,c).*dt).^spike_r(c,t).*exp(-lambda_non_I0(:,c).*dt));
    end
    totnorm=sum(onestep_I1.*L_I1)+sum(onestep_I0.*L_I0);
    postx_I1=onestep_I1.*L_I1./totnorm;
    postx_I0=onestep_I0.*L_I0./totnorm;
    pI1_vec(t)=sum(postx_I1);
    pI0_vec(t)=sum(postx_I0);
    clear onestep_I1 onestep_I2
end

subplot(3,3,6);
plot(1:numSteps,pI1_vec,'r.',1:numSteps,pI0_vec,'b.');
set(gca,'FontSize',14);title('Observation');
%legend('Outbound','Inbound','Location','SouthEast');
xlim([0 numSteps]);
ylim([0 1]);


%%
%%%%3.5 all three

clear onestep_I1 onestep_I0 postx_I1 postx_I0 postxM_r;
C=size(spike_r,1);
dt=1/33.4;
n=length(stateV);numSteps=size(spike_r,2);
Px_I1=exp(-stateV.^2./(2*0.2^2));Px_I1=Px_I1./sum(Px_I1);
Px_I0=max(Px_I1)*ones(1,n)-Px_I1; Px_I0=Px_I0./sum(Px_I0);
postx_I1=0.5*Px_I1';postx_I0=0.5*Px_I0';
pI1_vec=zeros(numSteps,1);pI0_vec=zeros(numSteps,1);
for t=1:numSteps
    onestep_I1=stateM_I1_gausnorm*postx_I1;
    onestep_I0=stateM_I0_gausnorm*postx_I0;
    L_I1=ones(n,1);L_I0=ones(n,1);
    for c=1:C
        L_I1=L_I1.*((lambda_non_I1(:,c).*dt).^spike_r(c,t).*exp(-lambda_non_I1(:,c).*dt));
        L_I0=L_I0.*((lambda_non_I0(:,c).*dt).^spike_r(c,t).*exp(-lambda_non_I0(:,c).*dt));
    end
    totnorm=sum(onestep_I1.*L_I1)+sum(onestep_I0.*L_I0);
    postx_I1=onestep_I1.*L_I1./totnorm;
    postx_I0=onestep_I0.*L_I0./totnorm;
    pI1_vec(t)=sum(postx_I1);
    pI0_vec(t)=sum(postx_I0);
    clear onestep_I1 onestep_I2
end

subplot(3,3,9);
plot(1:numSteps,pI1_vec,'r.',1:numSteps,pI0_vec,'b.');
set(gca,'FontSize',14);
title('Initial+State+Observation');
%legend('Outbound','Inbound','Location','SouthEast');
xlim([0 numSteps]);
ylim([0 1]);

%%
%%%%3.55 state+obs

clear onestep_I1 onestep_I0 postx_I1 postx_I0 postxM_r;
C=size(spike_r,1);
dt=1/33.4;
n=length(stateV);numSteps=size(spike_r,2);
postx_I1=0.5*ones(n,1)./n;postx_I0=0.5*ones(n,1)./n; %random prior
pI1_vec=zeros(numSteps,1);pI0_vec=zeros(numSteps,1);
for t=1:numSteps
    onestep_I1=stateM_I1_gausnorm*postx_I1;
    onestep_I0=stateM_I0_gausnorm*postx_I0;
    L_I1=ones(n,1);L_I0=ones(n,1);
    for c=1:C
        L_I1=L_I1.*((lambda_non_I1(:,c).*dt).^spike_r(c,t).*exp(-lambda_non_I1(:,c).*dt));
        L_I0=L_I0.*((lambda_non_I0(:,c).*dt).^spike_r(c,t).*exp(-lambda_non_I0(:,c).*dt));
    end
    totnorm=sum(onestep_I1.*L_I1)+sum(onestep_I0.*L_I0);
    postx_I1=onestep_I1.*L_I1./totnorm;
    postx_I0=onestep_I0.*L_I0./totnorm;
    pI1_vec(t)=sum(postx_I1);
    pI0_vec(t)=sum(postx_I0);
    clear onestep_I1 onestep_I2
end

subplot(3,3,8);
plot(1:numSteps,pI1_vec,'r.',1:numSteps,pI0_vec,'b.');
set(gca,'FontSize',14);
title('State+Observation');
%legend('Outbound','Inbound','Location','SouthEast');
xlim([0 numSteps]);
ylim([0 1]);



%%
clear onstep postx postxM_r;
C=size(spike_r,1);
%C=size(lambda_non,2);
dt=1/33.4;
n=size(stateM_gausnorm,1);numSteps=size(spike_r,2);
postx=ones(n,1)./n; %random prior
postxM_r=zeros(n,numSteps);
tic;
for t=1:numSteps
    onestep=stateM_gausnorm*postx;
    L=ones(size(postx));
    for c=1:C
        L=L.*((lambda_non(:,c).*dt).^spike_r(c,t).*exp(-lambda_non(:,c).*dt));
    end
    postx=onestep.*L./sum(onestep.*L);
    postxM_r(:,t)=postx;
    clear onestep
end
toc



subplot(3,3,3);
imagesc(1:numSteps,stateV,postxM_r);
ylabel('Linearized position','FontSize',14); xlabel('Time (ms)','FontSize',14);
%set(gca,'FontSize',14,'YTick',[-6 -4 -2 0 2 4 6],'YTickLabel',[-6 -4 -2 0 2 4 6]);
set(gca,'FontSize',14,'YTick',[-3 -2 -1 0 1 2 3],'YTickLabel',[-3 -2 -1 0 1 2 3]);
colormap(flipud(hot(256)));
caxis([0 0.2]);
xlim([0 numSteps]);


