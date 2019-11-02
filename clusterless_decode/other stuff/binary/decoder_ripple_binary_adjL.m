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
spike=[];
betaMEAN_non=[];betaSD_non=[];
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

    [betahat,~,stats]=glmfit(X,sptrain0,'poisson','constant','off');
    [lambda,dylo,dyhi]=glmval(betahat,newX,'log',stats,'constant','off');
    
    betaSD=(log(lambda+dyhi)-log(lambda-dylo))./4;
    betaMEAN=log(lambda);
    
    %lambda_pause: place field for non-ripples for individual cell
    betaSD_non=[betaSD_non betaSD];
    betaMEAN_non=[betaMEAN_non betaMEAN];
    clear sptrain0 lambda_pause betaSD betaMEAN;
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
    
    lambda=[lambdaI1 lambdaI0];
    dylo=[dyloI1 dyloI0];
    dyhi=[dyhiI1 dyhiI0];
    
    lambda_cell{kk}=lambda;
    dylo_cell{kk}=dylo;
    dyhi_cell{kk}=dyhi;
    
    clear betahat lambda dylo dyhi
end
toc

C=size(a,2);
for c=1:C
    lambda=lambda_cell{c};
    dylo=dylo_cell{c};
    dyhi=dyhi_cell{c};
    
    betaSD_non_I1(:,c)=(log(lambda(:,1)+dyhi(:,1))-log(lambda(:,1)-dylo(:,1)))./4;
    betaSD_non_I0(:,c)=(log(lambda(:,2)+dyhi(:,2))-log(lambda(:,2)-dylo(:,2)))./4;

    betaMEAN_non_I1(:,c)=log(lambda(:,1));
    betaMEAN_non_I0(:,c)=log(lambda(:,2));
end



%%

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

minSD=12;
quantV=linspace(-1.3,1.3,100);
pdfV=exp(-quantV.^2)./sqrt(2*pi);




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
    L=ones(n,1);
    for c=1:C
        for i=1:n
            ave=betaMEAN_non(i,c);
            stdev=min(betaSD_non(i,c),minSD);
            lam=exp(ave+quantV.*stdev);
            L0a=(lam.*dt).^spike_r(c,t).*exp(-lam*dt);
            L0=L0a.*pdfV;
            L1=sum(L0);
            L2(i)=L1;
        end
        L3(c,:)=L2./sum(L2);
    end
    L=prod(L3,1)./sum(prod(L3,1));L=L';
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
        for i=1:n
            ave=betaMEAN_non(i,c);
            stdev=min(betaSD_non(i,c),minSD);
            lam=exp(ave+quantV.*stdev);
            L0a=(lam.*dt).^spike_r(c,t).*exp(-lam*dt);
            L0=L0a.*pdfV;
            L1=sum(L0);
            L2(i)=L1;
        end
        L3(c,:)=L2./sum(L2);
    end
    L=prod(L3,1)./sum(prod(L3,1));L=L';
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
    for c=1:C
        for i=1:n
            ave=betaMEAN_non(i,c);
            stdev=min(betaSD_non(i,c),minSD);
            lam=exp(ave+quantV.*stdev);
            L0a=(lam.*dt).^spike_r(c,t).*exp(-lam*dt);
            L0=L0a.*pdfV;
            L1=sum(L0);
            L2(i)=L1;
        end
        L3(c,:)=L2./sum(L2);
    end
    L=prod(L3,1)./sum(prod(L3,1));L=L';
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
        for i=1:n
            ave_I1=betaMEAN_non_I1(i,c);
            stdev_I1=min(betaSD_non_I1(i,c),minSD);
            lam_I1=exp(ave_I1+quantV.*stdev_I1);
            L0a_I1=(lam_I1.*dt).^spike_r(c,t).*exp(-lam_I1*dt);
            L0_I1=L0a_I1.*pdfV;
            L1_I1=sum(L0_I1);
            L2_I1(i)=L1_I1;
            
            ave_I0=betaMEAN_non_I0(i,c);
            stdev_I0=min(betaSD_non_I0(i,c),minSD);
            lam_I0=exp(ave_I0+quantV.*stdev_I0);
            L0a_I0=(lam_I0.*dt).^spike_r(c,t).*exp(-lam_I0*dt);
            L0_I0=L0a_I0.*pdfV;
            L1_I0=sum(L0_I0);
            L2_I0(i)=L1_I0;
        end
        L3_I1(c,:)=L2_I1./sum(L2_I1);
        L3_I0(c,:)=L2_I0./sum(L2_I0);
    end
    L_I1=prod(L3_I1,1)./sum(prod(L3_I1,1));L_I1=L_I1';
    L_I0=prod(L3_I0,1)./sum(prod(L3_I0,1));L_I0=L_I0';
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
        for i=1:n
            ave_I1=betaMEAN_non_I1(i,c);
            stdev_I1=min(betaSD_non_I1(i,c),minSD);
            lam_I1=exp(ave_I1+quantV.*stdev_I1);
            L0a_I1=(lam_I1.*dt).^spike_r(c,t).*exp(-lam_I1*dt);
            L0_I1=L0a_I1.*pdfV;
            L1_I1=sum(L0_I1);
            L2_I1(i)=L1_I1;
            
            ave_I0=betaMEAN_non_I0(i,c);
            stdev_I0=min(betaSD_non_I0(i,c),minSD);
            lam_I0=exp(ave_I0+quantV.*stdev_I0);
            L0a_I0=(lam_I0.*dt).^spike_r(c,t).*exp(-lam_I0*dt);
            L0_I0=L0a_I0.*pdfV;
            L1_I0=sum(L0_I0);
            L2_I0(i)=L1_I0;
        end
        L3_I1(c,:)=L2_I1./sum(L2_I1);
        L3_I0(c,:)=L2_I0./sum(L2_I0);
    end
    L_I1=prod(L3_I1,1)./sum(prod(L3_I1,1));L_I1=L_I1';
    L_I0=prod(L3_I0,1)./sum(prod(L3_I0,1));L_I0=L_I0';
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
        for i=1:n
            ave_I1=betaMEAN_non_I1(i,c);
            stdev_I1=min(betaSD_non_I1(i,c),minSD);
            lam_I1=exp(ave_I1+quantV.*stdev_I1);
            L0a_I1=(lam_I1.*dt).^spike_r(c,t).*exp(-lam_I1*dt);
            L0_I1=L0a_I1.*pdfV;
            L1_I1=sum(L0_I1);
            L2_I1(i)=L1_I1;
            
            ave_I0=betaMEAN_non_I0(i,c);
            stdev_I0=min(betaSD_non_I0(i,c),minSD);
            lam_I0=exp(ave_I0+quantV.*stdev_I0);
            L0a_I0=(lam_I0.*dt).^spike_r(c,t).*exp(-lam_I0*dt);
            L0_I0=L0a_I0.*pdfV;
            L1_I0=sum(L0_I0);
            L2_I0(i)=L1_I0;
        end
        L3_I1(c,:)=L2_I1./sum(L2_I1);
        L3_I0(c,:)=L2_I0./sum(L2_I0);
    end
    L_I1=prod(L3_I1,1)./sum(prod(L3_I1,1));L_I1=L_I1';
    L_I0=prod(L3_I0,1)./sum(prod(L3_I0,1));L_I0=L_I0';
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
        for i=1:n
            ave=betaMEAN_non(i,c);
            stdev=min(betaSD_non(i,c),minSD);
            lam=exp(ave+quantV.*stdev);
            L0a=(lam.*dt).^spike_r(c,t).*exp(-lam*dt);
            L0=L0a.*pdfV;
            L1=sum(L0);
            L2(i)=L1;
        end
        L3(c,:)=L2./sum(L2);
    end
    L=prod(L3,1)./sum(prod(L3,1));L=L';
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


