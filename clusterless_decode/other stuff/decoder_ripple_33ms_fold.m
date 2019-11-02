%%% pre-processing folding %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% vecL : "folded" version of vecL0, eliminating directionality
    % while maintaining L vs. R
vecL = [];
for i=1:length(vecL0)
    if vecL0(i)>=3
        vecL(i)=6-vecL0(i);
    elseif vecL0(i)<=-3
        vecL(i)=-6-vecL0(i);
    else
        vecL(i)=vecL0(i);
    end
end

% vecLF : col 1: ascending time; col 2: sorted linearized folded position
vecLF0 = [vecT vecL'];
[d1,d2] = sort(vecLF0(:,1));  
vecLF = vecLF0(d2,:); % 

A=pos{ex}{ep}.data(:,1);                    % time stamps for animal's trajectory
ti=round(A(1)*1000):1:round(A(end)*1000);   % t: 1-ms binned time stamps


%%
%%%% encoder 33 ms %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X=[];
s=0.5;
c_pt_lin=linspace(-3, 3.0001, 15);     % c_pt_lin: [-3 +3] in 15 spatial bins
%13 bins for frank, bond03; 15 for bond04(26,fold1);
c_pt_lin_all=[-2.8 c_pt_lin 2.8];
num_c_pts = length(c_pt_lin_all);
vecInd=find(vecLF(:,1)~=0);

% iterate through vecLF (i.e. each sequential linearized folded position)
X = [];   
for k0=1:length(vecInd)
    
    % identify the linear bin
    k=vecInd(k0);   
    nearest_c_pt_index = find(c_pt_lin<=vecLF(k,2), 1, 'last')+1;
    nearest_c_pt_lin = c_pt_lin_all(nearest_c_pt_index);
    next_c_pt_lin = c_pt_lin_all(nearest_c_pt_index+1);
    
    % u is the distance to the nearest linear bin divided by the bin size
    u = (vecLF(k,2)-nearest_c_pt_lin)/(next_c_pt_lin-nearest_c_pt_lin);
    p=[u^3 u^2 u 1]*[-s 2-s s-2 s;2*s s-3 3-2*s -s;-s 0 s 0;0 1 0 0];
    X(k0,:) = [  zeros(1,nearest_c_pt_index-2)     p      zeros(1,num_c_pts-4-(nearest_c_pt_index-2))  ];
end

newX=[];
stateV=linspace(-3,3,121);           % folded linear bins, 121 bins
%stateV=[linspace(-3,0,61) linspace(1,3,40)];
num_c_pts = length(c_pt_lin_all);
for k=1:size(stateV,2)
    nearest_c_pt_index = find(c_pt_lin<=stateV(k), 1, 'last')+1;
    nearest_c_pt_lin = c_pt_lin_all(nearest_c_pt_index);
    next_c_pt_lin = c_pt_lin_all(nearest_c_pt_index+1);
    u = (stateV(k)-nearest_c_pt_lin)/(next_c_pt_lin-nearest_c_pt_lin);
    p=[u^3 u^2 u 1]*[-s 2-s s-2 s;2*s s-3 3-2*s -s;-s 0 s 0;0 1 0 0];
    newX(k,:) = [zeros(1,nearest_c_pt_index-2) p zeros(1,num_c_pts-4-(nearest_c_pt_index-2))];
end

if 1
   figure
   imagesc(X);
    figure
   imagesc(newX);
end


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
%%%%decoder; full trajectory %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear onstep postx postxM_r;
tic;
C=length(a);dt=1;n=size(stateM_gausnorm,1);T=size(vecLF,1);
postx=ones(n,1)./n; %random prior
for t=1:T
    onestep=stateM_gausnorm*postx; %statematrix is n*n; onestep is n*1;
    L=ones(size(postx));
    for c=1:C %per cell
        L=L.*((lambda_non(:,c)*dt).^spike(c,t).*exp(-lambda_non(:,c).*dt));
            %spike is C*numSteps; lambda_non is n*C
    end
    L=L./sum(L);
    postx=onestep.*L./sum(onestep.*L);
    postxM(:,t)=postx;
end
toc

figure
imagesc(vecLF(:,1),stateV,postxM);
set(gca,'FontSize',12);
%title('Posterior Density for Non-ripples','FontSize',18); %forward empirical state transition matrix
ylabel('Linearized position','FontSize',18); xlabel('Time (sec)','FontSize',18);
colormap(flipud(hot(256)));
c=colorbar('location','eastoutside');
%caxis([0 1]);
caxis([0 0.3]);
ylabel(c,'Probability density','FontSize',18);
hold on
plot(vecLF(:,1),vecLF(:,2),'ks','MarkerSize',2,'MarkerFaceColor','k','MarkerEdgeColor','none');
hold off



%%





%%
%%%%ripples
%%%%%%Create indicator function of ripples for 1ms binned spike trains
ripple_all=[];
ind_ripple=zeros(size(ti,2),size(ind_t,2));

for ii=1:size(ind_t,2)
    k=ind_t(ii);
    ripple_s=round(ripples{ex}{ep}{k}.starttime*1000);
    ripple_e=round(ripples{ex}{ep}{k}.endtime*1000);
    for l=1:size(ripple_s,1)
        ripple=[ripple_s(l) ripple_e(l)];
        ripple_all=[ripple_all;ripple];
    end
    clear ripple_s ripple_e ripple;
    for s=1:size(ripple_all,1)
        ind_ripple(find(ti>=ripple_all(s,1) & ti<=ripple_all(s,2)),ii) = 1;
    end
    ripple_all=[];
end

ripple_acr=zeros(size(ind_ripple,1),1);
%ripple_acr=sum(ind_ripple,2);
ripple_acr(find(sum(ind_ripple,2)))=1;
%%%indicator for ripples pooled across tetrodes

%ripple_acr(find(sum(ind_ripple,2)>=6))=1;
%%%indicator for ripples pooled across tetrodes where more than 6/9
%%%tetrodes have spiking



clear startT endT
for k=1:size(ti,2)-1
    if ripple_acr(k)==0 && ripple_acr(k+1)==1
        startT(k+1)=k+1;
    elseif ripple_acr(k)==1 && ripple_acr(k+1)==0
        endT(k)=k;
    end
end
if ripple_acr(size(ti,2))==1
    endT(size(ti,2))=size(ti,2);
end

ripple_seg=[find(startT)' find(endT)']; %index for ripple segments


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
    sumC(k)=length(find(sum(spike_r,2)));
    sumR(k)=sum(spike_r(:));
end

kk=find(sumR>25);
length(kk)


%%

%now need to specify which ripple to decode, for example ripple #41:
len_ripple=ripple_seg(:,2)-ripple_seg(:,1);

rInd=2;
%spike_r=[spike_r_all{kk(123)} spike_r_all{kk(124)}];
%spike_r=[spike_r_all{910} spike_r_all{911}];
spike_r=spike_r_all{kk(rInd)};

clear onstep postx postxM_r;
%C=size(spike_r,1);
C=size(lambda_non,2);
dt=1/33.4;
n=size(stateM_gausnorm,1);numSteps=size(spike_r,2);
postx=ones(n,1)./n; %random prior
%postx=[1/3;zeros(59,1);1/3;zeros(59,1);1/3]; %delta prior
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


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%      Figure      %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
imagesc(1:numSteps,stateV,postxM_r);
rloc_Ind=find(A*1000>ti(ripple_seg(kk(rInd),1))&A*1000<ti(ripple_seg(kk(rInd),2)));
set(gca,'XTick',1:50:numSteps,'XTickLabel',vecLF(rloc_Ind(1)-1,1):0.05:vecLF(rloc_Ind(end),1),'FontSize',14);
%title('Posterior Density, forward empirical state transition matrix','FontSize',10);
ylabel('Linearized position','FontSize',20); xlabel('Time (ms)','FontSize',20);
colormap(flipud(hot(256)));
c=colorbar('location','eastoutside');
caxis([0 0.3]);ylabel(c,'Probability density','FontSize',14);
%title('r1240 Sorted','FontSize',20);

figure;
imagesc(spike_r);colormap(flipud(gray));freezeColors
set(gca,'XTick',1:50:numSteps,'XTickLabel',vecLF(rloc_Ind(1)-1,1):0.05:vecLF(rloc_Ind(end),1),'FontSize',14);
xlabel('Time (ms)','FontSize',20);
ylabel('Cell index','FontSize',20);



