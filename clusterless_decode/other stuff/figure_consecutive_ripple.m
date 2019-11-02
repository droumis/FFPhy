%%%four ripples
figure;
subplot(3,4,1:4);
rloc_min=min(find(A*1000>ti(ripple_seg(kk(rIndV(1)),1))&A*1000<ti(ripple_seg(kk(rIndV(1)),2))))-5;
rloc_max=max(find(A*1000>ti(ripple_seg(kk(rIndV(end)),1))&A*1000<ti(ripple_seg(kk(rIndV(end)),2))))+5;
rloc_Ind=linspace(rloc_min,rloc_max,rloc_max-rloc_min+1);
plot(vecLF(rloc_Ind,1),vecLF(rloc_Ind,2),'.','MarkerSize',5);
ylim([-3.2 3.2]);
set(gca,'YDir','reverse','FontSize',14,'YTick',[-3 -1.5 0 1.5 3]);
xlabel 'Time (sec)';
hold on
for rInd=rIndV
spike_r=spike_r_all{kk(rInd)};
rloc_Ind=find(A*1000>ti(ripple_seg(kk(rInd),1))&A*1000<ti(ripple_seg(kk(rInd),2)));
plot(vecLF([rloc_Ind(1)-1;rloc_Ind],1),vecLF([rloc_Ind(1)-1;rloc_Ind],2),'.r','MarkerSize',5);
end
hold off

subplot(3,4,5)
spike_r=spike_r_all{kk(rIndV(1))};
imagesc(spike_r);colormap(flipud(gray));freezeColors
title(['ripple = ',num2str(kk(rIndV(1)))]);
set(gca,'YTick',1:13);

subplot(3,4,6)
spike_r=spike_r_all{kk(rIndV(2))};
imagesc(spike_r);colormap(flipud(gray));freezeColors
title(['ripple = ',num2str(kk(rIndV(2)))]);
set(gca,'YTick',1:13);

subplot(3,4,7)
spike_r=spike_r_all{kk(rIndV(3))};
imagesc(spike_r);colormap(flipud(gray));freezeColors
title(['ripple = ',num2str(kk(rIndV(3)))]);
set(gca,'YTick',1:13);

subplot(3,4,8)
spike_r=spike_r_all{kk(rIndV(4))};
imagesc(spike_r);colormap(flipud(gray));freezeColors
title(['ripple = ',num2str(kk(rIndV(4)))]);
set(gca,'YTick',1:13);

subplot(3,4,9)
spike_r=spike_r_all{kk(rIndV(1))};
clear onstep postx postxM_r;
C=size(lambda_non,2);
dt=1/33.4;
n=size(stateM_gausnorm,1);numSteps=size(spike_r,2);
postx=ones(n,1)./n; %random prior
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
stateV=linspace(-3,3,121);
imagesc(1:numSteps,stateV,postxM_r);
ylabel('Linearized position','FontSize',14); xlabel('Time (ms)','FontSize',14);
colormap(flipud(hot(256)));
caxis([0 0.1]);

subplot(3,4,10)
spike_r=spike_r_all{kk(rIndV(2))};
clear onstep postx postxM_r;
C=size(lambda_non,2);
dt=1/33.4;
n=size(stateM_gausnorm,1);numSteps=size(spike_r,2);
postx=ones(n,1)./n; %random prior
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
stateV=linspace(-3,3,121);
imagesc(1:numSteps,stateV,postxM_r);
ylabel('Linearized position','FontSize',14); xlabel('Time (ms)','FontSize',14);
colormap(flipud(hot(256)));
caxis([0 0.1]);

subplot(3,4,11)
spike_r=spike_r_all{kk(rIndV(3))};
clear onstep postx postxM_r;
C=size(lambda_non,2);
dt=1/33.4;
n=size(stateM_gausnorm,1);numSteps=size(spike_r,2);
postx=ones(n,1)./n; %random prior
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
stateV=linspace(-3,3,121);
imagesc(1:numSteps,stateV,postxM_r);
ylabel('Linearized position','FontSize',14); xlabel('Time (ms)','FontSize',14);
colormap(flipud(hot(256)));
caxis([0 0.1]);

subplot(3,4,12)
spike_r=spike_r_all{kk(rIndV(4))};
clear onstep postx postxM_r;
C=size(lambda_non,2);
dt=1/33.4;
n=size(stateM_gausnorm,1);numSteps=size(spike_r,2);
postx=ones(n,1)./n; %random prior
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
stateV=linspace(-3,3,121);
imagesc(1:numSteps,stateV,postxM_r);
ylabel('Linearized position','FontSize',14); xlabel('Time (ms)','FontSize',14);
colormap(flipud(hot(256)));
caxis([0 0.1]);





%%
%%%a pair
figure;
subplot(3,2,1:2);
rloc_min=min(find(A*1000>ti(ripple_seg(kk(rIndV(1)),1))&A*1000<ti(ripple_seg(kk(rIndV(1)),2))))-5;
rloc_max=max(find(A*1000>ti(ripple_seg(kk(rIndV(end)),1))&A*1000<ti(ripple_seg(kk(rIndV(end)),2))))+5;
rloc_Ind=linspace(rloc_min,rloc_max,rloc_max-rloc_min+1);
plot(vecLF(rloc_Ind,1),vecLF(rloc_Ind,2),'.','MarkerSize',5);
ylim([-3.2 3.2]);
set(gca,'YDir','reverse','FontSize',14,'YTick',[-3 -1.5 0 1.5 3]);
xlabel 'Time (sec)';
hold on
for rInd=rIndV
spike_r=spike_r_all{kk(rInd)};
rloc_Ind=find(A*1000>ti(ripple_seg(kk(rInd),1))&A*1000<ti(ripple_seg(kk(rInd),2)));
plot(vecLF([rloc_Ind(1)-1;rloc_Ind],1),vecLF([rloc_Ind(1)-1;rloc_Ind],2),'.r','MarkerSize',5);
end
hold off

subplot(3,2,3)
spike_r=spike_r_all{kk(rIndV(1))};
imagesc(spike_r);colormap(flipud(gray));freezeColors
title(['ripple = ',num2str(kk(rIndV(1)))]);
set(gca,'YTick',1:13);

subplot(3,2,4)
spike_r=spike_r_all{kk(rIndV(2))};
imagesc(spike_r);colormap(flipud(gray));freezeColors
title(['ripple = ',num2str(kk(rIndV(2)))]);
set(gca,'YTick',1:13);


subplot(3,2,5)
spike_r=spike_r_all{kk(rIndV(1))};
clear onstep postx postxM_r;
C=size(lambda_non,2);
dt=1/33.4;
n=size(stateM_gausnorm,1);numSteps=size(spike_r,2);
postx=ones(n,1)./n; %random prior
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
stateV=linspace(-3,3,121);
imagesc(1:numSteps,stateV,postxM_r);
ylabel('Linearized position','FontSize',14); xlabel('Time (ms)','FontSize',14);
colormap(flipud(hot(256)));
caxis([0 0.1]);

subplot(3,2,6)
spike_r=spike_r_all{kk(rIndV(2))};
clear onstep postx postxM_r;
C=size(lambda_non,2);
dt=1/33.4;
n=size(stateM_gausnorm,1);numSteps=size(spike_r,2);
postx=ones(n,1)./n; %random prior
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
stateV=linspace(-3,3,121);
imagesc(1:numSteps,stateV,postxM_r);
ylabel('Linearized position','FontSize',14); xlabel('Time (ms)','FontSize',14);
colormap(flipud(hot(256)));
caxis([0 0.1]);








%%
%%%sorted spiking
% subplot(3,4,5)
% spike_r=spike_r_all{kk(rIndV(1))};
% a0=zeros(size(spike_r,1),1);
% for i=1:size(spike_r,1)
%     if isempty(find(spike_r(i,:)))==0
%     a0(i)=min(find(spike_r(i,:)));
%     end
% end
% [b0,c0]=sort(a0,'descend');
% imagesc(spike_r(c0,:));colormap(flipud(gray));freezeColors
% title(['ripple = ',num2str(kk(rIndV(1)))]);
% set(gca,'YTickLabel',c0,'YTick',1:13);
% 
% subplot(3,4,6)
% spike_r=spike_r_all{kk(rIndV(2))};
% a0=zeros(size(spike_r,1),1);
% for i=1:size(spike_r,1)
%     if isempty(find(spike_r(i,:)))==0
%     a0(i)=min(find(spike_r(i,:)));
%     end
% end
% [b0,c0]=sort(a0,'descend');
% imagesc(spike_r(c0,:));colormap(flipud(gray));freezeColors
% title(['ripple = ',num2str(kk(rIndV(2)))]);
% set(gca,'YTickLabel',c0,'YTick',1:13);
% 
% subplot(3,4,7)
% spike_r=spike_r_all{kk(rIndV(3))};
% a0=zeros(size(spike_r,1),1);
% for i=1:size(spike_r,1)
%     if isempty(find(spike_r(i,:)))==0
%     a0(i)=min(find(spike_r(i,:)));
%     end
% end
% [b0,c0]=sort(a0,'descend');
% imagesc(spike_r(c0,:));colormap(flipud(gray));freezeColors
% title(['ripple = ',num2str(kk(rIndV(3)))]);
% set(gca,'YTickLabel',c0,'YTick',1:13);
% 
% subplot(3,4,8)
% spike_r=spike_r_all{kk(rIndV(4))};
% a0=zeros(size(spike_r,1),1);
% for i=1:size(spike_r,1)
%     if isempty(find(spike_r(i,:)))==0
%     a0(i)=min(find(spike_r(i,:)));
%     end
% end
% [b0,c0]=sort(a0,'descend');
% imagesc(spike_r(c0,:));colormap(flipud(gray));freezeColors
% title(['ripple = ',num2str(kk(rIndV(4)))]);
% set(gca,'YTickLabel',c0,'YTick',1:13);
%%





%%









%%
%%%%unfolded
figure;
subplot(3,4,1:4);
rloc_min=min(find(A*1000>ti(ripple_seg(kk(rIndV(1)),1))&A*1000<ti(ripple_seg(kk(rIndV(1)),2))))-5;
rloc_max=max(find(A*1000>ti(ripple_seg(kk(rIndV(end)),1))&A*1000<ti(ripple_seg(kk(rIndV(end)),2))))+5;
rloc_Ind=linspace(rloc_min,rloc_max,rloc_max-rloc_min+1);
plot(vecLF(rloc_Ind,1),vecLF(rloc_Ind,2),'.','MarkerSize',5);
ylim([-6.2 6.2]);
set(gca,'YDir','reverse','FontSize',14,'YTick',[-6 -3 0 3 6]);
xlabel 'Time (sec)';
hold on
for rInd=rIndV
spike_r=spike_r_all{kk(rInd)};
rloc_Ind=find(A*1000>ti(ripple_seg(kk(rInd),1))&A*1000<ti(ripple_seg(kk(rInd),2)));
plot(vecLF([rloc_Ind(1)-1;rloc_Ind],1),vecLF([rloc_Ind(1)-1;rloc_Ind],2),'.r','MarkerSize',5);
end
hold off

subplot(3,4,5)
spike_r=spike_r_all{kk(rIndV(1))};
imagesc(spike_r);colormap(flipud(gray));freezeColors
title(['ripple = ',num2str(kk(rIndV(1)))]);
set(gca,'YTick',1:13);

subplot(3,4,6)
spike_r=spike_r_all{kk(rIndV(2))};
imagesc(spike_r);colormap(flipud(gray));freezeColors
title(['ripple = ',num2str(kk(rIndV(2)))]);
set(gca,'YTick',1:13);

subplot(3,4,7)
spike_r=spike_r_all{kk(rIndV(3))};
imagesc(spike_r);colormap(flipud(gray));freezeColors
title(['ripple = ',num2str(kk(rIndV(3)))]);
set(gca,'YTick',1:13);

subplot(3,4,8)
spike_r=spike_r_all{kk(rIndV(4))};
imagesc(spike_r);colormap(flipud(gray));freezeColors
title(['ripple = ',num2str(kk(rIndV(4)))]);
set(gca,'YTick',1:13);



subplot(3,4,9)
spike_r=spike_r_all{kk(rIndV(1))};
clear onstep postx postxM_r;
C=size(lambda_non,2);
dt=1/33.4;
n=size(stateM_gausnorm,1);numSteps=size(spike_r,2);
postx=ones(n,1)./n; %random prior
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
imagesc(1:numSteps,stateV,postxM_r);
set(gca,'YTick',[-6 -3 0 3 6]);
ylabel('Linearized position','FontSize',14); xlabel('Time (ms)','FontSize',14);
colormap(flipud(hot(256)));
caxis([0 0.1]);

subplot(3,4,10)
spike_r=spike_r_all{kk(rIndV(2))};
clear onstep postx postxM_r;
C=size(lambda_non,2);
n=size(stateM_gausnorm,1);numSteps=size(spike_r,2);
postx=ones(n,1)./n; %random prior
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
imagesc(1:numSteps,stateV,postxM_r);
ylabel('Linearized position','FontSize',14); xlabel('Time (ms)','FontSize',14);
colormap(flipud(hot(256)));
caxis([0 0.1]);
set(gca,'YTick',[-6 -3 0 3 6]);

subplot(3,4,11)
spike_r=spike_r_all{kk(rIndV(3))};
clear onstep postx postxM_r;
C=size(lambda_non,2);
n=size(stateM_gausnorm,1);numSteps=size(spike_r,2);
postx=ones(n,1)./n; %random prior
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
imagesc(1:numSteps,stateV,postxM_r);
ylabel('Linearized position','FontSize',14); xlabel('Time (ms)','FontSize',14);
colormap(flipud(hot(256)));
caxis([0 0.1]);
set(gca,'YTick',[-6 -3 0 3 6]);

subplot(3,4,12)
spike_r=spike_r_all{kk(rIndV(4))};
clear onstep postx postxM_r;
C=size(lambda_non,2);
n=size(stateM_gausnorm,1);numSteps=size(spike_r,2);
postx=ones(n,1)./n; %random prior
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
imagesc(1:numSteps,stateV,postxM_r);
ylabel('Linearized position','FontSize',14); xlabel('Time (ms)','FontSize',14);
colormap(flipud(hot(256)));
caxis([0 0.1]);
set(gca,'YTick',[-6 -3 0 3 6]);






%%
%%%a pair
figure;
subplot(3,2,1:2);
rloc_min=min(find(A*1000>ti(ripple_seg(kk(rIndV(1)),1))&A*1000<ti(ripple_seg(kk(rIndV(1)),2))))-5;
rloc_max=max(find(A*1000>ti(ripple_seg(kk(rIndV(end)),1))&A*1000<ti(ripple_seg(kk(rIndV(end)),2))))+5;
rloc_Ind=linspace(rloc_min,rloc_max,rloc_max-rloc_min+1);
plot(vecLF(rloc_Ind,1),vecLF(rloc_Ind,2),'.','MarkerSize',5);
ylim([-6.2 6.2]);
set(gca,'YDir','reverse','FontSize',14,'YTick',[-6 -3 0 3 6]);
xlabel 'Time (sec)';
hold on
for rInd=rIndV
spike_r=spike_r_all{kk(rInd)};
rloc_Ind=find(A*1000>ti(ripple_seg(kk(rInd),1))&A*1000<ti(ripple_seg(kk(rInd),2)));
plot(vecLF([rloc_Ind(1)-1;rloc_Ind],1),vecLF([rloc_Ind(1)-1;rloc_Ind],2),'.r','MarkerSize',5);
end
hold off

subplot(3,2,3)
spike_r=spike_r_all{kk(rIndV(1))};
imagesc(spike_r);colormap(flipud(gray));freezeColors
title(['ripple = ',num2str(kk(rIndV(1)))]);
set(gca,'YTick',1:13);

subplot(3,2,4)
spike_r=spike_r_all{kk(rIndV(2))};
imagesc(spike_r);colormap(flipud(gray));freezeColors
title(['ripple = ',num2str(kk(rIndV(2)))]);
set(gca,'YTick',1:13);


subplot(3,2,5)
spike_r=spike_r_all{kk(rIndV(1))};
clear onstep postx postxM_r;
C=size(lambda_non,2);
dt=1/33.4;
n=size(stateM_gausnorm,1);numSteps=size(spike_r,2);
postx=ones(n,1)./n; %random prior
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
imagesc(1:numSteps,stateV,postxM_r);
ylabel('Linearized position','FontSize',14); xlabel('Time (ms)','FontSize',14);
colormap(flipud(hot(256)));
caxis([0 0.1]);

subplot(3,2,6)
spike_r=spike_r_all{kk(rIndV(2))};
clear onstep postx postxM_r;
C=size(lambda_non,2);
dt=1/33.4;
n=size(stateM_gausnorm,1);numSteps=size(spike_r,2);
postx=ones(n,1)./n; %random prior
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
imagesc(1:numSteps,stateV,postxM_r);
ylabel('Linearized position','FontSize',14); xlabel('Time (ms)','FontSize',14);
colormap(flipud(hot(256)));
caxis([0 0.1]);

