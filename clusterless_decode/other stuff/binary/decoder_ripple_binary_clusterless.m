%%
%%%%run decoder_ripple_binary.m first
%%%%clusterless

spike_tim=ripple_seg(rippleI(rIndV),1):ripple_seg(rippleI(rIndV),2); %from 1 to 90000~
numSteps=length(spike_tim);


%%
xi=round(time/10);

%%
clear onestep_I1 onestep_I0 postx_I1 postx_I0 postxM_r;
C=size(spike_r,1);
dt=1/33.4;
n=length(stateV);
%numSteps=size(spike_r,2);
%P(x0|I);
Px_I1=exp(-stateV.^2./(2*0.2^2));Px_I1=Px_I1./sum(Px_I1);
Px_I0=max(Px_I1)*ones(1,n)-Px_I1; Px_I0=Px_I0./sum(Px_I0);
%P(x0)=P(x0|I)P(I);
postx_I1=0.5*Px_I1';postx_I0=0.5*Px_I0';
pI1_vec=zeros(numSteps,1);pI0_vec=zeros(numSteps,1);
postxM_r_I1=zeros(n,numSteps);postxM_r_I0=zeros(n,numSteps);

tic;
for t=1:numSteps
    tt=spike_tim(t);
    aa=find(xi==ti(tt));
    
    onestep_I1=stateM_I1_gausnorm*postx_I1;
    onestep_I0=stateM_I0_gausnorm*postx_I0;
    L=ones(n,1);

    if isempty(aa)==1 %if no spike occurs at time t
        L=exp(-Lint.*dt);
    elseif isempty(aa)==0 %if spikes
        l=zeros(n,length(aa));
        for j=1:length(aa)
            jj=aa(j);
            tetVec=tet_ind(jj,:);
            
            if tetVec(1)==1 %tet1
                spike_r(1,t)=1;
                i=tet_sum(jj,1);
                l0=normpdf(markAll_t1(i,2)*ones(1,length(spikeT0_t1)),markAll_t1(:,2)',smker).*normpdf(markAll_t1(i,3)*ones(1,length(spikeT0_t1)),markAll_t1(:,3)',smker).*normpdf(markAll_t1(i,4)*ones(1,length(spikeT0_t1)),markAll_t1(:,4)',smker).*normpdf(markAll_t1(i,5)*ones(1,length(spikeT0_t1)),markAll_t1(:,5)',smker);
                l1=Xnum_t1*l0'./occ(:,1)./dt;
                l2=l1.*dt.*exp(-Lint_t1.*dt);
                l2=l2./sum(l2);
                l(:,j)=l2;
            elseif tetVec(2)==1 %tet2
                spike_r(2,t)=1;
                i=tet_sum(jj,2);
                l0=normpdf(markAll_t2(i,2)*ones(1,length(spikeT0_t2)),markAll_t2(:,2)',smker).*normpdf(markAll_t2(i,3)*ones(1,length(spikeT0_t2)),markAll_t2(:,3)',smker).*normpdf(markAll_t2(i,4)*ones(1,length(spikeT0_t2)),markAll_t2(:,4)',smker).*normpdf(markAll_t2(i,5)*ones(1,length(spikeT0_t2)),markAll_t2(:,5)',smker);
                l1=Xnum_t2*l0'./occ(:,1)./dt;
                l2=l1.*dt.*exp(-Lint_t2.*dt);
                l2=l2./sum(l2);
                l(:,j)=l2;
            elseif tetVec(3)==1 %tet4
                spike_r(3,t)=1;
                i=tet_sum(jj,3);
                l0=normpdf(markAll_t4(i,2)*ones(1,length(spikeT0_t4)),markAll_t4(:,2)',smker).*normpdf(markAll_t4(i,3)*ones(1,length(spikeT0_t4)),markAll_t4(:,3)',smker).*normpdf(markAll_t4(i,4)*ones(1,length(spikeT0_t4)),markAll_t4(:,4)',smker).*normpdf(markAll_t4(i,5)*ones(1,length(spikeT0_t4)),markAll_t4(:,5)',smker);
                l1=Xnum_t4*l0'./occ(:,1)./dt;
                l2=l1.*dt.*exp(-Lint_t4.*dt);
                l2=l2./sum(l2);
                l(:,j)=l2;
            elseif tetVec(4)==1 %tet5
                spike_r(4,t)=1;
                i=tet_sum(jj,4);
                l0=normpdf(markAll_t5(i,2)*ones(1,length(spikeT0_t5)),markAll_t5(:,2)',smker).*normpdf(markAll_t5(i,3)*ones(1,length(spikeT0_t5)),markAll_t5(:,3)',smker).*normpdf(markAll_t5(i,4)*ones(1,length(spikeT0_t5)),markAll_t5(:,4)',smker).*normpdf(markAll_t5(i,5)*ones(1,length(spikeT0_t5)),markAll_t5(:,5)',smker);
                l1=Xnum_t5*l0'./occ(:,1)./dt;
                l2=l1.*dt.*exp(-Lint_t5.*dt);
                l2=l2./sum(l2);
                l(:,j)=l2;
            elseif tetVec(5)==1 %tet7
                spike_r(5,t)=1;
                i=tet_sum(jj,5);
                l0=normpdf(markAll_t7(i,2)*ones(1,length(spikeT0_t7)),markAll_t7(:,2)',smker).*normpdf(markAll_t7(i,3)*ones(1,length(spikeT0_t7)),markAll_t7(:,3)',smker).*normpdf(markAll_t7(i,4)*ones(1,length(spikeT0_t7)),markAll_t7(:,4)',smker).*normpdf(markAll_t7(i,5)*ones(1,length(spikeT0_t7)),markAll_t7(:,5)',smker);
                l1=Xnum_t7*l0'./occ(:,1)./dt;
                l2=l1.*dt.*exp(-Lint_t7.*dt);
                l2=l2./sum(l2);
                l(:,j)=l2;
            elseif tetVec(6)==1 %tet10
                spike_r(6,t)=1;
                i=tet_sum(jj,6);
                l0=normpdf(markAll_t10(i,2)*ones(1,length(spikeT0_t10)),markAll_t10(:,2)',smker).*normpdf(markAll_t10(i,3)*ones(1,length(spikeT0_t10)),markAll_t10(:,3)',smker).*normpdf(markAll_t10(i,4)*ones(1,length(spikeT0_t10)),markAll_t10(:,4)',smker).*normpdf(markAll_t10(i,5)*ones(1,length(spikeT0_t10)),markAll_t10(:,5)',smker);
                l1=Xnum_t10*l0'./occ(:,1)./dt;
                l2=l1.*dt.*exp(-Lint_t10.*dt);
                l2=l2./sum(l2);
                l(:,j)=l2;
            elseif tetVec(7)==1 %tet11
                spike_r(7,t)=1;
                i=tet_sum(jj,7);
                l0=normpdf(markAll_t11(i,2)*ones(1,length(spikeT0_t11)),markAll_t11(:,2)',smker).*normpdf(markAll_t11(i,3)*ones(1,length(spikeT0_t11)),markAll_t11(:,3)',smker).*normpdf(markAll_t11(i,4)*ones(1,length(spikeT0_t11)),markAll_t11(:,4)',smker).*normpdf(markAll_t11(i,5)*ones(1,length(spikeT0_t11)),markAll_t11(:,5)',smker);
                l1=Xnum_t11*l0'./occ(:,1)./dt;
                l2=l1.*dt.*exp(-Lint_t11.*dt);
                l2=l2./sum(l2);
                l(:,j)=l2;
            elseif tetVec(8)==1 %tet12
                spike_r(8,t)=1;
                i=tet_sum(jj,8);
                l0=normpdf(markAll_t12(i,2)*ones(1,length(spikeT0_t12)),markAll_t12(:,2)',smker).*normpdf(markAll_t12(i,3)*ones(1,length(spikeT0_t12)),markAll_t12(:,3)',smker).*normpdf(markAll_t12(i,4)*ones(1,length(spikeT0_t12)),markAll_t12(:,4)',smker).*normpdf(markAll_t12(i,5)*ones(1,length(spikeT0_t12)),markAll_t12(:,5)',smker);
                l1=Xnum_t12*l0'./occ(:,1)./dt;
                l2=l1.*dt.*exp(-Lint_t12.*dt);
                l2=l2./sum(l2);
                l(:,j)=l2;
            elseif tetVec(9)==1 %tet13
                spike_r(9,t)=1;
                i=tet_sum(jj,9);
                l0=normpdf(markAll_t13(i,2)*ones(1,length(spikeT0_t13)),markAll_t13(:,2)',smker).*normpdf(markAll_t13(i,3)*ones(1,length(spikeT0_t13)),markAll_t13(:,3)',smker).*normpdf(markAll_t13(i,4)*ones(1,length(spikeT0_t13)),markAll_t13(:,4)',smker).*normpdf(markAll_t13(i,5)*ones(1,length(spikeT0_t13)),markAll_t13(:,5)',smker);
                l1=Xnum_t13*l0'./occ(:,1)./dt;
                l2=l1.*dt.*exp(-Lint_t13.*dt);
                l2=l2./sum(l2);
                l(:,j)=l2;
            elseif tetVec(10)==1 %tet14
                spike_r(10,t)=1;
                i=tet_sum(jj,10);
                l0=normpdf(markAll_t14(i,2)*ones(1,length(spikeT0_t14)),markAll_t14(:,2)',smker).*normpdf(markAll_t14(i,3)*ones(1,length(spikeT0_t14)),markAll_t14(:,3)',smker).*normpdf(markAll_t14(i,4)*ones(1,length(spikeT0_t14)),markAll_t14(:,4)',smker).*normpdf(markAll_t14(i,5)*ones(1,length(spikeT0_t14)),markAll_t14(:,5)',smker);
                l1=Xnum_t14*l0'./occ(:,1)./dt;
                l2=l1.*dt.*exp(-Lint_t14.*dt);
                l2=l2./sum(l2);
                l(:,j)=l2;
            elseif tetVec(11)==1 %tet17
                spike_r(11,t)=1;
                i=tet_sum(jj,11);
                l0=normpdf(markAll_t17(i,2)*ones(1,length(spikeT0_t17)),markAll_t17(:,2)',smker).*normpdf(markAll_t17(i,3)*ones(1,length(spikeT0_t17)),markAll_t17(:,3)',smker).*normpdf(markAll_t17(i,4)*ones(1,length(spikeT0_t17)),markAll_t17(:,4)',smker).*normpdf(markAll_t17(i,5)*ones(1,length(spikeT0_t17)),markAll_t17(:,5)',smker);
                l1=Xnum_t17*l0'./occ(:,1)./dt;
                l2=l1.*dt.*exp(-Lint_t17.*dt);
                l2=l2./sum(l2);
                l(:,j)=l2;
            elseif tetVec(12)==1 %tet18
                spike_r(12,t)=1;
                i=tet_sum(jj,12);
                l0=normpdf(markAll_t18(i,2)*ones(1,length(spikeT0_t18)),markAll_t18(:,2)',smker).*normpdf(markAll_t18(i,3)*ones(1,length(spikeT0_t18)),markAll_t18(:,3)',smker).*normpdf(markAll_t18(i,4)*ones(1,length(spikeT0_t18)),markAll_t18(:,4)',smker).*normpdf(markAll_t18(i,5)*ones(1,length(spikeT0_t18)),markAll_t18(:,5)',smker);
                l1=Xnum_t18*l0'./occ(:,1)./dt;
                l2=l1.*dt.*exp(-Lint_t18.*dt);
                l2=l2./sum(l2);
                l(:,j)=l2;
            elseif tetVec(13)==1 %tet19
                spike_r(13,t)=1;
                i=tet_sum(jj,13);
                l0=normpdf(markAll_t19(i,2)*ones(1,length(spikeT0_t19)),markAll_t19(:,2)',smker).*normpdf(markAll_t19(i,3)*ones(1,length(spikeT0_t19)),markAll_t19(:,3)',smker).*normpdf(markAll_t19(i,4)*ones(1,length(spikeT0_t19)),markAll_t19(:,4)',smker).*normpdf(markAll_t19(i,5)*ones(1,length(spikeT0_t19)),markAll_t19(:,5)',smker);
                l1=Xnum_t19*l0'./occ(:,1)./dt;
                l2=l1.*dt.*exp(-Lint_t19.*dt);
                l2=l2./sum(l2);
                l(:,j)=l2;
            elseif tetVec(14)==1 %tet20
                spike_r(14,t)=1;
                i=tet_sum(jj,14);
                l0=normpdf(markAll_t20(i,2)*ones(1,length(spikeT0_t20)),markAll_t20(:,2)',smker).*normpdf(markAll_t20(i,3)*ones(1,length(spikeT0_t20)),markAll_t20(:,3)',smker).*normpdf(markAll_t20(i,4)*ones(1,length(spikeT0_t20)),markAll_t20(:,4)',smker).*normpdf(markAll_t20(i,5)*ones(1,length(spikeT0_t20)),markAll_t20(:,5)',smker);
                l1=Xnum_t20*l0'./occ(:,1)./dt;
                l2=l1.*dt.*exp(-Lint_t20.*dt);
                l2=l2./sum(l2);
                l(:,j)=l2;
            elseif tetVec(15)==1 %tet22
                spike_r(15,t)=1;
                i=tet_sum(jj,15);
                l0=normpdf(markAll_t22(i,2)*ones(1,length(spikeT0_t22)),markAll_t22(:,2)',smker).*normpdf(markAll_t22(i,3)*ones(1,length(spikeT0_t22)),markAll_t22(:,3)',smker).*normpdf(markAll_t22(i,4)*ones(1,length(spikeT0_t22)),markAll_t22(:,4)',smker).*normpdf(markAll_t22(i,5)*ones(1,length(spikeT0_t22)),markAll_t22(:,5)',smker);
                l1=Xnum_t22*l0'./occ(:,1)./dt;
                l2=l1.*dt.*exp(-Lint_t22.*dt);
                l2=l2./sum(l2);
                l(:,j)=l2;
            elseif tetVec(16)==1 %tet23
                spike_r(16,t)=1;
                i=tet_sum(jj,16);
                l0=normpdf(markAll_t23(i,2)*ones(1,length(spikeT0_t23)),markAll_t23(:,2)',smker).*normpdf(markAll_t23(i,3)*ones(1,length(spikeT0_t23)),markAll_t23(:,3)',smker).*normpdf(markAll_t23(i,4)*ones(1,length(spikeT0_t23)),markAll_t23(:,4)',smker).*normpdf(markAll_t23(i,5)*ones(1,length(spikeT0_t23)),markAll_t23(:,5)',smker);
                l1=Xnum_t23*l0'./occ(:,1)./dt;
                l2=l1.*dt.*exp(-Lint_t23.*dt);
                l2=l2./sum(l2);
                l(:,j)=l2;
            elseif tetVec(17)==1 %tet27
                spike_r(17,t)=1;
                i=tet_sum(jj,17);
                l0=normpdf(markAll_t27(i,2)*ones(1,length(spikeT0_t27)),markAll_t27(:,2)',smker).*normpdf(markAll_t27(i,3)*ones(1,length(spikeT0_t27)),markAll_t27(:,3)',smker).*normpdf(markAll_t27(i,4)*ones(1,length(spikeT0_t27)),markAll_t27(:,4)',smker).*normpdf(markAll_t27(i,5)*ones(1,length(spikeT0_t27)),markAll_t27(:,5)',smker);
                l1=Xnum_t27*l0'./occ(:,1)./dt;
                l2=l1.*dt.*exp(-Lint_t27.*dt);
                l2=l2./sum(l2);
                l(:,j)=l2;
            elseif tetVec(18)==1 %tet29
                spike_r(18,t)=1;
                i=tet_sum(jj,18);
                l0=normpdf(markAll_t29(i,2)*ones(1,length(spikeT0_t29)),markAll_t29(:,2)',smker).*normpdf(markAll_t29(i,3)*ones(1,length(spikeT0_t29)),markAll_t29(:,3)',smker).*normpdf(markAll_t29(i,4)*ones(1,length(spikeT0_t29)),markAll_t29(:,4)',smker).*normpdf(markAll_t29(i,5)*ones(1,length(spikeT0_t29)),markAll_t29(:,5)',smker);
                l1=Xnum_t29*l0'./occ(:,1)./dt;
                l2=l1.*dt.*exp(-Lint_t29.*dt);
                l2=l2./sum(l2);
                l(:,j)=l2;
            end
        end
        L=prod(l,2);
        L=L./sum(L);
    end
    totnorm=sum(onestep_I1.*L)+sum(onestep_I0.*L);
    postx_I1=onestep_I1.*L./totnorm;
    postx_I0=onestep_I0.*L./totnorm;
    pI1_vec(t)=sum(postx_I1);
    pI0_vec(t)=sum(postx_I0);
    postxM_r_I1(:,t)=postx_I1;
    postxM_r_I0(:,t)=postx_I0;
    clear onestep_I1 onestep_I2 a l L;
end
toc

figure;
plot(1:numSteps,pI1_vec,'r.',1:numSteps,pI0_vec,'b.');
title('Initial+State, clusterless','FontSize',14);
xlim([0 numSteps]);
ylim([0 1]);
xlabel('Time (ms)','FontSize',14);
set(gca,'FontSize',20);
legend('I=1','I=0','Location','SouthEast');



figure;
subplot(3,2,1); ind_seg1=1:21;
imagesc(1:numSteps,stateV(ind_seg1),postxM_r_I1(ind_seg1,:));
title('postx|I=1, clusterless','FontSize',18);
ylabel('Linearized position','FontSize',14);
set(gca,'FontSize',14,'YTick',[-3 -2.5 -2 -1.5 -1 -0.5 0 0.5 1 1.5 2 2.5 3],'YTickLabel',[-3 -2.5 -2 -1.5 -1 -0.5 0 0.5 1 1.5 2 2.5 3]);
colormap(flipud(hot(256)));
caxis([0 0.3]);
xlim([0 numSteps]);
subplot(3,2,3);
imagesc(1:numSteps,stateV(21:31),[(postxM_r_I1(21:30,:)+flipud(postxM_r_I1(32:41,:)))./2;postxM_r_I1(31,:)]);
ylabel('Linearized position','FontSize',14);
set(gca,'FontSize',14,'YTick',[-3 -2.5 -2 -1.5 -1 -0.5 0 0.5 1 1.5 2 2.5 3],'YTickLabel',[-3 -2.5 -2 -1.5 -1 -0.5 0 0.5 1 1.5 2 2.5 3]);
colormap(flipud(hot(256)));
caxis([0 0.3]);
xlim([0 numSteps]);
subplot(3,2,5); ind_seg3=41:61;
imagesc(1:numSteps,stateV(ind_seg3),postxM_r_I1(ind_seg3,:));
ylabel('Linearized position','FontSize',14); xlabel('Time (ms)','FontSize',14);
set(gca,'FontSize',14,'YTick',[-3 -2.5 -2 -1.5 -1 -0.5 0 0.5 1 1.5 2 2.5 3],'YTickLabel',[-3 -2.5 -2 -1.5 -1 -0.5 0 0.5 1 1.5 2 2.5 3]);
colormap(flipud(hot(256)));
caxis([0 0.3]);
xlim([0 numSteps]);

subplot(3,2,2); ind_seg1=1:21;
imagesc(1:numSteps,stateV(ind_seg1),postxM_r_I0(ind_seg1,:));
title('postx|I=0, clusterless','FontSize',18);
ylabel('Linearized position','FontSize',14);
set(gca,'FontSize',14,'YTick',[-3 -2.5 -2 -1.5 -1 -0.5 0 0.5 1 1.5 2 2.5 3],'YTickLabel',[-3 -2.5 -2 -1.5 -1 -0.5 0 0.5 1 1.5 2 2.5 3]);
colormap(flipud(hot(256)));
caxis([0 0.3]);
xlim([0 numSteps]);
subplot(3,2,4);
imagesc(1:numSteps,stateV(21:31),[(postxM_r_I0(21:30,:)+flipud(postxM_r_I0(32:41,:)))./2;postxM_r_I0(31,:)]);
ylabel('Linearized position','FontSize',14);
set(gca,'FontSize',14,'YTick',[-3 -2.5 -2 -1.5 -1 -0.5 0 0.5 1 1.5 2 2.5 3],'YTickLabel',[-3 -2.5 -2 -1.5 -1 -0.5 0 0.5 1 1.5 2 2.5 3]);
colormap(flipud(hot(256)));
caxis([0 0.3]);
xlim([0 numSteps]);
subplot(3,2,6); ind_seg3=41:61;
imagesc(1:numSteps,stateV(ind_seg3),postxM_r_I0(ind_seg3,:));
ylabel('Linearized position','FontSize',14); xlabel('Time (ms)','FontSize',14);
set(gca,'FontSize',14,'YTick',[-3 -2.5 -2 -1.5 -1 -0.5 0 0.5 1 1.5 2 2.5 3],'YTickLabel',[-3 -2.5 -2 -1.5 -1 -0.5 0 0.5 1 1.5 2 2.5 3]);
colormap(flipud(hot(256)));
caxis([0 0.3]);
xlim([0 numSteps]);
