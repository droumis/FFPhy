%load('bond04_vecLF_fold.mat');
%load('bond04_vecLF.mat');

load('/opt/data13/kkay/Bon/bonlinpos04.mat');

%%% basic data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
postimevec=linpos{4}{2}.statematrix.time;            
poslin=vecLF(:,2);  
    poslin2 = poslin';

%%% basic values:   xdel, xbins, dt %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xdel=0.1; 
xbins=min(poslin):xdel:max(poslin);        % -3 : 0.1 : 3
    numlinbins = length(xbins);
dt=postimevec(2)-postimevec(1);         

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% make State M %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stateM=[];
[~,binind]=histc(poslin,xbins);        % histogram of linearized positions
xbinstate=[binind(1:end-1) binind(2:end)];  % [ <current step lin pos>  <next step lin pos>  ]
%by column = departuring
for bb = 1:numlinbins
    % nextstates is the next states;
    nextstates = xbinstate( xbinstate(:,1)==bb , 2  ); 
    if ~isempty(nextstates)
        % histograms the nextstates and normalizes the distribution (why?)
        stateM(:,bb) = histc( nextstates , linspace(1,numlinbins,numlinbins))./size(nextstates,1 );
    elseif isempty(nextstates)
        stateM(:,bb) = zeros(1,numlinbins);
    end
end
K = inline('exp(-(x.^2+y.^2)/2/sig^2)');                        %gaussian
[dx,dy] = meshgrid([-1:1]);
sig = 0.5;
weight = K(sig,dx,dy)/sum(sum(K(sig,dx,dy)));                   %normalizing weights
stateM_gaus = conv2(stateM,weight,'same');                      %gaussian smoothed
stateM_gausnorm = stateM_gaus*diag(1./sum(stateM_gaus,1));      %normalized to confine probability to 1




%%%%%%  encode per tetrode
load('/mnt/vortex_data/mkarlsso/backup/bond/bond04/01-149/bond04-01_params.mat');
%ind_t1=find(filedata.params(:,1)/10000>=postimevec(1)&filedata.params(:,1)/10000<=postimevec(end)&filedata.params(:,8)>10&(filedata.params(:,2)>100|filedata.params(:,3)>100|filedata.params(:,4)>100|filedata.params(:,5)>100));
	% identify spikes that occur in transcribed epoch times
ind_t1= find( filedata.params(:,1)/10000  >=  postimevec(1) &  filedata.params(:,1)/10000 <= postimevec(end) );
    % spike timestamps
time_t1=filedata.params(ind_t1,1);
    % mark vector (4 channel amplitudes)
mark0_t1=[filedata.params(ind_t1,2) filedata.params(ind_t1,3) filedata.params(ind_t1,4) filedata.params(ind_t1,5)];
    % spike timestamps in sec
time2_t1=time_t1/10000;
    
spikeT0_t1=time2_t1;
[procInd0_t1,procInd1_t1]=histc(spikeT0_t1,postimevec);
procInd_t1=find(procInd0_t1);
spikeT_t1=postimevec(procInd_t1);
spike_t1=procInd0_t1';
[~,rawInd0_t1]=histc(spikeT0_t1,time2_t1);
markAll_t1(:,1)=procInd1_t1;markAll_t1(:,2:5)=mark0_t1(rawInd0_t1(rawInd0_t1~=0),:);
mdel=20; ms=min(min(markAll_t1(:,2:5))):mdel:max(max(markAll_t1(:,2:5))); 
sxker=2*xdel; smker = mdel; T=size(postimevec,1);
occ=normpdf(xbins'*ones(1,T),ones(length(xbins),1)*poslin2,sxker)*ones(T,length(ms));
%occ: columns are identical; occupancy based on position; denominator
Xnum_t1=normpdf(xbins'*ones(1,length(spikeT0_t1)),ones(length(xbins),1)*poslin2(procInd1_t1),sxker);
%Xnum: Gaussian kernel estimators for position
Lint_t1=sum(Xnum_t1,2)./occ(:,1)./dt; %integral
Lint_t1=Lint_t1./sum(Lint_t1);

load('/mnt/vortex_data/mkarlsso/backup/bond/bond04/02-151/bond04-02_params.mat');
%ind_t2=find(filedata.params(:,1)/10000>=postimevec(1)&filedata.params(:,1)/10000<=postimevec(end)&filedata.params(:,8)>10&(filedata.params(:,2)>100|filedata.params(:,3)>100|filedata.params(:,4)>100|filedata.params(:,5)>100));
ind_t2=find(filedata.params(:,1)/10000>=postimevec(1)&filedata.params(:,1)/10000<=postimevec(end));
time_t2=filedata.params(ind_t2,1);
mark0_t2=[filedata.params(ind_t2,2) filedata.params(ind_t2,3) filedata.params(ind_t2,4) filedata.params(ind_t2,5)];
time2_t2=time_t2/10000;
spikeT0_t2=time2_t2;
[procInd0_t2,procInd1_t2]=histc(spikeT0_t2,postimevec);
procInd_t2=find(procInd0_t2);
spikeT_t2=postimevec(procInd_t2);
spike_t2=procInd0_t2';
[~,rawInd0_t2]=histc(spikeT0_t2,time2_t2);
markAll_t2(:,1)=procInd1_t2;markAll_t2(:,2:5)=mark0_t2(rawInd0_t2(rawInd0_t2~=0),:);
mdel=20; ms=min(min(markAll_t2(:,2:5))):mdel:max(max(markAll_t2(:,2:5))); 
sxker=2*xdel; smker = mdel; T=size(postimevec,1);
occ=normpdf(xbins'*ones(1,T),ones(length(xbins),1)*poslin2,sxker)*ones(T,length(ms));
%occ: columns are identical; occupancy based on position; denominator
Xnum_t2=normpdf(xbins'*ones(1,length(spikeT0_t2)),ones(length(xbins),1)*poslin2(procInd1_t2),sxker);
%Xnum: Gaussian kernel estimators for position
Lint_t2=sum(Xnum_t2,2)./occ(:,1)./dt; %integral
Lint_t2=Lint_t2./sum(Lint_t2);

load('/mnt/vortex_data/mkarlsso/backup/bond/bond04/04-107/bond04-04_params.mat');
%ind_t4=find(filedata.params(:,1)/10000>=postimevec(1)&filedata.params(:,1)/10000<=postimevec(end)&filedata.params(:,8)>10&(filedata.params(:,2)>100|filedata.params(:,3)>100|filedata.params(:,4)>100|filedata.params(:,5)>100));
ind_t4=find(filedata.params(:,1)/10000>=postimevec(1)&filedata.params(:,1)/10000<=postimevec(end));
time_t4=filedata.params(ind_t4,1);
mark0_t4=[filedata.params(ind_t4,2) filedata.params(ind_t4,3) filedata.params(ind_t4,4) filedata.params(ind_t4,5)];
time2_t4=time_t4/10000;
spikeT0_t4=time2_t4;
[procInd0_t4,procInd1_t4]=histc(spikeT0_t4,postimevec);
procInd_t4=find(procInd0_t4);
spikeT_t4=postimevec(procInd_t4);
spike_t4=procInd0_t4';
[~,rawInd0_t4]=histc(spikeT0_t4,time2_t4);
markAll_t4(:,1)=procInd1_t4;markAll_t4(:,2:5)=mark0_t4(rawInd0_t4(rawInd0_t4~=0),:);
mdel=20; ms=min(min(markAll_t4(:,2:5))):mdel:max(max(markAll_t4(:,2:5))); 
sxker=2*xdel; smker = mdel; T=size(postimevec,1);
occ=normpdf(xbins'*ones(1,T),ones(length(xbins),1)*poslin2,sxker)*ones(T,length(ms));
%occ: columns are identical; occupancy based on position; denominator
Xnum_t4=normpdf(xbins'*ones(1,length(spikeT0_t4)),ones(length(xbins),1)*poslin2(procInd1_t4),sxker);
%Xnum: Gaussian kernel estimators for position
Lint_t4=sum(Xnum_t4,2)./occ(:,1)./dt; %integral
Lint_t4=Lint_t4./sum(Lint_t4);

load('/mnt/vortex_data/mkarlsso/backup/bond/bond04/05-105/bond04-05_params.mat');
%ind_t5=find(filedata.params(:,1)/10000>=postimevec(1)&filedata.params(:,1)/10000<=postimevec(end)&filedata.params(:,8)>10&(filedata.params(:,2)>100|filedata.params(:,3)>100|filedata.params(:,4)>100|filedata.params(:,5)>100));
ind_t5=find(filedata.params(:,1)/10000>=postimevec(1)&filedata.params(:,1)/10000<=postimevec(end));
time_t5=filedata.params(ind_t5,1);
mark0_t5=[filedata.params(ind_t5,2) filedata.params(ind_t5,3) filedata.params(ind_t5,4) filedata.params(ind_t5,5)];
time2_t5=time_t5/10000;
spikeT0_t5=time2_t5;
[procInd0_t5,procInd1_t5]=histc(spikeT0_t5,postimevec);
procInd_t5=find(procInd0_t5);
spikeT_t5=postimevec(procInd_t5);
spike_t5=procInd0_t5';
[~,rawInd0_t5]=histc(spikeT0_t5,time2_t5);
markAll_t5(:,1)=procInd1_t5;markAll_t5(:,2:5)=mark0_t5(rawInd0_t5(rawInd0_t5~=0),:);
mdel=20; ms=min(min(markAll_t5(:,2:5))):mdel:max(max(markAll_t5(:,2:5))); 
sxker=2*xdel; smker = mdel; T=size(postimevec,1);
occ=normpdf(xbins'*ones(1,T),ones(length(xbins),1)*poslin2,sxker)*ones(T,length(ms));
%occ: columns are identical; occupancy based on position; denominator
Xnum_t5=normpdf(xbins'*ones(1,length(spikeT0_t5)),ones(length(xbins),1)*poslin2(procInd1_t5),sxker);
%Xnum: Gaussian kernel estimators for position
Lint_t5=sum(Xnum_t5,2)./occ(:,1)./dt; %integral
Lint_t5=Lint_t5./sum(Lint_t5);

load('/mnt/vortex_data/mkarlsso/backup/bond/bond04/07-162/bond04-07_params.mat');
%ind_t7=find(filedata.params(:,1)/10000>=postimevec(1)&filedata.params(:,1)/10000<=postimevec(end)&filedata.params(:,8)>10&(filedata.params(:,2)>100|filedata.params(:,3)>100|filedata.params(:,4)>100|filedata.params(:,5)>100));
ind_t7=find(filedata.params(:,1)/10000>=postimevec(1)&filedata.params(:,1)/10000<=postimevec(end));
time_t7=filedata.params(ind_t7,1);
mark0_t7=[filedata.params(ind_t7,2) filedata.params(ind_t7,3) filedata.params(ind_t7,4) filedata.params(ind_t7,5)];
time2_t7=time_t7/10000;
spikeT0_t7=time2_t7;
[procInd0_t7,procInd1_t7]=histc(spikeT0_t7,postimevec);
procInd_t7=find(procInd0_t7);
spikeT_t7=postimevec(procInd_t7);
spike_t7=procInd0_t7';
[~,rawInd0_t7]=histc(spikeT0_t7,time2_t7);
markAll_t7(:,1)=procInd1_t7;markAll_t7(:,2:5)=mark0_t7(rawInd0_t7(rawInd0_t7~=0),:);
mdel=20; ms=min(min(markAll_t7(:,2:5))):mdel:max(max(markAll_t7(:,2:5))); 
sxker=2*xdel; smker = mdel; T=size(postimevec,1);
occ=normpdf(xbins'*ones(1,T),ones(length(xbins),1)*poslin2,sxker)*ones(T,length(ms));
%occ: columns are identical; occupancy based on position; denominator
Xnum_t7=normpdf(xbins'*ones(1,length(spikeT0_t7)),ones(length(xbins),1)*poslin2(procInd1_t7),sxker);
%Xnum: Gaussian kernel estimators for position
Lint_t7=sum(Xnum_t7,2)./occ(:,1)./dt; %integral
Lint_t7=Lint_t7./sum(Lint_t7);

load('/mnt/vortex_data/mkarlsso/backup/bond/bond04/10-146/bond04-10_params.mat');
%ind_t10=find(filedata.params(:,1)/10000>=postimevec(1)&filedata.params(:,1)/10000<=postimevec(end)&filedata.params(:,8)>10&(filedata.params(:,2)>100|filedata.params(:,3)>100|filedata.params(:,4)>100|filedata.params(:,5)>100));
ind_t10=find(filedata.params(:,1)/10000>=postimevec(1)&filedata.params(:,1)/10000<=postimevec(end));
time_t10=filedata.params(ind_t10,1);
mark0_t10=[filedata.params(ind_t10,2) filedata.params(ind_t10,3) filedata.params(ind_t10,4) filedata.params(ind_t10,5)];
time2_t10=time_t10/10000;
spikeT0_t10=time2_t10;
[procInd0_t10,procInd1_t10]=histc(spikeT0_t10,postimevec);
procInd_t10=find(procInd0_t10);
spikeT_t10=postimevec(procInd_t10);
spike_t10=procInd0_t10';
[~,rawInd0_t10]=histc(spikeT0_t10,time2_t10);
markAll_t10(:,1)=procInd1_t10;markAll_t10(:,2:5)=mark0_t10(rawInd0_t10(rawInd0_t10~=0),:);
mdel=20; ms=min(min(markAll_t10(:,2:5))):mdel:max(max(markAll_t10(:,2:5))); 
sxker=2*xdel; smker = mdel; T=size(postimevec,1);
occ=normpdf(xbins'*ones(1,T),ones(length(xbins),1)*poslin2,sxker)*ones(T,length(ms));
%occ: columns are identical; occupancy based on position; denominator
Xnum_t10=normpdf(xbins'*ones(1,length(spikeT0_t10)),ones(length(xbins),1)*poslin2(procInd1_t10),sxker);
%Xnum: Gaussian kernel estimators for position
Lint_t10=sum(Xnum_t10,2)./occ(:,1)./dt; %integral
Lint_t10=Lint_t10./sum(Lint_t10);

load('/mnt/vortex_data/mkarlsso/backup/bond/bond04/11-096/bond04-11_params.mat');
%ind_t11=find(filedata.params(:,1)/10000>=postimevec(1)&filedata.params(:,1)/10000<=postimevec(end)&filedata.params(:,8)>10&(filedata.params(:,2)>100|filedata.params(:,3)>100|filedata.params(:,4)>100|filedata.params(:,5)>100));
ind_t11=find(filedata.params(:,1)/10000>=postimevec(1)&filedata.params(:,1)/10000<=postimevec(end));
time_t11=filedata.params(ind_t11,1);
mark0_t11=[filedata.params(ind_t11,2) filedata.params(ind_t11,3) filedata.params(ind_t11,4) filedata.params(ind_t11,5)];
time2_t11=time_t11/10000;
spikeT0_t11=time2_t11;
[procInd0_t11,procInd1_t11]=histc(spikeT0_t11,postimevec);
procInd_t11=find(procInd0_t11);
spikeT_t11=postimevec(procInd_t11);
spike_t11=procInd0_t11';
[~,rawInd0_t11]=histc(spikeT0_t11,time2_t11);
markAll_t11(:,1)=procInd1_t11;markAll_t11(:,2:5)=mark0_t11(rawInd0_t11(rawInd0_t11~=0),:);
mdel=20; ms=min(min(markAll_t11(:,2:5))):mdel:max(max(markAll_t11(:,2:5)));
sxker=2*xdel; smker = mdel; T=size(postimevec,1);
occ=normpdf(xbins'*ones(1,T),ones(length(xbins),1)*poslin2,sxker)*ones(T,length(ms));
%occ: columns are identical; occupancy based on position; denominator
Xnum_t11=normpdf(xbins'*ones(1,length(spikeT0_t11)),ones(length(xbins),1)*poslin2(procInd1_t11),sxker);
%Xnum: Gaussian kernel estimators for position
Lint_t11=sum(Xnum_t11,2)./occ(:,1)./dt; %integral
Lint_t11=Lint_t11./sum(Lint_t11);

load('/mnt/vortex_data/mkarlsso/backup/bond/bond04/12-136/bond04-12_params.mat');
%ind_t12=find(filedata.params(:,1)/10000>=postimevec(1)&filedata.params(:,1)/10000<=postimevec(end)&filedata.params(:,8)>10&(filedata.params(:,2)>100|filedata.params(:,3)>100|filedata.params(:,4)>100|filedata.params(:,5)>100));
ind_t12=find(filedata.params(:,1)/10000>=postimevec(1)&filedata.params(:,1)/10000<=postimevec(end));
time_t12=filedata.params(ind_t12,1);
mark0_t12=[filedata.params(ind_t12,2) filedata.params(ind_t12,3) filedata.params(ind_t12,4) filedata.params(ind_t12,5)];
time2_t12=time_t12/10000;
spikeT0_t12=time2_t12;
[procInd0_t12,procInd1_t12]=histc(spikeT0_t12,postimevec);
procInd_t12=find(procInd0_t12);
spikeT_t12=postimevec(procInd_t12);
spike_t12=procInd0_t12';
[~,rawInd0_t12]=histc(spikeT0_t12,time2_t12);
markAll_t12(:,1)=procInd1_t12;markAll_t12(:,2:5)=mark0_t12(rawInd0_t12(rawInd0_t12~=0),:);
mdel=20; ms=min(min(markAll_t12(:,2:5))):mdel:max(max(markAll_t12(:,2:5)));
sxker=2*xdel; smker = mdel; T=size(postimevec,1);
occ=normpdf(xbins'*ones(1,T),ones(length(xbins),1)*poslin2,sxker)*ones(T,length(ms));
%occ: columns are identical; occupancy based on position; denominator
Xnum_t12=normpdf(xbins'*ones(1,length(spikeT0_t12)),ones(length(xbins),1)*poslin2(procInd1_t12),sxker);
%Xnum: Gaussian kernel estimators for position
Lint_t12=sum(Xnum_t12,2)./occ(:,1)./dt; %integral
Lint_t12=Lint_t12./sum(Lint_t12);

load('/mnt/vortex_data/mkarlsso/backup/bond/bond04/13-094/bond04-13_params.mat');
%ind_t13=find(filedata.params(:,1)/10000>=postimevec(1)&filedata.params(:,1)/10000<=postimevec(end)&filedata.params(:,8)>10&(filedata.params(:,2)>100|filedata.params(:,3)>100|filedata.params(:,4)>100|filedata.params(:,5)>100));
ind_t13=find(filedata.params(:,1)/10000>=postimevec(1)&filedata.params(:,1)/10000<=postimevec(end));
time_t13=filedata.params(ind_t13,1);
mark0_t13=[filedata.params(ind_t13,2) filedata.params(ind_t13,3) filedata.params(ind_t13,4) filedata.params(ind_t13,5)];
time2_t13=time_t13/10000;
spikeT0_t13=time2_t13;
[procInd0_t13,procInd1_t13]=histc(spikeT0_t13,postimevec);
procInd_t13=find(procInd0_t13);
spikeT_t13=postimevec(procInd_t13);
spike_t13=procInd0_t13';
[~,rawInd0_t13]=histc(spikeT0_t13,time2_t13);
markAll_t13(:,1)=procInd1_t13;markAll_t13(:,2:5)=mark0_t13(rawInd0_t13(rawInd0_t13~=0),:);
mdel=20; ms=min(min(markAll_t13(:,2:5))):mdel:max(max(markAll_t13(:,2:5)));
sxker=2*xdel; smker = mdel; T=size(postimevec,1);
occ=normpdf(xbins'*ones(1,T),ones(length(xbins),1)*poslin2,sxker)*ones(T,length(ms));
%occ: columns are identical; occupancy based on position; denominator
Xnum_t13=normpdf(xbins'*ones(1,length(spikeT0_t13)),ones(length(xbins),1)*poslin2(procInd1_t13),sxker);
%Xnum: Gaussian kernel estimators for position
Lint_t13=sum(Xnum_t13,2)./occ(:,1)./dt; %integral
Lint_t13=Lint_t13./sum(Lint_t13);

load('/mnt/vortex_data/mkarlsso/backup/bond/bond04/14-096/bond04-14_params.mat');
%ind_t14=find(filedata.params(:,1)/10000>=postimevec(1)&filedata.params(:,1)/10000<=postimevec(end)&filedata.params(:,8)>10&(filedata.params(:,2)>100|filedata.params(:,3)>100|filedata.params(:,4)>100|filedata.params(:,5)>100));
ind_t14=find(filedata.params(:,1)/10000>=postimevec(1)&filedata.params(:,1)/10000<=postimevec(end));
time_t14=filedata.params(ind_t14,1);
mark0_t14=[filedata.params(ind_t14,2) filedata.params(ind_t14,3) filedata.params(ind_t14,4) filedata.params(ind_t14,5)];
time2_t14=time_t14/10000;
spikeT0_t14=time2_t14;
[procInd0_t14,procInd1_t14]=histc(spikeT0_t14,postimevec);
procInd_t14=find(procInd0_t14);
spikeT_t14=postimevec(procInd_t14);
spike_t14=procInd0_t14';
[~,rawInd0_t14]=histc(spikeT0_t14,time2_t14);
markAll_t14(:,1)=procInd1_t14;markAll_t14(:,2:5)=mark0_t14(rawInd0_t14(rawInd0_t14~=0),:);
mdel=20; ms=min(min(markAll_t14(:,2:5))):mdel:max(max(markAll_t14(:,2:5)));
sxker=2*xdel; smker = mdel; T=size(postimevec,1);
occ=normpdf(xbins'*ones(1,T),ones(length(xbins),1)*poslin2,sxker)*ones(T,length(ms));
%occ: columns are identical; occupancy based on position; denominator
Xnum_t14=normpdf(xbins'*ones(1,length(spikeT0_t14)),ones(length(xbins),1)*poslin2(procInd1_t14),sxker);
%Xnum: Gaussian kernel estimators for position
Lint_t14=sum(Xnum_t14,2)./occ(:,1)./dt; %integral
Lint_t14=Lint_t14./sum(Lint_t14);

load('/mnt/vortex_data/mkarlsso/backup/bond/bond04/17-108/bond04-17_params.mat');
%ind_t17=find(filedata.params(:,1)/10000>=postimevec(1)&filedata.params(:,1)/10000<=postimevec(end)&filedata.params(:,8)>10&(filedata.params(:,2)>100|filedata.params(:,3)>100|filedata.params(:,4)>100|filedata.params(:,5)>100));
ind_t17=find(filedata.params(:,1)/10000>=postimevec(1)&filedata.params(:,1)/10000<=postimevec(end));
time_t17=filedata.params(ind_t17,1);
mark0_t17=[filedata.params(ind_t17,2) filedata.params(ind_t17,3) filedata.params(ind_t17,4) filedata.params(ind_t17,5)];
time2_t17=time_t17/10000;
spikeT0_t17=time2_t17;
[procInd0_t17,procInd1_t17]=histc(spikeT0_t17,postimevec);
procInd_t17=find(procInd0_t17);
spikeT_t17=postimevec(procInd_t17);
spike_t17=procInd0_t17';
[~,rawInd0_t17]=histc(spikeT0_t17,time2_t17);
markAll_t17(:,1)=procInd1_t17;markAll_t17(:,2:5)=mark0_t17(rawInd0_t17(rawInd0_t17~=0),:);
mdel=20; ms=min(min(markAll_t17(:,2:5))):mdel:max(max(markAll_t17(:,2:5)));
sxker=2*xdel; smker = mdel; T=size(postimevec,1);
occ=normpdf(xbins'*ones(1,T),ones(length(xbins),1)*poslin2,sxker)*ones(T,length(ms));
%occ: columns are identical; occupancy based on position; denominator
Xnum_t17=normpdf(xbins'*ones(1,length(spikeT0_t17)),ones(length(xbins),1)*poslin2(procInd1_t17),sxker);
%Xnum: Gaussian kernel estimators for position
Lint_t17=sum(Xnum_t17,2)./occ(:,1)./dt; %integral
Lint_t17=Lint_t17./sum(Lint_t17);

load('/mnt/vortex_data/mkarlsso/backup/bond/bond04/18-138/bond04-18_params.mat');
%ind_t18=find(filedata.params(:,1)/10000>=postimevec(1)&filedata.params(:,1)/10000<=postimevec(end)&filedata.params(:,8)>10&(filedata.params(:,2)>100|filedata.params(:,3)>100|filedata.params(:,4)>100|filedata.params(:,5)>100));
ind_t18=find(filedata.params(:,1)/10000>=postimevec(1)&filedata.params(:,1)/10000<=postimevec(end));
time_t18=filedata.params(ind_t18,1);
mark0_t18=[filedata.params(ind_t18,2) filedata.params(ind_t18,3) filedata.params(ind_t18,4) filedata.params(ind_t18,5)];
time2_t18=time_t18/10000;
spikeT0_t18=time2_t18;
[procInd0_t18,procInd1_t18]=histc(spikeT0_t18,postimevec);
procInd_t18=find(procInd0_t18);
spikeT_t18=postimevec(procInd_t18);
spike_t18=procInd0_t18';
[~,rawInd0_t18]=histc(spikeT0_t18,time2_t18);
markAll_t18(:,1)=procInd1_t18;markAll_t18(:,2:5)=mark0_t18(rawInd0_t18(rawInd0_t18~=0),:);
mdel=20; ms=min(min(markAll_t18(:,2:5))):mdel:max(max(markAll_t18(:,2:5)));
sxker=2*xdel; smker = mdel; T=size(postimevec,1);
occ=normpdf(xbins'*ones(1,T),ones(length(xbins),1)*poslin2,sxker)*ones(T,length(ms));
%occ: columns are identical; occupancy based on position; denominator
Xnum_t18=normpdf(xbins'*ones(1,length(spikeT0_t18)),ones(length(xbins),1)*poslin2(procInd1_t18),sxker);
%Xnum: Gaussian kernel estimators for position
Lint_t18=sum(Xnum_t18,2)./occ(:,1)./dt; %integral
Lint_t18=Lint_t18./sum(Lint_t18);

load('/mnt/vortex_data/mkarlsso/backup/bond/bond04/19-130/bond04-19_params.mat');
%ind_t19=find(filedata.params(:,1)/10000>=postimevec(1)&filedata.params(:,1)/10000<=postimevec(end)&filedata.params(:,8)>10&(filedata.params(:,2)>100|filedata.params(:,3)>100|filedata.params(:,4)>100|filedata.params(:,5)>100));
ind_t19=find(filedata.params(:,1)/10000>=postimevec(1)&filedata.params(:,1)/10000<=postimevec(end));
time_t19=filedata.params(ind_t19,1);
mark0_t19=[filedata.params(ind_t19,2) filedata.params(ind_t19,3) filedata.params(ind_t19,4) filedata.params(ind_t19,5)];
time2_t19=time_t19/10000;
spikeT0_t19=time2_t19;
[procInd0_t19,procInd1_t19]=histc(spikeT0_t19,postimevec);
procInd_t19=find(procInd0_t19);
spikeT_t19=postimevec(procInd_t19);
spike_t19=procInd0_t19';
[~,rawInd0_t19]=histc(spikeT0_t19,time2_t19);
markAll_t19(:,1)=procInd1_t19;markAll_t19(:,2:5)=mark0_t19(rawInd0_t19(rawInd0_t19~=0),:);
mdel=20; ms=min(min(markAll_t19(:,2:5))):mdel:max(max(markAll_t19(:,2:5)));
sxker=2*xdel; smker = mdel; T=size(postimevec,1);
occ=normpdf(xbins'*ones(1,T),ones(length(xbins),1)*poslin2,sxker)*ones(T,length(ms));
%occ: columns are identical; occupancy based on position; denominator
Xnum_t19=normpdf(xbins'*ones(1,length(spikeT0_t19)),ones(length(xbins),1)*poslin2(procInd1_t19),sxker);
%Xnum: Gaussian kernel estimators for position
Lint_t19=sum(Xnum_t19,2)./occ(:,1)./dt; %integral
Lint_t19=Lint_t19./sum(Lint_t19);

load('/mnt/vortex_data/mkarlsso/backup/bond/bond04/20-132/bond04-20_params.mat');
%ind_t20=find(filedata.params(:,1)/10000>=postimevec(1)&filedata.params(:,1)/10000<=postimevec(end)&filedata.params(:,8)>10&(filedata.params(:,2)>100|filedata.params(:,3)>100|filedata.params(:,4)>100|filedata.params(:,5)>100));
ind_t20=find(filedata.params(:,1)/10000>=postimevec(1)&filedata.params(:,1)/10000<=postimevec(end));
time_t20=filedata.params(ind_t20,1);
mark0_t20=[filedata.params(ind_t20,2) filedata.params(ind_t20,3) filedata.params(ind_t20,4) filedata.params(ind_t20,5)];
time2_t20=time_t20/10000;
spikeT0_t20=time2_t20;
[procInd0_t20,procInd1_t20]=histc(spikeT0_t20,postimevec);
procInd_t20=find(procInd0_t20);
spikeT_t20=postimevec(procInd_t20);
spike_t20=procInd0_t20';
[~,rawInd0_t20]=histc(spikeT0_t20,time2_t20);
markAll_t20(:,1)=procInd1_t20;markAll_t20(:,2:5)=mark0_t20(rawInd0_t20(rawInd0_t20~=0),:);
mdel=20; ms=min(min(markAll_t20(:,2:5))):mdel:max(max(markAll_t20(:,2:5)));
sxker=2*xdel; smker = mdel; T=size(postimevec,1);
occ=normpdf(xbins'*ones(1,T),ones(length(xbins),1)*poslin2,sxker)*ones(T,length(ms));
%occ: columns are identical; occupancy based on position; denominator
Xnum_t20=normpdf(xbins'*ones(1,length(spikeT0_t20)),ones(length(xbins),1)*poslin2(procInd1_t20),sxker);
%Xnum: Gaussian kernel estimators for position
Lint_t20=sum(Xnum_t20,2)./occ(:,1)./dt; %integral
Lint_t20=Lint_t20./sum(Lint_t20);

load('/mnt/vortex_data/mkarlsso/backup/bond/bond04/22-171/bond04-22_params.mat');
%ind_t22=find(filedata.params(:,1)/10000>=postimevec(1)&filedata.params(:,1)/10000<=postimevec(end)&filedata.params(:,8)>10&(filedata.params(:,2)>100|filedata.params(:,3)>100|filedata.params(:,4)>100|filedata.params(:,5)>100));
ind_t22=find(filedata.params(:,1)/10000>=postimevec(1)&filedata.params(:,1)/10000<=postimevec(end));
time_t22=filedata.params(ind_t22,1);
mark0_t22=[filedata.params(ind_t22,2) filedata.params(ind_t22,3) filedata.params(ind_t22,4) filedata.params(ind_t22,5)];
time2_t22=time_t22/10000;
spikeT0_t22=time2_t22;
[procInd0_t22,procInd1_t22]=histc(spikeT0_t22,postimevec);
procInd_t22=find(procInd0_t22);
spikeT_t22=postimevec(procInd_t22);
spike_t22=procInd0_t22';
[~,rawInd0_t22]=histc(spikeT0_t22,time2_t22);
markAll_t22(:,1)=procInd1_t22;markAll_t22(:,2:5)=mark0_t22(rawInd0_t22(rawInd0_t22~=0),:);
mdel=20; ms=min(min(markAll_t22(:,2:5))):mdel:max(max(markAll_t22(:,2:5)));
sxker=2*xdel; smker = mdel; T=size(postimevec,1);
occ=normpdf(xbins'*ones(1,T),ones(length(xbins),1)*poslin2,sxker)*ones(T,length(ms));
%occ: columns are identical; occupancy based on position; denominator
Xnum_t22=normpdf(xbins'*ones(1,length(spikeT0_t22)),ones(length(xbins),1)*poslin2(procInd1_t22),sxker);
%Xnum: Gaussian kernel estimators for position
Lint_t22=sum(Xnum_t22,2)./occ(:,1)./dt; %integral
Lint_t22=Lint_t22./sum(Lint_t22);

load('/mnt/vortex_data/mkarlsso/backup/bond/bond04/23-152/bond04-23_params.mat');
%ind_t23=find(filedata.params(:,1)/10000>=postimevec(1)&filedata.params(:,1)/10000<=postimevec(end)&filedata.params(:,8)>10&(filedata.params(:,2)>100|filedata.params(:,3)>100|filedata.params(:,4)>100|filedata.params(:,5)>100));
ind_t23=find(filedata.params(:,1)/10000>=postimevec(1)&filedata.params(:,1)/10000<=postimevec(end));
time_t23=filedata.params(ind_t23,1);
mark0_t23=[filedata.params(ind_t23,2) filedata.params(ind_t23,3) filedata.params(ind_t23,4) filedata.params(ind_t23,5)];
time2_t23=time_t23/10000;
spikeT0_t23=time2_t23;
[procInd0_t23,procInd1_t23]=histc(spikeT0_t23,postimevec);
procInd_t23=find(procInd0_t23);
spikeT_t23=postimevec(procInd_t23);
spike_t23=procInd0_t23';
[~,rawInd0_t23]=histc(spikeT0_t23,time2_t23);
markAll_t23(:,1)=procInd1_t23;markAll_t23(:,2:5)=mark0_t23(rawInd0_t23(rawInd0_t23~=0),:);
mdel=20; ms=min(min(markAll_t23(:,2:5))):mdel:max(max(markAll_t23(:,2:5)));
sxker=2*xdel; smker = mdel; T=size(postimevec,1);
occ=normpdf(xbins'*ones(1,T),ones(length(xbins),1)*poslin2,sxker)*ones(T,length(ms));
%occ: columns are identical; occupancy based on position; denominator
Xnum_t23=normpdf(xbins'*ones(1,length(spikeT0_t23)),ones(length(xbins),1)*poslin2(procInd1_t23),sxker);
%Xnum: Gaussian kernel estimators for position
Lint_t23=sum(Xnum_t23,2)./occ(:,1)./dt; %integral
Lint_t23=Lint_t23./sum(Lint_t23);

load('/mnt/vortex_data/mkarlsso/backup/bond/bond04/27-158/bond04-27_params.mat');
%ind_t27=find(filedata.params(:,1)/10000>=postimevec(1)&filedata.params(:,1)/10000<=postimevec(end)&filedata.params(:,8)>10&(filedata.params(:,2)>100|filedata.params(:,3)>100|filedata.params(:,4)>100|filedata.params(:,5)>100));
ind_t27=find(filedata.params(:,1)/10000>=postimevec(1)&filedata.params(:,1)/10000<=postimevec(end));
time_t27=filedata.params(ind_t27,1);
mark0_t27=[filedata.params(ind_t27,2) filedata.params(ind_t27,3) filedata.params(ind_t27,4) filedata.params(ind_t27,5)];
time2_t27=time_t27/10000;
spikeT0_t27=time2_t27;
[procInd0_t27,procInd1_t27]=histc(spikeT0_t27,postimevec);
procInd_t27=find(procInd0_t27);
spikeT_t27=postimevec(procInd_t27);
spike_t27=procInd0_t27';
[~,rawInd0_t27]=histc(spikeT0_t27,time2_t27);
markAll_t27(:,1)=procInd1_t27;markAll_t27(:,2:5)=mark0_t27(rawInd0_t27(rawInd0_t27~=0),:);
mdel=20; ms=min(min(markAll_t27(:,2:5))):mdel:max(max(markAll_t27(:,2:5)));
sxker=2*xdel; smker = mdel; T=size(postimevec,1);
occ=normpdf(xbins'*ones(1,T),ones(length(xbins),1)*poslin2,sxker)*ones(T,length(ms));
%occ: columns are identical; occupancy based on position; denominator
Xnum_t27=normpdf(xbins'*ones(1,length(spikeT0_t27)),ones(length(xbins),1)*poslin2(procInd1_t27),sxker);
%Xnum: Gaussian kernel estimators for position
Lint_t27=sum(Xnum_t27,2)./occ(:,1)./dt; %integral
Lint_t27=Lint_t27./sum(Lint_t27);

load('/mnt/vortex_data/mkarlsso/backup/bond/bond04/29-118/bond04-29_params.mat');
%ind_t29=find(filedata.params(:,1)/10000>=postimevec(1)&filedata.params(:,1)/10000<=postimevec(end)&filedata.params(:,8)>10&(filedata.params(:,2)>100|filedata.params(:,3)>100|filedata.params(:,4)>100|filedata.params(:,5)>100));
ind_t29=find(filedata.params(:,1)/10000>=postimevec(1)&filedata.params(:,1)/10000<=postimevec(end));
time_t29=filedata.params(ind_t29,1);
mark0_t29=[filedata.params(ind_t29,2) filedata.params(ind_t29,3) filedata.params(ind_t29,4) filedata.params(ind_t29,5)];
time2_t29=time_t29/10000;
spikeT0_t29=time2_t29;
[procInd0_t29,procInd1_t29]=histc(spikeT0_t29,postimevec);
procInd_t29=find(procInd0_t29);
spikeT_t29=postimevec(procInd_t29);
spike_t29=procInd0_t29';
[~,rawInd0_t29]=histc(spikeT0_t29,time2_t29);
markAll_t29(:,1)=procInd1_t29;markAll_t29(:,2:5)=mark0_t29(rawInd0_t29(rawInd0_t29~=0),:);
mdel=20; ms=min(min(markAll_t29(:,2:5))):mdel:max(max(markAll_t29(:,2:5)));
sxker=2*xdel; smker = mdel; T=size(postimevec,1);
occ=normpdf(xbins'*ones(1,T),ones(length(xbins),1)*poslin2,sxker)*ones(T,length(ms));
%occ: columns are identical; occupancy based on position; denominator
Xnum_t29=normpdf(xbins'*ones(1,length(spikeT0_t29)),ones(length(xbins),1)*poslin2(procInd1_t29),sxker);
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
sxker=2*xdel; smker = mdel; T=size(postimevec,1);
occ=normpdf(xbins'*ones(1,T),ones(length(xbins),1)*poslin2,sxker)*ones(T,length(ms));
%occ: columns are identical; occupancy based on position; denominator
Xnum=normpdf(xbins'*ones(1,length(time0)),ones(length(xbins),1)*poslin2(procInd1),sxker);
%Xnum: Gaussian kernel estimators for position
Lint=sum(Xnum,2)./occ(:,1)./dt; %integral
Lint=Lint./sum(Lint);
%Lint: conditional intensity function for the unmarked case



%%%% Ripples %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ripoption = 2;   % == 1 for old ripples, == 2 for ripplescons


%%
if ripoption == 1
    
    load('/opt/data13/kkay/Bon/bonripples04.mat');
    %ind_t=[1 2 4 5 7 10 11 12 13 14 17 18 19 20 22 23 27 29];
    ex=4;ep=2;
    
    ti=round(postimevec(1)*1000):1:round(postimevec(end)*1000);   % 1-ms time vector
    %position time; from 2461011 to 3404979; length(ti)~~max(spike_tim)
    
    % Make ind_ripple 
    ind_ripple=zeros(size(ti,2),size(ind_t,2));  % <time bin> x <tetnum>
    % iterate tetrodes
    ripple_all = [];
    for ii=1:size(ind_t,2)
        k=ind_t(ii);  % tet
        ripple_s=round(ripples{ex}{ep}{k}.starttime*1000);
        ripple_e=round(ripples{ex}{ep}{k}.endtime*1000);
        for l=1:size(ripple_s,1)
            ripple=[ripple_s(l) ripple_e(l)];
            ripple_all=[ripple_all;ripple];
        end
        clear ripple_s ripple_e ripple;
        for s=1:size(ripple_all,1)
            ind_ripple( find(  ti>=ripple_all(s,1)   &   ti<=ripple_all(s,2)  ),ii )=1;
        end
        ripple_all=[];
    end
    
    ripple_acr=zeros(size(ind_ripple,1),1);
    %ripple_acr=sum(ind_ripple,2);
    ripple_acr(   find(  sum(ind_ripple,2)  )   )=1;
    %%%indicator for ripples pooled across tetrodes
    
    %ripple_acr(find(sum(ind_ripple,2)>=6))=1;
    %%%indicator for ripples pooled across tetrodes where more than 6/9
    %%%tetrodes have spiking   ************** erroneous
    
    clear startT endT
    % iterate through 1-ms bin times
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


elseif ripoption == 2

    day=4;ep=2;    
    %ind_t=[1 2 4 5 7 10 11 12 13 14 17 18 19 20 22 23 27 29];
    animalname = 'Bond';
    animalinfo = animaldef(animalname);
    
    cons_name1 = 'ripplescons';
    tetfilter1 = 1;     % 1 for validripples CA1, 2 for CA3-DG
    consensus_numtets_ripc = 3;
    minthresh_ripc = 2;
    exclusion_dur_ripc = 0;
    minvelocity_ripc = 0;
    maxvelocity_ripc = 4;    
    output1 = kk_getconstimes(animalinfo{2},animalinfo{3},[day ep],cons_name1,tetfilter1,...
        'consensus_numtets',consensus_numtets_ripc,'minthresh',minthresh_ripc,...
        'exclusion_dur',exclusion_dur_ripc,'minvelocity',minvelocity_ripc,'maxvelocity',maxvelocity_ripc);
    consvec = output1{day}{ep}.cons;
    consvectimes = output1{day}{ep}.time;
    riplist = vec2list(consvec,consvectimes);
    
    ti=round(postimevec(1)*1000):1:round(postimevec(end)*1000);   % 1-ms time vector (from position time)
    
    ind_ripple=list2vec(riplist,ti/1000);  % <time bin> x <consvec>
    
    
    ripple_acr=zeros(size(ind_ripple,1),1);
    %ripple_acr=sum(ind_ripple,2);
    ripple_acr(   find(  sum(ind_ripple,2)  )   )=1;
    %%%indicator for ripples pooled across tetrodes
    
    %ripple_acr(find(sum(ind_ripple,2)>=6))=1;
    %%%indicator for ripples pooled across tetrodes where more than 6/9
    %%%tetrodes have spiking   ************** erroneous
    
    clear startT endT
    % iterate through 1-ms bin times
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
    
end


%%
ripple_dur=ripple_seg(:,2)-ripple_seg(:,1);
%1203, 1237

xi=round(time/10);  %spike times
% tot_ripple=zeros(length(ripple_dur),1);
% for i=1:length(ripple_dur)
%     numSteps=ripple_dur(i);
%     spike_tim=ripple_seg(i,1):ripple_seg(i,2); 
%     a1=0;
%     for t=1:numSteps
%         tt=spike_tim(t);
%         a=find(xi==ti(tt));
%         a1=a1+length(a);
%     end
%     tot_ripple(i)=a1;
% end


numrips = size(ripple_seg,1);

for ripnum = 1:numrips

%%
%spike_tim=ripple_seg(548,1):ripple_seg(548,2); %from 1 to 90000~
spike_tim=ripple_seg(ripnum,1):ripple_seg(ripnum,2); %from 1 to 90000~

%%
clear onstep postx postxM_r;
n=size(stateM_gausnorm,1); %# of position; length(xbins)
postx=ones(n,1)./n; %random prior
numSteps=length(spike_tim);
postxM_r=zeros(n,numSteps);

spike_r=zeros(length(ind_t),numSteps);
xi=round(time/10);  %spike times

tic;
for t=1:numSteps
    tt=spike_tim(t);
    onestep=stateM_gausnorm*postx;
    L=ones(size(postx));
    a=find(xi==ti(tt));
    
    if isempty(a)==1 %if no spike occurs at time t
        L=exp(-Lint.*dt);
    elseif isempty(a)==0 %if spikes
        l=zeros(n,length(a));
        for j=1:length(a)
            jj=a(j);
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
    postx=onestep.*L./sum(onestep.*L);
    postxM_r(:,t)=postx;
    clear onestep a l L;
end
toc



spike_r_raw=spike_r;



%%


figure;
subplot(3,2,[1 3 5]);
imagesc(spike_r);colormap(flipud(gray));freezeColors
set(gca,'YTick',1:5:30,'YTickLabel',1:5:30,'FontSize',14);
xlabel('Time (ms)');ylabel('Tetrode');
title(['ripple = ',num2str(ripnum)]);
box off

subplot(3,2,2); ind_seg1=1:21;
imagesc(1:numSteps,xbins(ind_seg1),postxM_r(ind_seg1,:));
title('postx, clusterless','FontSize',18);
ylabel('Linearized position','FontSize',14);
set(gca,'FontSize',14,'YTick',[-3 -2.5 -2 -1.5 -1 -0.5 0 0.5 1 1.5 2 2.5 3],'YTickLabel',[-3 -2.5 -2 -1.5 -1 -0.5 0 0.5 1 1.5 2 2.5 3]);
colormap(flipud(hot(256)));
caxis([0 0.3]);
xlim([0 numSteps]);

subplot(3,2,4);
imagesc(1:numSteps,xbins(21:31),[postxM_r(21:30,:)+flipud(postxM_r(32:41,:));postxM_r(31,:)]);
ylabel('Linearized position','FontSize',14);
set(gca,'FontSize',14,'YTick',[-3 -2.5 -2 -1.5 -1 -0.5 0 0.5 1 1.5 2 2.5 3],'YTickLabel',[-3 -2.5 -2 -1.5 -1 -0.5 0 0.5 1 1.5 2 2.5 3]);
colormap(flipud(hot(256)));
caxis([0 0.3]);
xlim([0 numSteps]);
ylim([-1 1])

subplot(3,2,6); ind_seg3=41:61;
imagesc(1:numSteps,xbins(ind_seg3),postxM_r(ind_seg3,:));
ylabel('Linearized position','FontSize',14); xlabel('Time (ms)','FontSize',14);
set(gca,'FontSize',14,'YTick',[-3 -2.5 -2 -1.5 -1 -0.5 0 0.5 1 1.5 2 2.5 3],'YTickLabel',[-3 -2.5 -2 -1.5 -1 -0.5 0 0.5 1 1.5 2 2.5 3]);
colormap(flipud(hot(256)));
caxis([0 0.3]);
xlim([0 numSteps]);

pause
close all
end
