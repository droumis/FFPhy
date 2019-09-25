% % ran first with 0 (real), then jitter=1 which mean looking at spikes one
% % second after ripples.
%
% for jitter=[0 1]
%     ripTrigResps=[];
%
% for day=8:17
%     for epoch=1:10
%         for tet=1:21
%             for cell=1:8
%
%                     try
%                     ['day ' num2str(day) ' epoch ' num2str(epoch) ' tet ' num2str(tet) ' cell ' num2str(cell)]
%                     ripTrigResp=getRippleTriggeredResponsiveness('Ndl', day,epoch, [9 10 11 13 14], tet, cell,[],[],jitter);
%                     ripTrigResps=[ripTrigResps ripTrigResp];
%                 catch
%                 end
%         end
%       end
%     end
% end
%
% if jitter==0
%     save ripTrigRespsReal ripTrigResps
% else
%         save ripTrigRespsJitter ripTrigResps
%
% end
% end

%save ripTrigResps ripTrigResps

%% calculating rates of rip-responsive neurons per region

load ripTrigRespsJitter
ripTrigRespsJitter=ripTrigResps;
load ripTrigRespsReal

calcjitter=1;


ripTrigRespsMat=[];
for i=1:length(ripTrigResps)
    tmp=[ ripTrigResps(i).ind ripTrigResps(i).h ripTrigResps(i).p];
    ripTrigRespsMat=[ripTrigRespsMat;tmp];
end



acrossEpochInds=unique(ripTrigRespsMat(:,[1,3,4]),'rows');
acrossEpochInds=[acrossEpochInds zeros(length(acrossEpochInds),1)];
for i=1:length(acrossEpochInds)
    % Lines in which this day-tet-cell appeared
    curInds=find(ismember(ripTrigRespsMat(:,[1,3,4]),acrossEpochInds(i,1:3),'rows'));
    fullCurInds=ripTrigRespsMat(curInds,:);
    % current responsiveness criteria: if the cell is significantly
    % responsive in more than 20% of epochs
    if nanmean(fullCurInds(:,5))>0.2
        acrossEpochInds(i,end)=1;
    else
        %      keyboard
        
    end
end

% now acrossEpochInds has: day,tet,cell,0/1 where 0/1 indicates if this
% neuron was ripple-responsive in at least one epoch during that day

responsiveness=acrossEpochInds(:,4);
tetrodes=acrossEpochInds(:,2);
rateACripresponsiveneuorons=nanmean(responsiveness(tetrodes<8));
ratePFCripresponsiveneuorons=nanmean(responsiveness(tetrodes>14));
rateCA1ripresponsiveneuorons=nanmean(responsiveness(tetrodes>7&tetrodes<15));


%% calculating rip-evoked noise correlation
regions={'AC','HC','PFC'}
valsRes=[];
valsResJ=[];

plotFigures=0;
% Identifying rip-responsive neurons that were recorded simultaneously
% First, get rip-responsive neurons per epoch
% This gives the same number for all (rip-responsive) neurons recorded simultaneously
ripResp=find(ripTrigRespsMat(:,5)==1);
ripRespInds=ripTrigRespsMat(ripResp,1:4);
[xx yy ripRespInd1]=unique(ripRespInds(:,1:2),'rows');
% each kk is a group of neurons recorded simultaneously
for kk=1:max(ripRespInd1)
    curRipRespInd=ripResp(find(ripRespInd1==kk));
    perTrialResps=[];
    perTrialRespsJitter=[];
    
    for pp=1:length(curRipRespInd)
        
        perTrialRespsJitter=[perTrialRespsJitter ripTrigRespsJitter(curRipRespInd(pp)).trialResps ];
        
        perTrialResps=[perTrialResps ripTrigResps(curRipRespInd(pp)).trialResps ];
        
    end
    numCells=size(perTrialResps,2);
    
    
    tetSorts=ceil((ripTrigRespsMat(curRipRespInd,3))/7);
    
    %plotting correlations of ripple-triggered activity of the different
    %cells
    plotFigures=0;
    if plotFigures==1
        scrsz = get(0,'ScreenSize');
        figure('Position',[scrsz(3)/3 scrsz(4)/8 scrsz(3)/3 scrsz(4)/1.2])
        
        subplot(4,1,1:3)
        imagesc(perTrialResps)
        subplot(4,1,4)
        imagesc(corrcoef(perTrialResps));
        caxis([-1 1])
        colorbar
        set(gca,'YTick',1:numCells)
        set(gca,'XTick',1:numCells)
        set(gca,'XTickLabel',regions(tetSorts))
        set(gca,'YTickLabel',regions(tetSorts))
        daystr=num2str(ripTrigRespsMat(curRipRespInd(1),1));
        epochstr=num2str(ripTrigRespsMat(curRipRespInd(1),2));
        title(['Day ' daystr ' Epoch ' epochstr])
        figure
        % plotting corr patterns for parts of the data to assess stability
%         numTrigsX=size(perTrialResps,1);
%         trialRes=100;
%         nd=floor(numTrigsX/trialRes);
%         for i=1:nd
%             subplot(nd,1,i);
%             imagesc(corrcoef(perTrialResps((i-1)*trialRes+1:i*trialRes,:)));
%         end
        keyboard
    end
    day1=ripTrigRespsMat(curRipRespInd(1),1);
    epoch1=ripTrigRespsMat(curRipRespInd(1),2);
    
    m=ones(numCells);
    m2=tril(m,-1);
    c=corrcoef(perTrialResps);
    cj=corrcoef(perTrialRespsJitter);
    vals=c(m2==1);
    valsJ=cj(m2==1);
    [x y]=find(m2==1);
    
    % Keeps corr values in the following form: 1/2/3 1/2/3 val,
    % where 1=AC, 2=HC, 3=PFC
    valsRes=[valsRes;tetSorts(x) tetSorts(y) vals];
    valsResJ=[valsResJ;tetSorts(x) tetSorts(y) valsJ];
    
end
%% plot distribution of correlation values per region pairs

ACACvals=valsRes(valsRes(:,1)==1&valsRes(:,2)==1,3);
HCHCvals=valsRes(valsRes(:,1)==2&valsRes(:,2)==2,3);
HCACvals=valsRes(valsRes(:,1)==2&valsRes(:,2)==1,3);
HCPFCvals=valsRes(valsRes(:,1)==3&valsRes(:,2)==2,3);
PFCACvals=valsRes(valsRes(:,1)==3&valsRes(:,2)==1,3);
PFCPFCvals=valsRes(valsRes(:,1)==3&valsRes(:,2)==3,3);

histres=-1:0.05:1;
[xACAC y]=hist(ACACvals,histres);
[xHCHC y]=hist(HCHCvals,histres);
[xHCAC y]=hist(HCACvals,histres);
[xHCPFC y]=hist(HCPFCvals,histres);
[xPFCAC y]=hist(PFCACvals,histres);
[xPFCPFC y]=hist(PFCPFCvals,histres);

figure;
subplot(2,1,1)
plot(y,xACAC./sum(xACAC),'linewidth',3)
hold on
plot(y,xPFCAC./sum(xPFCAC),'k','linewidth',3)
plot(y,xPFCPFC./sum(xPFCPFC),'g','linewidth',3)
plot(y,xHCAC./sum(xHCAC),'m','linewidth',3)
plot(y,xHCPFC./sum(xHCPFC),'r','linewidth',3)
plot(y,xHCHC./sum(xHCHC),'k--','linewidth',3)
xlabel('Correlation')
ylabel('Probability')

legend('AC-AC','PFC-AC','PFC-PFC','HC-AC','HC-PFC','HC-HC')
title('Real data')
% jitter

ACACvalsJ=valsResJ(valsRes(:,1)==1&valsRes(:,2)==1,3);
HCHCvalsJ=valsResJ(valsRes(:,1)==2&valsRes(:,2)==2,3);
HCACvalsJ=valsResJ(valsRes(:,1)==2&valsRes(:,2)==1,3);
HCPFCvalsJ=valsResJ(valsRes(:,1)==3&valsRes(:,2)==2,3);
PFCACvalsJ=valsResJ(valsRes(:,1)==3&valsRes(:,2)==1,3);
PFCPFCvalsJ=valsResJ(valsRes(:,1)==3&valsRes(:,2)==3,3);

[xACACJ yj]=hist(ACACvalsJ,histres);
[xHCHCJ yj]=hist(HCHCvalsJ,histres);
[xHCACJ yj]=hist(HCACvalsJ,histres);
[xHCPFCJ yj]=hist(HCPFCvalsJ,histres);
[xPFCACJ yj]=hist(PFCACvalsJ,histres);
[xPFCPFCJ yj]=hist(PFCPFCvalsJ,histres);

subplot(2,1,2)
plot(y,xACACJ./sum(xACACJ),'linewidth',3)
hold on
plot(y,xPFCACJ./sum(xPFCACJ),'k','linewidth',3)
plot(y,xPFCPFCJ./sum(xPFCPFCJ),'g','linewidth',3)
plot(y,xHCACJ./sum(xHCACJ),'m','linewidth',3)
plot(y,xHCPFCJ./sum(xHCPFCJ),'r','linewidth',3)
plot(y,xHCHCJ./sum(xHCHCJ),'k--','linewidth',3)
xlabel('Correlation')
ylabel('Probability')

legend('AC-AC','PFC-AC','PFC-PFC','HC-AC','HC-PFC','HC-HC')
title('jittered')

figure;

subplot(3,2,1)
plot(ACACvalsJ,ACACvals,'ok','markerfacecolor','k')
hold on;plot(histres,histres)
xlabel('Corr during jittered')
ylabel('Corr during ripples')
title('AC-AC')

subplot(3,2,2)
plot(PFCACvalsJ,PFCACvals,'ok','markerfacecolor','k')
hold on;plot(histres,histres)
xlabel('Corr during jittered')
ylabel('Corr during ripples')
title('PFC-AC')


subplot(3,2,3)
plot(PFCPFCvalsJ,PFCPFCvals,'ok','markerfacecolor','k')
hold on;plot(histres,histres)
xlabel('Corr during jittered')
ylabel('Corr during ripples')
title('PFC-PFC')



subplot(3,2,4)
plot(HCPFCvalsJ,HCPFCvals,'ok','markerfacecolor','k')
hold on;plot(histres,histres)
xlabel('Corr during jittered')
ylabel('Corr during ripples')
title('HC-PFC')



subplot(3,2,5)
plot(HCACvalsJ,HCACvals,'ok','markerfacecolor','k')
hold on;plot(histres,histres)
xlabel('Corr during jittered')
ylabel('Corr during ripples')
title('HC-AC')



subplot(3,2,6)
plot(HCHCvalsJ,HCHCvals,'ok','markerfacecolor','k')
hold on;plot(histres,histres)
xlabel('Corr during jittered')
ylabel('Corr during ripples')
title('HC-HC')
