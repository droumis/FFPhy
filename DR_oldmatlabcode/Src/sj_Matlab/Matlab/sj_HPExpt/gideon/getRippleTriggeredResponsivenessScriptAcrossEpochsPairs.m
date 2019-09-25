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
% 0 is sleep, 1 is run, -1 is no epoch
epochTypes=[-1*ones(7,10);0 1 0 1 0 1 0 -1 -1 -1;0 1 0 0 1 0 -1 -1 -1 -1;0 1 0 1 0 1 0 -1 -1 -1;repmat([0 1 0 1 0 1 0 1 0 -1],4,1);0 0 1 0 1 0 1 0 1 0;repmat([0 1 0 1 0 1 0 1 0 -1],2,1)];


regions={'AC','HC','PFC'};
valsRes=[];

plotFigures=0;
% Identifying rip-responsive neurons that were recorded simultaneously
% First, get rip-responsive neurons per epoch
% This gives the same number for all (rip-responsive) neurons recorded simultaneously
ripResp=find(ripTrigRespsMat(:,5)==1);
ripRespInds=ripTrigRespsMat(ripResp,1:4);
[xx yy ripRespInd1]=unique(ripRespInds(:,1:2),'rows');
% each kk is a group of neurons recorded simultaneously
for kk=1:max(ripRespInd1)
    close all
    curRipRespInd=ripResp(find(ripRespInd1==kk));
    if length(curRipRespInd)>1
        tetSorts=ceil((ripTrigRespsMat(curRipRespInd,3))/7);
        regionLabels= regions(tetSorts);
        perTrialResps=[];
        perTrialRespsJitter=[];
        allPairs=nchoosek(curRipRespInd,2);
        
        allPairsInd=nchoosek(1:length(curRipRespInd),2);
        
        binProbMat=[];
        corMat=[];
        pMat=[];
        cellLocations=[];
        ensemblePerTrialResps=[];
        for pp=1:length(curRipRespInd)
            ensemblePerTrialResps=[ensemblePerTrialResps ripTrigResps(curRipRespInd(pp)).trialResps ];
        end
        
        if length(allPairs)>1
            for pp=1:size(allPairs,1)
                
                perTrialRespsJitter=[ripTrigRespsJitter(allPairs(pp,1)).trialResps ripTrigRespsJitter(allPairs(pp,2)).trialResps];
                trialsCell1= ripTrigResps(allPairs(pp,1)).trialResps;
                trialsCell2= ripTrigResps(allPairs(pp,2)).trialResps;
                
                % Discretizing responses by converting each spike-count response to 1/0 indicating higher/lower than average
                % trialsCell1D=(trialsCell1-mean(trialsCell1))>0;
                % trialsCell2D=(trialsCell2-mean(trialsCell2))>0;
                
                %
                % P12=mean(trialsCell1D+trialsCell2D==2)/(mean(trialsCell1D)*mean(trialsCell2D))
                % this is an easy extension if we want to have bi-directional
                %     P21=sum(trialsCell1D+trialsCell2D==2)/sum(trialsCell2D)
                
                %binProbMat(allPairsInd(pp,2),allPairsInd(pp,1))=P12;
                
                %binProbMat(allPairsInd(pp,2),allPairsInd(pp,1))=P21;
                
                
                perTrialResps=[ trialsCell1 trialsCell2];
                perTrialRespsShuffled=[trialsCell1(randperm(length(trialsCell1))) trialsCell2(randperm(length(trialsCell2)))];
                
                [RcorValJitter PcorValJitter]=corrcoef(perTrialRespsJitter);
                [RcorVal PcorVal]=corrcoef(perTrialResps);
                [RcorValShuff PcorValShuff]=corrcoef(perTrialRespsShuffled);
                
                
                RcorValJitter=RcorValJitter(2,1);
                RcorVal=RcorVal(2,1);
                RcorValShuff=RcorValShuff(2,1);
                
                PcorValJitter=PcorValJitter(2,1);
                PcorVal=PcorVal(2,1);
                PcorValShuff=PcorValShuff(2,1);
                
                corMat(allPairsInd(pp,2),allPairsInd(pp,1))=RcorVal;
                pMat(allPairsInd(pp,2),allPairsInd(pp,1))=PcorVal;
                
                dayX=ripTrigResps(allPairs(pp,1)).ind(1);
                epochX=ripTrigResps(allPairs(pp,1)).ind(2);
                tet1X=ripTrigResps(allPairs(pp,1)).ind(3);
                cell1X=ripTrigResps(allPairs(pp,1)).ind(4);
                tet2X=ripTrigResps(allPairs(pp,2)).ind(3);
                cell2X=ripTrigResps(allPairs(pp,2)).ind(4);
                sleepOrRun=epochTypes(dayX,epochX);
                
                
                valsRes=[valsRes;dayX epochX tet1X cell1X tet2X cell2X RcorVal PcorVal RcorValJitter PcorValJitter RcorValShuff PcorValShuff sleepOrRun];
            end
        end
        
        numCells=size(corMat,1);
        positions=[cos(0:2*pi/(numCells):2*pi)' sin(0:2*pi/(numCells):2*pi)'];
        
        if 0 % plotting or not
            
            scrsz = get(0,'ScreenSize');
            figure('Position',[scrsz(3)/3 scrsz(4)/8 scrsz(3)/3 scrsz(4)/1.2])
            subplot(2,1,1);
            imagesc(ensemblePerTrialResps)
            subplot(2,1,2)
            for i=1:numCells
                for j=1:i-1
                    curVal=corMat(i,j);
                    curP=pMat(i,j);
                    if curVal>0 & curP<0.05
                        plot([positions(i,1) positions(j,1)],[positions(i,2) positions(j,2)],'linewidth',abs(50*corMat(i,j)),'color','k');
                    elseif curVal<0 & curP<0.05
                        plot([positions(i,1) positions(j,1)],[positions(i,2) positions(j,2)],'linewidth',abs(50*corMat(i,j)),'color','r');
                    else
                        plot([positions(i,1) positions(j,1)],[positions(i,2) positions(j,2)],'--','linewidth',1,'color','k');
                        %
                    end
                    hold on;
                end;
            end
            
            for i=1:numCells
                if positions(i,1)>0
                    text(positions(i,1)*1.5, positions(i,2)*1.5, [regionLabels(i)],'color','b','fontsize',20)
                else
                    text(positions(i,1)*1.5, positions(i,2)*1.5, [regionLabels(i)],'color','b','fontsize',20)
                end
            end
            axis([-2 2 -2 2])
            title(['Day ' num2str(dayX) ' Epoch ' num2str(epochX)])
            keyboard
            
        end
        
    end
end
valsResSleep=valsRes(valsRes(:,end)==0,:);
valsResWake=valsRes(valsRes(:,end)==1,:);
% example, can insert into loop:
%figure;plot(ensemblePerTrialResps(:,1)+0.1*randn(size(ensemblePerTrialResps(:,6))),ensemblePerTrialResps(:,5)+0.1*randn(size(ensemblePerTrialResps(:,5))),'x')
%%

% valsRes columns: 1-day 2-epoch 3-tet1 4-cell1 5-tet2 6-cell2 7-RcorVal 8-PcorVal 9-RcorValJitter 10-PcorValJitter 11-RcorValShuff 12-PcorVaShuff 13-sleepOrRun

% There are significantly more correlated pairs during ripples with real data as compared
% to trial-shuffled data.
corrRateRipples=nanmean(valsRes(:,8)<0.05)
corrRateAfterRipples=nanmean(valsRes(:,10)<0.05)
corrRateTrialShuffled=nanmean(valsRes(:,12)<0.05)

[r p]=ttest(valsRes(:,8)<0.05,valsRes(:,12)<0.05)

% This is also true for after-ripple time window
[r p]=ttest(valsRes(:,10)<0.05,valsRes(:,12)<0.05)

% Significantly more pairs of neurons are correlated during ripples compared to after ripples
[r p]=ttest(valsRes(:,8),valsRes(:,10))

ACACrippleCorrRate=mean(valsRes(:,3)<8 & valsRes(:,5)<8 & valsRes(:,8)<0.05 )/mean(valsRes(:,3)<8 & valsRes(:,5)<8);
ACACafterRippleCorrRate=mean(valsRes(:,3)<8 & valsRes(:,5)<8 & valsRes(:,10)<0.05 )/mean(valsRes(:,3)<8 & valsRes(:,5)<8)
ACACshuffCorrRate=mean(valsRes(:,3)<8 & valsRes(:,5)<8 & valsRes(:,12)<0.05 )/mean(valsRes(:,3)<8 & valsRes(:,5)<8)

PFCPFCrippleCorrRate=mean(valsRes(:,3)>14 & valsRes(:,5)>14 & valsRes(:,8)<0.05 )/mean(valsRes(:,3)>14 & valsRes(:,5)>14)
PFCPFCafterRippleCorrRate=mean(valsRes(:,3)>14 & valsRes(:,5)>14 & valsRes(:,10)<0.05 )/mean(valsRes(:,3)>14 & valsRes(:,5)>14)
PFCPFCshufCorrRate=mean(valsRes(:,3)>14 & valsRes(:,5)>14 & valsRes(:,12)<0.05 )/mean(valsRes(:,3)>14 & valsRes(:,5)>14)

ACPFCrippleCorrRate=mean(valsRes(:,3)<8 & valsRes(:,5)>14 & valsRes(:,8)<0.05 )/mean(valsRes(:,3)<8 & valsRes(:,5)>14)
ACPFCafterRippleCorrRate=mean(valsRes(:,3)<8 & valsRes(:,5)>14 & valsRes(:,10)<0.05 )/mean(valsRes(:,3)<8 & valsRes(:,5)>14)
ACPFCshufCorrRate=mean(valsRes(:,3)<8 & valsRes(:,5)>14 & valsRes(:,12)<0.05 )/mean(valsRes(:,3)<8 & valsRes(:,5)>14)

ACHCrippleCorrRate=mean(valsRes(:,3)<8 & valsRes(:,5)>7 & valsRes(:,5)<15 & valsRes(:,8)<0.05 )/mean(valsRes(:,3)<8 & valsRes(:,5)>7 & valsRes(:,5)<15)
ACHCafterRippleCorrRate=mean(valsRes(:,3)<8 & valsRes(:,5)>7 & valsRes(:,5)<15 & valsRes(:,10)<0.05 )/mean(valsRes(:,3)<8 & valsRes(:,5)>7 & valsRes(:,5)<15)
ACHCshufCorrRate=mean(valsRes(:,3)<8 & valsRes(:,5)>7 & valsRes(:,5)<15 & valsRes(:,12)<0.05 )/mean(valsRes(:,3)<8 & valsRes(:,5)>7 & valsRes(:,5)<15)

HCPFCrippleCorrRate=mean(valsRes(:,3)>7 & valsRes(:,3)<15 & valsRes(:,5)>14 & valsRes(:,8)<0.05 )/mean(valsRes(:,3)>7 & valsRes(:,3)<15 & valsRes(:,5)>14 )
HCPFCafterRippleCorrRate=mean(valsRes(:,3)>7 & valsRes(:,3)<15 & valsRes(:,5)>14 & valsRes(:,10)<0.05 )/mean(valsRes(:,3)>7 & valsRes(:,3)<15 & valsRes(:,5)>14 )
HCPFCshufCorrRate=mean(valsRes(:,3)>7 & valsRes(:,3)<15 & valsRes(:,5)>14 & valsRes(:,12)<0.05 )/mean(valsRes(:,3)>7 & valsRes(:,3)<15 & valsRes(:,5)>14 )

HCPFCrippleCorrRate=mean(valsRes(:,3)>7 & valsRes(:,3)<15 & valsRes(:,5)>7 & valsRes(:,5)<15 & valsRes(:,8)<0.05 )/mean(valsRes(:,3)>7 & valsRes(:,3)<15 & valsRes(:,5)>7 & valsRes(:,5)<15 )
HCPFCafterRippleCorrRate=mean(valsRes(:,3)>7 & valsRes(:,3)<15 & valsRes(:,5)>7 & valsRes(:,5)<15 & valsRes(:,10)<0.05 )/mean(valsRes(:,3)>7 & valsRes(:,3)<15 & valsRes(:,5)>7 & valsRes(:,5)<15 )
HCPFCshufCorrRate=mean(valsRes(:,3)>7 & valsRes(:,3)<15 & valsRes(:,5)>7 & valsRes(:,5)<15 & valsRes(:,12)<0.05 )/mean(valsRes(:,3)>7 & valsRes(:,3)<15 & valsRes(:,5)>7 & valsRes(:,5)<15 )

ACACcorrVals=valsRes(valsRes(:,3)<8 & valsRes(:,5)<8 & valsRes(:,8)<0.05 ,7 );
PFCPFCcorrVals=valsRes(valsRes(:,3)>14 & valsRes(:,5)>14 & valsRes(:,8)<0.05 ,7 );
ACPFCcorrVals=valsRes(valsRes(:,3)<8 & valsRes(:,5)>14 & valsRes(:,8)<0.05 ,7 );
HCPFCcorrVals=valsRes(valsRes(:,3)>7 & valsRes(:,3)<15 & valsRes(:,5)>14 & valsRes(:,8)<0.05 ,7 );
HCACcorrVals=valsRes(valsRes(:,5)>7 & valsRes(:,5)<15 & valsRes(:,3)<8 & valsRes(:,8)<0.05 ,7 );
HCHCcorrVals=valsRes(valsRes(:,5)>7 & valsRes(:,5)<15 & valsRes(:,3)>7 & valsRes(:,3)<15 & valsRes(:,8)<0.05 ,7 );
%% Taking the std of correlations of the same pair across different epochs, and comparing to shuffled epochs
% std of the same pair across epochs is significantly lower than when
% choosing 3 epochs randomly.
% See in next sessions that this becomes even more significant if
% considering only sleep. This is not true for only wake, maybe because
% there are few cases where there are 3 or more wake epochs for same
% pair. But when I lower the threshold to 2 epochs, the trend is the
% opposite- that corr values change more in consecutive wake epochs in
% comparison to random

% This line finds rows (=pairs of cells) recorded on same day, tetrodes and cells, across different epochs
[xx yy acrossEpochInds2]=unique(valsRes(:,[1 3:6]),'rows');
realstds=[];
randstds=[];
for kk=1:max(acrossEpochInds2)
    curPairs=valsRes(find(acrossEpochInds2==kk),:)
    numepochsK=size(curPairs,1)
    if numepochsK>=3
        realstds=[realstds std(curPairs(:,7))];
        
        randomeps=valsRes(round(rand(1,numepochsK)*size(valsRes,1)),:);
        randstds=[randstds std(randomeps(:,7))];
    end
end
nanmean(randstds)
nanmean(realstds)
[r p]=ttest(randstds,realstds)
%% same for only sleep epochs
[xx yy acrossEpochInds2]=unique(valsResSleep(:,[1 3:6]),'rows');
realstds=[];
randstds=[];
for kk=1:max(acrossEpochInds2)
    curPairs=valsResSleep(find(acrossEpochInds2==kk),:)
    numepochsK=size(curPairs,1)
    if numepochsK>=3
        realstds=[realstds std(curPairs(:,7))];
        
        randomeps=valsRes(round(rand(1,numepochsK)*size(valsRes,1)),:);
        randstds=[randstds std(randomeps(:,7))];
    end
end
nanmean(randstds)
nanmean(realstds)
[r p]=ttest(randstds,realstds)

%% same for only wake epochs
[xx yy acrossEpochInds2]=unique(valsResWake(:,[1 3:6]),'rows');
realstds=[];
randstds=[];
for kk=1:max(acrossEpochInds2)
    curPairs=valsResWake(find(acrossEpochInds2==kk),:)
    numepochsK=size(curPairs,1)
    if numepochsK>=3
        realstds=[realstds std(curPairs(:,7))]
        
        randomeps=valsRes(round(rand(1,numepochsK)*size(valsRes,1)),:);
        randstds=[randstds std(randomeps(:,7))]
    end
end
nanmean(randstds)
nanmean(realstds)
[r p]=ttest(randstds,realstds)
%%
% example of AC neuron positively noise correlated with one place cell and
% negatively with another
xmplInd1=find(ismember(ripTrigRespsMat(:,1:4),[16 3 1 1],'rows'));
xmplInd2=find(ismember(ripTrigRespsMat(:,1:4),[16 3 11 1],'rows'));
xmplInd3=find(ismember(ripTrigRespsMat(:,1:4),[16 3 14 1],'rows'));
figure;plot(ripTrigResps(xmplInd1).trialResps)
hold on
plot(ripTrigResps(xmplInd2).trialResps,'r')
plot(ripTrigResps(xmplInd3).trialResps,'g')

