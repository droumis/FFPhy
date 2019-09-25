% Ver4 : Starting 10Feb2014 - Sync codes with everyone

% GLM model for theta times. Similar to DFSsj_HPexpt_glm_ripalign.m
% Use initial filter settings from "DFSsj_HPexpt_ThetacorrAndRipresp_ver3.m"
% In the function DFAsj_glm_theta.m (based on DFAsj_glm_ripalign), get the time series with theta times
% similar to the beginnin of the function"DFAsj_getthetacrosscov_timecondition.m"

% -----------------
% Glm model for PFC rip-align response as a function of rip-aligned CA1 popln responses
% Goal is to get coefficients and significance
% Get ripple aligned responses directly from ripplemod structure
% Also get pairwise corrlns for now, without shuffle sign for speed. You cal always compare to saved apirwise corrln files later


clear; %close all;
runscript = 1;
savedata = 1; % save data option - only works if runscript is also on
figopt1 = 0; % Figure Options - Individual cells
plotGraphs=1;

savedir = '/data25/sjadhav/HPExpt/ProcessedData/';
%savedir = '/data15/gideon/ProcessedData/';
[y, m, d] = datevec(date);

val=1; savefile = [savedir 'HP_ripmod_glmfit_theta',num2str(m),'-',num2str(d),'-',num2str(y)]; area = 'PFC'; clr = 'b'; % PFC


savefig1=0;


% Plot options
plotanimidx =  []; % To pick animals for plotting
plotdays = []; % If you only load data when runscript=0 and savedata=0, then this field will supplant days


% If runscript, run Datafilter and save data
if runscript == 1
    
    %Animal selection
    %-----------------------------------------------------
    animals = {'HPa','HPb','HPc','Ndl'};
    %animals = {'HPa','HPb','HPc','Nadal'};
     %   animals = {'Nadal'};

    %Filter creation
    %-----------------------------------------------------
    
    % Epoch filter
    % -------------
    % dayfilter = '1:8'; % Shantanu - I am adding day filter to parse out epoch filter
    % Either Only do 1st w-track. 2 or 1 epochs per day
    % Or do Wtr1 and Wtr2, 2 epochs per day
    runepochfilter = 'isequal($environment, ''wtr1'') || isequal($environment, ''wtr2'') || isequal($environment, ''ytr'')';
    %sleepepochfilter = 'isequal($type, ''sleep'')'; % Only pre and post sleep marked as sleep
    
    % %Cell filter
    % %-----------
    % %PFC
    % %----
    % && strcmp($thetamodtag, ''y'')
    switch val
        case 1
            % All Ca1 cells and PFC Ripple modulated cells. Function will parse them out.
            cellfilter = '((strcmp($area, ''CA1'') || strcmp($area, ''iCA1'')) && ($numspikes > 100) && ($meanrate < 7))  ||  (strcmp($area, ''PFC'') && ($numspikes > 100) && strcmp($ripmodtag2, ''y'')) ';
    end
    
    % Time filter - none.
    % -----------
    riptetfilter = '(isequal($descrip, ''riptet''))';
    
    % Use absvel instead of linearvel
    timefilter_place_new = { {'DFTFsj_getvelpos', '(($absvel >= 5))'},...
        {'DFTFsj_getriptimes','($nripples == 0)','tetfilter',riptetfilter,'minthresh',3} };
    
    % Iterator
    % --------
    iterator = 'multicellanal';
    
    % Filter creation
    % ----------------
    modf = createfilter('animal',animals,'epochs',runepochfilter, 'cells',...
        cellfilter,  'excludetime', timefilter_place_new, 'iterator', iterator);
    
    modg = createfilter('animal',animals,'epochs',runepochfilter, 'cells',...
        cellfilter, 'iterator', iterator);
    
    disp('Done Filter Creation');
    
    % Set analysis function
    % ----------------------
    
    % % Need spikes to get time series
    % % If it is across regions, you will need cellinfo to get area where cells are recorded from
    modg = setfilterfunction(modg,'DFAsj_glm_ripalign_dataForTheta4',{'ripplemod','cellinfo'}); %
    modf = setfilterfunction(modf,'DFAsj_glm_thetaGR4',{'spikes','cellinfo'},'thrstime',1); % With includetime condition
    
    
    % Run analysis
    % ------------
    modg = runfilter(modg);
    modf = runfilter(modf);
    disp('Finished running filter script');
    %--------------------- Finished Filter Function Run -------------------
    
    if savedata == 1
        clear figopt1 runscript plotdays plotanimidx savedata
        save(savefile);
    end
    
else
    
    load(savefile);
    
end % end runscript

if ~exist('savedata')
    return
end


% -------------------------  Filter Format Done -------------------------

% ----------------------------------
% Whether to gather data or to load previously gathered data
% --------------------------------------------------------------------
gatherdata = 1; savegatherdata = 1;
[y, m, d] = datevec(date);

switch val
    case 1
        gatherdatafile = [savedir 'HP_ripmod_glmfit_theta_gather',num2str(m),'-',num2str(d),'-',num2str(y)];
end



if gatherdata
    %----- Getting theta model parameters
    cnt=0;
    %Glm
    allglmidxstheta=[];    allglmidxstheta2=[];
    allmodelbtheta=[]; allmodelbtheta2=[];allmodelptheta=[]; allmodelfitstheta=[];
    allnsigtheta=[]; allnsigpostheta=[]; allnsignegtheta=[];  allfracsigtheta=[]; allfracsigpostheta=[]; allfracsignegtheta=[];
    
    for an = 1:length(modf)
        for i=1:length(modf(an).output{1})
            % Check for empty output - If Cell defined in rpoch and Nspks in ripple response wndow > 0
            if ~isempty(modf(an).output{1}(i).nsig)
                cnt=cnt+1;
                anim_index{an}{cnt} = modf(an).output{1}(i).indices;
                % Only indexes
                %animindex=[an modf(an).output{1}(i).indices]; % Put animal index in front
                %allanimindex = [allanimindex; animindex]; % Collect all Anim Day Epoch Tet Cell Index
                %Sub indices for glm and corrcoef
                allglmidxstheta = [allglmidxstheta; modf(an).output{1}(i).glmidxs]; % Need to put "an" in there
                
                indNoAnim=modf(an).output{1}(i).glmidxs2;
                indNoAnim=indNoAnim(find(nansum(indNoAnim')~=0),:);
                indWAnim=[an*ones(size(indNoAnim,1),1) indNoAnim];
                allglmidxstheta2 = [allglmidxstheta2; indWAnim];
                
                
                % Glm data
                allmodelbtheta = [allmodelbtheta; modf(an).output{1}(i).allmodelb];
                allmodelbtheta2 = [allmodelbtheta2; modf(an).output{1}(i).allmodelb2];
                allmodelptheta = [allmodelptheta; modf(an).output{1}(i).allmodelp];
                allmodelfitstheta{cnt} = modf(an).output{1}(i).allmodelfits;
                allnsigtheta = [allnsigtheta; modf(an).output{1}(i).nsig'];
                allnsigpostheta = [allnsigpostheta; modf(an).output{1}(i).nsigpos'];
                allnsignegtheta = [allnsignegtheta; modf(an).output{1}(i).nsigneg'];
                allfracsigtheta = [allfracsigtheta; modf(an).output{1}(i).fracsig'];
                allfracsigpostheta = [allfracsigpostheta; modf(an).output{1}(i).fracsigpos'];
                allfracsignegtheta = [allfracsignegtheta; modf(an).output{1}(i).fracsigneg'];
                % Corr data
                
            end
        end
        
    end
    %removing rows and columns that are all nan
    allglmidxstheta2=allglmidxstheta2(find(nansum(allglmidxstheta2')~=0),1:find(nansum(allglmidxstheta2)==0,1));
    allmodelbtheta2=allmodelbtheta2(find(nansum(allmodelbtheta2')~=0),1:find(nansum(allmodelbtheta2)==0,1));
    
    cnt=0;
    %---- Getting ripple model data
    allglmidxsrip2=[];
    %     allmodelbtheta=[]; allmodelbtheta2=[];allmodelptheta=[]; allmodelfitstheta=[];
    %     allnsigtheta=[]; allnsigpostheta=[]; allnsignegtheta=[];  allfracsigtheta=[]; allfracsigpostheta=[]; allfracsignegtheta=[];
    Xmats={};
    Ymats={};
    for an = 1:length(modg)
        for i=1:length(modg(an).output{1})
            % Check for empty output - If Cell defined in rpoch and Nspks in ripple response wndow > 0
            %    if ~isempty(modg(an).output{1}(i).nsig)
            cnt=cnt+1;
            anim_index{an}{cnt} = modf(an).output{1}(i).indices;
            % Only indexes
            %animindex=[an modf(an).output{1}(i).indices]; % Put animal index in front
            %allanimindex = [allanimindex; animindex]; % Collect all Anim Day Epoch Tet Cell Index
            %Sub indices for glm and corrcoef
            % allglmidxstheta = [allglmidxstheta; modf(an).output{1}(i).glmidxs]; % Need to put "an" in there
            
            
            indNoAnim=modg(an).output{1}(i).glmidxs2;
            indNoAnim=indNoAnim(find(nansum(indNoAnim')~=0),:);
            indWAnim=[an*ones(size(indNoAnim,1),1) indNoAnim];
            % Indices in new form. Each line is one ensemble of cells in the
            % following form:
            % animal day,epoch,hctet,hccell,hctet,hccell...,hctet,hccell,pfctet,pfccell
            allglmidxsrip2 = [allglmidxsrip2; indWAnim]; % Need to put "an" in there
            Xmats{cnt}=modg(an).output{1}(i).Xmat;
            Ymats{cnt}=modg(an).output{1}(i).Ymat;
            
            
            %    end
        end
        
    end
    
    allglmidxsrip2=allglmidxsrip2(find(nansum(allglmidxsrip2')~=0),1:find(nansum(allglmidxsrip2)==0,1));
    
    % sanity check
    numPFCCells=0;for i=1:size(Ymats,2);numPFCCells=numPFCCells+size(Ymats{i},2);end
    if numPFCCells~=size(allglmidxsrip2)|numPFCCells~=size(allglmidxstheta2)
        'Numbers dont add up'
        keyboard
    end
    allErrReal_theta=[];
    allErrShuf_theta=[];
    allErrReal=[];
    allErrShuf=[];
    allPs_theta=[];
    allPs=[];
    allPsK=[];
    nsig=[];
    fracsig=[];
    corrPlast=[];
    counter=1;
    allRealPlast=[];
    allShufPlast=[];
    figdir = '/data15/gideon/Figs/';
    plotGraphs=1;
    for ii=1:size(Xmats,2)
        currY=Ymats{ii};
        currX=Xmats{ii};
        %plotting graphs
        nPFCcells=size(currY,2);
        nCA1cells=size(currX,2);
        if plotGraphs
            [rr pp]=corrcoef([currX currY]);
            regionsLabels=[repmat('CA1',nCA1cells,1);repmat('PFC',nPFCcells,1)];
            % plotting only PFC-CA1 edges (removing CA1-CA1 and PFC-PFC)
            for i=1:length(rr),for j=1:length(rr), if i<=nCA1cells&j<=nCA1cells,rr(i,j)=NaN;end,if i>nCA1cells&j>nCA1cells,rr(i,j)=NaN;end;end;end
            plotNiceCorrGraph(rr,pp,regionsLabels)
            animStr=num2str(allglmidxsrip2(counter,1));
         
            dayStr=num2str(allglmidxsrip2(counter,2));
            epochStr=num2str(allglmidxsrip2(counter,3));
            title(['Rip- corrs, anim=' animStr ' day=' dayStr ' ep=' epochStr])
            saveas(gcf,[figdir 'graph' animStr dayStr epochStr 'B'],'jpg')
            close all
        end
        for jj=1:size(currY,2)
            y=currY(:,jj);
            %-----------
            [btrall, ~, statsall] = glmfit(currX,y,'poisson');
            currsig = find(statsall.p(2:end) < 0.05);
            nsig = [nsig length(currsig)];
            fracsig = [fracsig length(currsig)/nCA1cells];
            
            numRips=size(currX,1);
            
            allErrReal1=[];
            allErrShuf1=[];
            for ii=1:1000
                ripidxs=randperm(numRips);
                dataPercentForTrain=0.9;
                Ntrain=ripidxs(1:round(numRips*dataPercentForTrain));
                Ntest=ripidxs(round(numRips*dataPercentForTrain)+1:numRips);
                [btr, ~, statstr] = glmfit(currX(Ntrain,:),y(Ntrain),'poisson');
                % btr(statstr.p > 0.05)=0;
                yfit = glmval(btr, currX(Ntest,:),'log',statstr,0.95);
                Ntestshufd=Ntest(randperm(length(Ntest)));
                yfit(yfit>10)=NaN;
                errReal=nanmean(abs(yfit-y(Ntest)));
                errShuf=nanmean(abs(yfit-y(Ntestshufd)));
                
                allErrReal1=[allErrReal1 errReal];
                allErrShuf1=[allErrShuf1 errShuf];
            end
            [r1,p1]=ttest2(allErrReal1,allErrShuf1,0.05,'left');
            [kp1]=kruskalwallis([allErrReal1' allErrShuf1'],[],'off');
            allErrReal=[allErrReal nanmean(allErrReal1)];
            allErrShuf=[allErrShuf nanmean(allErrShuf1)];
            
            allPs=[allPs p1];
            allPsK=[allPsK kp1];
            %--------------
            currB=allmodelbtheta2(counter,1:find(~isnan(allmodelbtheta2(counter,:)),1,'last'));
            yfit_theta = glmval(currB', currX,'log');
            yfit_theta(yfit_theta>10)=NaN;
            errReal_theta=nanmean(abs(yfit_theta-y));
            allErrShufTmp=[];
            for pp=1:1000
                errShuf=nanmean(abs(yfit_theta-y(randperm(length(y)))));
                allErrShufTmp=[allErrShufTmp errShuf];
            end
            currP=nanmean(errReal_theta>allErrShufTmp);
            allPs_theta=[allPs_theta currP];
            allErrReal_theta=[allErrReal_theta errReal_theta];
            allErrShuf_theta=[allErrShuf_theta nanmean(allErrShufTmp)];
            counter=counter+1;
        end
    end
    
    % ----------
    % No combination across epochs
    % ---------
    % Save
    % -----
    if savegatherdata == 1
        save(gatherdatafile);
    end
    
else % gatherdata=0
    
    load(gatherdatafile);
    
end % end gather data


% ------------
% PLOTTING, ETC
% ------------


fractionSigPredictable=mean(allPs<0.05)
fractionSigPredictableTheta=mean(allPs_theta<0.05)
errDecrease=(allErrReal./allErrShuf);
errDecrease_theta=(allErrReal_theta./allErrShuf_theta);


t=[];for i=1:4,t=[t mean(allPs(nsig<=(i-1))<0.05)];end
terr=[];for i=1:4,terr=[terr std(allPs(nsig<=(i-1))<0.05)];end
tn=[];for i=1:4,tn=[tn sum(allPs(nsig<=(i-1))<0.05)];end
figure;errorbar(0:3,t,terr./sqrt(tn-1),'k','linewidth',2)
xlabel('Maximal number of significant CA1 units');
ylabel('Fraction of predictable PFC cells')

hold on

t=[];for i=1:4,t=[t mean(allPs_theta(allnsigtheta<=(i-1))<0.05)];end
terr=[];for i=1:4,terr=[terr std(allPs_theta(allnsigtheta<=(i-1))<0.05)];end
tn=[];for i=1:4,tn=[tn sum(allPs_theta(allnsigtheta<=(i-1))<0.05)];end
errorbar(0:3,t,terr./sqrt(tn-1),'r','linewidth',2)
xlabel('Maximal number of significant CA1 units');
ylabel('Fraction of predictable PFC cells')


t2=[];for i=1:4,t2=[t2 nanmean(errDecrease(nsig<=(i-1)))];end
terr2=[];for i=1:4,terr2=[terr2 nanstd(errDecrease(nsig<=(i-1)))];end
tn2=[];for i=1:4,tn2=[tn2 sum(errDecrease(nsig<=(i-1)))];end
figure;errorbar(0:3,t2,terr2./sqrt(tn2),'k','linewidth',2)
xlabel('Maximal number of significant CA1 units');
ylabel('Error relative to shuf')
hold on

t2=[];for i=1:4,t2=[t2 mean(errDecrease_theta(nsig<=(i-1)))];end
terr2=[];for i=1:4,terr2=[terr2 std(errDecrease_theta(nsig<=(i-1)))];end
tn2=[];for i=1:4,tn2=[tn2 sum(errDecrease_theta(nsig<=(i-1)))];end
errorbar(0:3,t2,terr2./sqrt(tn2),'r','linewidth',2)
xlabel('Maximal number of significant CA1 units');
ylabel('Error relative to shuf')
% ------------------
% Population Figures
% ------------------

% Idxs for allidxs and corrodxs are exactly matched. So you can just take values from glmfits and corrcef directly

% Skip bad fits, and corrlns where no. of co-occurences are <10
rem = find( ((allmodelb>1) | (allmodelb<-1)) & allmodelp>0.99); % Corresponding p will usually be 1
%rem = find(allmodelp ==1);
rem2 = find(allnsimul<10);
allrem = union(rem, rem2);
allmodelb(allrem)=[]; allmodelp(allrem)=[]; allrcorr(allrem)=[]; allpcorr(allrem)=[]; allglmidxs(allrem,:)=[]; allcorridxs(allrem,:)=[];
sigglm = find(allmodelp < 0.05);
sigcorr = find(allpcorr < 0.05);

% Sig GLM vs Sig Corr
glmvec = zeros(size(allmodelb)); glmvec(sigglm)=1;
corrvec = zeros(size(allmodelb)); corrvec(sigcorr)=1;
[rvec,pvec] = corrcoef(glmvec,corrvec)


forppr = 0;
% If yes, everything set to redimscreen_figforppr1
% If not, everything set to redimscreen_figforppt1
figdir = '/data25/sjadhav/HPExpt/Figures/ThetaMod/';
summdir = figdir;
set(0,'defaultaxesfontweight','normal'); set(0,'defaultaxeslinewidth',2);

if forppr==1
    set(0,'defaultaxesfontsize',16);
    tfont = 18; % title font
    xfont = 16;
    yfont = 16;
else
    set(0,'defaultaxesfontsize',24);
    tfont = 28;
    xfont = 20;
    yfont = 20;
end


figure; hold on;redimscreen_figforppt1;
plot(allrcorr, allmodelb, 'k.','MarkerSize',24);
plot(allrcorr(sigglm), allmodelb(sigglm), 'r.','MarkerSize',24);
plot(allrcorr(sigcorr), allmodelb(sigcorr), 'bo','MarkerSize',20);
title(sprintf('GLM fits vs Corr Coeff'),'FontSize',24,'Fontweight','normal');
xlabel(['Corr Coeff'],'FontSize',24,'Fontweight','normal');
ylabel(['GLM coeffs'],'FontSize',24,'Fontweight','normal');
legend('All Pairs','Sig GLM','Sig Corr');
text(-0.22,0.4,['Npairs:' num2str(length(allmodelb))],'FontSize',24,'Fontweight','normal');


[rall,pall] = corrcoef(allrcorr,allmodelb)
[rglm,pglm] = corrcoef(allrcorr(sigglm),allmodelb(sigglm))
[rc,pc] = corrcoef(allrcorr(sigcorr),allmodelb(sigcorr))

% figure; hold on;redimscreen_figforppt1;
% plot(abs(allrcorr), abs(allmodelb), 'k.','MarkerSize',24);
% plot(abs(allrcorr(sigglm)), abs(allmodelb(sigglm)), 'r.','MarkerSize',24);
% plot(abs(allrcorr(sigcorr)), abs(allmodelb(sigcorr)), 'bo','MarkerSize',20);
% title(sprintf('ABS GLM fits vs Corr Coeff'),'FontSize',24,'Fontweight','normal');
% xlabel(['Corr Coeff'],'FontSize',24,'Fontweight','normal');
% ylabel(['GLM coeffs'],'FontSize',24,'Fontweight','normal');
% legend('All Pairs','Sig GLM','Sig Corr');



figure; hold on;%redimscreen_figforppt1;
subplot(3,1,1)
h = bar(1+[0:10], hist(nsig,[0:10])); set(h(1),'facecolor','k'); set(h(1),'edgecolor','k');
title(sprintf('No. of sig CA1 predictors'),'FontSize',24,'Fontweight','normal');
xlabel(['No. of sig CA1 cells in epoch'],'FontSize',24,'Fontweight','normal');
ylabel(['Hist of PFC cells'],'FontSize',24,'Fontweight','normal');
legend('All');
subplot(3,1,2)
h = bar(1+[0:10], hist(nsigpos,[0:10])); set(h(1),'facecolor','r'); set(h(1),'edgecolor','r');
legend('Pos');
subplot(3,1,3)
h = bar(1+[0:10], hist(nsigneg,[0:10])); set(h(1),'facecolor','c'); set(h(1),'edgecolor','c');
legend('Neg');


figure; hold on;%redimscreen_figforppt1;
subplot(3,1,1)
h = bar([0:0.1:1], hist(allfracsig,[0:0.1:1])); set(h(1),'facecolor','k'); set(h(1),'edgecolor','k');
title(sprintf('No. of sig CA1 predictors'),'FontSize',24,'Fontweight','normal');
xlabel(['Frac of sig CA1 cells in epoch'],'FontSize',24,'Fontweight','normal');
ylabel(['Hist of PFC cells'],'FontSize',24,'Fontweight','normal');
legend('All');
subplot(3,1,2)
h = bar([0:0.1:1], hist(allfracsigpos,[0:0.1:1])); set(h(1),'facecolor','r'); set(h(1),'edgecolor','r');
legend('Pos');
subplot(3,1,3)
h = bar([0:0.1:1], hist(allfracsigneg,[0:0.1:1])); set(h(1),'facecolor','c'); set(h(1),'edgecolor','c');
legend('Neg');



keyboard;



