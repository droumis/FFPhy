% Ver4 : Starting 10Feb2014 - Sync codes with everyone

% Glm model for PFC rip-align response as a function of rip-aligned CA1 popln responses
% Goal is to get coefficients and significance
% Get ripple aligned responses directly from ripplemod structure
% Also get pairwise corrlns for now, without shuffle sign for speed. You cal always compare to saved apirwise corrln files later 


clear; %close all;
runscript = 0;
savedata = 0; % save data option - only works if runscript is also on
figopt1 = 0; % Figure Options - Individual cells

savedir = '/data25/sjadhav/HPExpt/ProcessedData/';
[y, m, d] = datevec(date);

% Version4
% --------
if runscript==1
    val=1; savefile = [savedir 'HP_ripmod_glmfit_',num2str(m),'-',num2str(d),'-',num2str(y)]; area = 'PFC'; clr = 'b'; % PFC
else
    val=1; savefile = [savedir 'HP_ripmod_glmfit_2-12-2014'];
end

% Pre version4
%val=1; savefile = [savedir 'HP_ripmod_glmfit2_feb14']; area = 'PFC'; clr = 'b'; % PFC
%val=1; savefile = [savedir 'HP_ripmod_glmfit2']; area = 'PFC'; clr = 'b'; % PFC


savefig1=0;


% Plot options
plotanimidx =  []; % To pick animals for plotting
plotdays = []; % If you only load data when runscript=0 and savedata=0, then this field will supplant days


% If runscript, run Datafilter and save data
if runscript == 1
    
    %Animal selection
    %-----------------------------------------------------
    animals = {'HPa','HPb','HPc','Ndl'};
    
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
    
    % Iterator
    % --------
    iterator = 'multicellanal';
    
    % Filter creation
    % ----------------
    modf = createfilter('animal',animals,'epochs',runepochfilter, 'cells',...
        cellfilter, 'iterator', iterator);
    
    disp('Done Filter Creation');
    
    % Set analysis function
    % ----------------------
    
    % % Need only the ripplemod structure for all of this. cellinfo needed to parse cell area 
    % % For calculating ripplealigned resp from scratch, you will need spikes, ripples, tetinfo, and pos
    % % If it is across regions, you will need cellinfo to get area where cells are recorded from
    %modf = setfilterfunction(modf,'DFAsj_glm_ripalign',{'ripplemod','cellinfo'}); %
    modf = setfilterfunction(modf,'DFAsj_glm_ripalign4',{'ripplemod','cellinfo'}); %

    
    % Run analysis
    % ------------
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

% Version 4 onwards
% Filename of gatherdata. If generating, put current date. If not, then load a file from previous date.
% -----------------------------------------------------------------------------------------------------
switch val
    case 1
        if gatherdata
            gatherdatafile = [savedir 'HP_ripmod_glmfit_gather_',num2str(m),'-',num2str(d),'-',num2str(y)],
        else
            gatherdatafile = [savedir 'HP_ripmod_glmfit_gather_2-12-2014'],
        end
        % Pre version4
        %gatherdatafile = [savedir 'HP_ripmod_glmfit_gather2'];
end


if gatherdata    
    % Parameters if any
    % -----------------
    % -------------------------------------------------------------
    
    cnt=0; % How many kept? Any condition?
    %Glm
    allglmidxs = []; allmodelb=[]; allmodelp=[]; allmodelfits=[]; 
    allnsig=[]; allnsigpos=[]; allnsigneg=[];  allfracsig=[]; allfracsigpos=[]; allfracsigneg=[];
    %Corr
    allrcorr=[];
    allpcorr=[];
    allnsimul=[];
    
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
                allglmidxs = [allglmidxs; modf(an).output{1}(i).glmidxs]; % Need to put "an" in there
                allcorridxs = [allglmidxs; modf(an).output{1}(i).corridxs];
                % Glm data
                allmodelb = [allmodelb; modf(an).output{1}(i).allmodelb];
                allmodelp = [allmodelp; modf(an).output{1}(i).allmodelp];
                allmodelfits{cnt} = modf(an).output{1}(i).allmodelfits;
                modf(an).output{1}(i).nsig;
                allnsig = [allnsig; modf(an).output{1}(i).nsig'];
                allnsigpos = [allnsigpos; modf(an).output{1}(i).nsigpos'];
                allnsigneg = [allnsigneg; modf(an).output{1}(i).nsigneg'];
                allfracsig = [allfracsig; modf(an).output{1}(i).fracsig'];
                allfracsigpos = [allfracsigpos; modf(an).output{1}(i).fracsigpos'];
                allfracsigneg = [allfracsigneg; modf(an).output{1}(i).fracsigneg'];
                % Corr data
                allrcorr = [allrcorr; modf(an).output{1}(i).rcorr];
                allpcorr = [allpcorr; modf(an).output{1}(i).pcorr];
                allnsimul = [allnsimul; modf(an).output{1}(i).nsimul];     
            end
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
h = bar(1+[0:10], hist(allnsig,[0:10])); set(h(1),'facecolor','k'); set(h(1),'edgecolor','k');
title(sprintf('No. of sig CA1 predictors'),'FontSize',24,'Fontweight','normal');
xlabel(['No. of sig CA1 cells in epoch'],'FontSize',24,'Fontweight','normal');
ylabel(['Hist of PFC cells'],'FontSize',24,'Fontweight','normal');
legend('All');
subplot(3,1,2)
h = bar(1+[0:10], hist(allnsigpos,[0:10])); set(h(1),'facecolor','r'); set(h(1),'edgecolor','r');
legend('Pos');
subplot(3,1,3)
h = bar(1+[0:10], hist(allnsigneg,[0:10])); set(h(1),'facecolor','c'); set(h(1),'edgecolor','c');
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



