
% From placefield1 and placefield_params
% Get iso qual for placefields for exp and con groups

clear; close all;
runscript = 0;
savedata = 0; % save data option - only works if runscript is also on
figopt1 = 0; % Figure Options - Plot individual cell place plots

savedir = '/data25/sjadhav/RippleInterruption/ProcessedData/';
savefile = [savedir 'PlaceFieldIso'];
%savefile = [savedir 'REGrp_PlaceFieldParamsn'];  % n is after std is changed for mapfields
%savefile = [savedir 'RCGrp_PlaceFieldParamsn'];

% Plot options
plotanimidx = []; % To pick animals for plotting
plotdays = []; % If you only load data when runscript=0 and savedata=0, then this field will supplant days


% If runscript, run Datafilter and save data
if runscript == 1
    
    
    %Animal selection
    %-----------------------------------------------------
    Expanimals = {'REc,''REd','REe','REf'};
    Conanimals = {'RCa','RCb','RCc','RCd'};
    %Expanimals = {'REd','REe','REf'};
    %Conanimals = {'RCb','RCc','RE1'};
    
    %Filter creation
    %-----------------------------------------------------
    
    % Epoch filter
    % -------------
    dayfilter = '1:8'; % Shantanu - I am adding day filter to parse out epoch filter
    epochfilter = 'isequal($type, ''run'') || isequal($type, ''sleep'')';
    %epochfiltersl = 'isequal($type, ''sleep'')'; 
    
    % Cell filter
    % -----------
    placecellfilter = 'strcmp($tag, ''PyrSR'')';
    
    % Iterator
    % --------
    iterator = 'singlecellanal';
    
    % Filter creation
    % ----------------
    Expplaf = createfilter('animal',Expanimals,'days',dayfilter,'epochs',epochfilter,'cells',placecellfilter,'iterator', iterator);
    Conplaf = createfilter('animal',Conanimals,'days',dayfilter,'epochs',epochfilter,'cells',placecellfilter,'iterator', iterator);
    
    % Set analysis function
    % ----------------------
    % Get iso qual from saved files
    Expplaf = setfilterfunction(Expplaf, 'DFAsj_getclustqual', {'clustqual'});
    Conplaf = setfilterfunction(Conplaf, 'DFAsj_getclustqual', {'clustqual'});
    
    % Run analysis
    % ------------
    Expplaf = runfilter(Expplaf);  
    Conplaf = runfilter(Conplaf);  
    
    disp('Finished running filter script');
    %--------------------- Finished Filter Function Run -------------------
    
    if savedata == 1
        clear figopt1 runscript plotdays plotanimidx savedata 
        save(savefile);
    end
    
else
    
    load(savefile);
    
end  % end runscript

if ~exist('savedata')
    return
end

% ------------------------------
% Figure and Font Sizes

forppr = 1;
% If yes, everything set to redimscreen_figforppr1
% If not, everything set to redimscreen_figforppt1

figdir = '/data25/sjadhav/RippleInterruption/Figures/01AugSep11_RippleDisFigs/PlaceFields/Params/';
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

clr = {'b','r','g','c','m','y','k','r'};

% Saving processed data and figures
savefig1 = 1;
%saveprocdatafile = [savedir,'PlaceFieldParamData'];
%saveprocdata = 0; saveproctag = 'exp';

% ---------------------------------------

str=['Con';'Exp'];
for g = 1:size(str,1)   % Do Exp and Con groups separately
    
  % Anim
    if ~isempty(plotanimidx)
        useanim = plotanimidx;
    else
        useanim = eval(['1:length(',str(g,:),'plaf);']); 
    end 
    
    % Days
    if ~isempty(plotdays)
        days = plotdays;
    else
        days = eval(['unique(',str(g,:),'plaf(1).epochs{1}(:,1));']);   % all days get from struct
    end
 
    allepochs = eval(['unique(',str(g,:),'plaf(1).epochs{1}(:,2));']);  % SHould be all epochs 1 to 5
    totanim = length(useanim);
 
    % saving across animals and days. Can keep track of animals and days as well if I want to
    cntallcells=0; 
    for anidx=1:totanim % Across animals
        an=useanim(anidx);
        
        % For run
        currindex=[]; 
        eval(['nidxs_an = length(',str(g,:),'plaf(an).output{1});']);
        for i=1:nidxs_an
            eval(['currindex(i,:) = ',str(g,:),'plaf(an).output{1}(i).index;']);
        end
        % should have same number for all epochs
        ep1s = find(currindex(:,2)==1); 
        ep2s = find(currindex(:,2)==2); 
        ep3s = find(currindex(:,2)==3); 
        ep4s = find(currindex(:,2)==4); 
        ep5s = find(currindex(:,2)==5); 
        
        for n=1:length(ep2s)
            cntallcells = cntallcells+1;
            
            eval([str(g,:),'isoldist(cntallcells,1) =',str(g,:),'plaf(an).output{1}(ep1s(n)).isoldist;']);
            eval([str(g,:),'isoldist(cntallcells,2) =',str(g,:),'plaf(an).output{1}(ep2s(n)).isoldist;']);
            eval([str(g,:),'isoldist(cntallcells,3) =',str(g,:),'plaf(an).output{1}(ep3s(n)).isoldist;']);
            eval([str(g,:),'isoldist(cntallcells,4) =',str(g,:),'plaf(an).output{1}(ep4s(n)).isoldist;']);
            eval([str(g,:),'isoldist(cntallcells,5) =',str(g,:),'plaf(an).output{1}(ep5s(n)).isoldist;']);
            
        end % end ep2s
    end % anim
end % grp

x=(nanmax(Expisoldist,[],2));
y=(nanmax(Conisoldist,[],2));
x=(nanmean(Expisoldist,2));
y=(nanmean(Conisoldist,2));

[h,p,c,st]=ttest2(x,y);


keyboard;

%****************************************
% Figures -
%*******************************************

% Mean across epochs: Mean rates / Peak rates, etc

% Expmeanratesl = mean(Expmeanratesl,2); % for sleep
% Expmeanrate = mean(Expmeanrate,2); % Hz - for runs only
% Exppeakrate = mean(Exppeakrate,2); % Hz - for runs only
% Exptotal_lthunder1 = mean(Exptotal_lthunder1,2)*binsize; % in cm
% Exptotal_lthunder3 = mean(Exptotal_lthunder3,2)*binsize;
% Exptotal_trajlth = mean(Exptotal_trajlth,2)*binsize; % Total Traj Lth - length of all 4 trajs. 
% Exptotal_fracunder1 = mean(Exptotal_fracunder1,2);
% Exptotal_fracunder3 = mean(Exptotal_fracunder3,2);
% 
% Conmeanratesl = mean(Conmeanratesl,2); % for sleep
% Conmeanrate = mean(Conmeanrate,2);
% Conpeakrate = mean(Conpeakrate,2);
% Contotal_lthunder1 = mean(Contotal_lthunder1,2)*binsize;
% Contotal_lthunder3 = mean(Contotal_lthunder3,2)*binsize;
% Contotal_trajlth = mean(Contotal_trajlth,2)*binsize; % Total Traj Lth - length of all 4 trajs. 
% Contotal_fracunder1 = mean(Contotal_fracunder1,2);
% Contotal_fracunder3 = mean(Contotal_fracunder3,2);
% 


if figopt1 == 1
 
            
    % Proportion Track Active under 1
    % --------------------------------
%     figure; hold on;
%     if forppr==1
%         redimscreen_figforppr1;
%     else
%         redimscreen_figforppt1;
%     end
%     bar(1,mean(Contotal_fracunder1),'b');
%     bar(2,mean(Exptotal_fracunder1),'r');
%     errorbar(1,mean(Contotal_fracunder1),sem(Contotal_fracunder1),'k','LineWidth',3);
%     errorbar(2,mean(Exptotal_fracunder1),sem(Exptotal_fracunder1),'k','LineWidth',3);
%     [h_fracunder1,p_fracunder1] = ttest2(Contotal_fracunder1,Exptotal_fracunder1); %h=0, p=0.76
%     title('Proportion track active under 1Hz (all 4 traj)');
%     ylabel('Proportion track active');
%     set(gca,'XTick',[1 2],'XTickLabel',{'Con';'Exp'});
%     if savefig1==1,
%         figfile = [figdir,'PlaceFieldFrac_1Hz'];
%         print('-dpdf', figfile); print('-djpeg', figfile); saveas(gcf,figfile,'fig');
%     end
 
end % end figopt



keyboard;



















