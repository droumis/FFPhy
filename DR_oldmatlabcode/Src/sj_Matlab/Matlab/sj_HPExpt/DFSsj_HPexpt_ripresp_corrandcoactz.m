
% For ripple modulated and unmodulated cells, get corrcoeff and coactivez using the ripple response
% Instead of aligning to ripples again, use the saved output in ripplemod files


clear; %close all;
runscript = 0;
savedata = 0; % save data option - only works if runscript is also on
figopt1 = 0; % Figure Options - Individual cells

savedir = '/data25/sjadhav/HPExpt/ProcessedData/';
%savefile = [savedir 'HP_ripmod_corrandcoactz_PFC']; area = 'PFC'; clr = 'b'; % PFC
%savefile = [savedir 'HP_ripunmod_corrandcoactz_PFC']; area = 'PFC'; clr = 'b'; % PFC
% savefile = [savedir 'HP_ripmod_corrandcoactz_CA1']; area = 'CA1'; clr = 'r'; 
% savefile = [savedir 'HP_ripunmod_corrandcoactz_CA1']; area = 'CA1'; clr = 'r'; 
% savefile = [savedir 'HP_ripmod_corrandcoactz_CA1PFC']; area = 'CA1PFC'; clr = 'c'; 
% savefile = [savedir 'HP_ripunmod_corrandcoactz_CA1PFC']; area = 'CA1PFC'; clr = 'c'; 

% savefile = [savedir 'HP_ripmod_corrandcoactz_CA1allPFC']; area = 'CA1allPFC'; clr = 'c'; 
% savefile = [savedir 'HP_ripunmod_corrandcoactz_CA1allPFC']; area = 'CA1allPFC'; clr = 'c'; 
% savefile = [savedir 'HP_ripmodsleep_corrandcoactz_CA1allPFC']; area = 'CA1allPFC'; clr = 'c'; 
% savefile = [savedir 'HP_ripunmodsleep_corrandcoactz_CA1allPFC']; area = 'CA1allPFC'; clr = 'c'; 


% Sleep
% ------
% savefile = [savedir 'HP_ripmodsleep_corrandcoactz_PFC']; area = 'PFC'; clr = 'b'; % PFC
% savefile = [savedir 'HP_ripunmodsleep_corrandcoactz_PFC']; area = 'PFC'; clr = 'b'; % PFC
% savefile = [savedir 'HP_ripmodsleep_corrandcoactz_CA1']; area = 'CA1'; clr = 'r'; 
% savefile = [savedir 'HP_ripunmodsleep_corrandcoactz_CA1']; area = 'CA1'; clr = 'r'; 
% savefile = [savedir 'HP_ripmodsleep_corrandcoactz_CA1PFC']; area = 'CA1PFC'; clr = 'c'; 
% savefile = [savedir 'HP_ripunmodsleep_corrandcoactz_CA1PFC']; area = 'CA1PFC'; clr = 'c'; 


savefig1=0;


% Plot options
plotanimidx =  []; % To pick animals for plotting
plotdays = []; % If you only load data when runscript=0 and savedata=0, then this field will supplant days


% If runscript, run Datafilter and save data
if runscript == 1
    
    %Animal selection
    %-----------------------------------------------------
    animals = {'HPa','HPb', 'HPc'};
    
    %Filter creation
    %-----------------------------------------------------
    
    % Epoch filter
    % -------------
    dayfilter = '1:8'; % Shantanu - I am adding day filter to parse out epoch filter
    % Either Only do 1st w-track. 2 or 1 epochs per day
    % Or do Wtr1 and Wtr2, 2 epochs per day
     runepochfilter = 'isequal($environment, ''wtr1'') || isequal($environment, ''wtr2'')';
    % sleepepochfilter = 'isequal($type, ''sleep'')'; % Only pre and post sleep marked as sleep
    % sleepepochfilter = 'isequal($environment, ''postsleep'')'; % Only pre and post sleep marked as sleep
    
    % %Cell filter
    % %-----------
    
    % %PFC
    % %----
    %cellpairfilter = {'allcomb','strcmp($area, ''PFC'') && strcmp($ripmodtag, ''y'')','strcmp($area, ''PFC'') && strcmp($ripmodtag, ''y'')'}; % Ripple mod
    cellpairfilter = {'allcomb','strcmp($area, ''PFC'') && strcmp($ripmodtag, ''n'')','strcmp($area, ''PFC'') && strcmp($ripmodtag, ''n'')'}; % Ripple unmod
    
    % %CA1
    % %----
    %cellpairfilter = {'allcomb','(strcmp($area, ''CA1'') || strcmp($area, ''iCA1'')) && strcmp($ripmodtag, ''y'')','(strcmp($area, ''CA1'') || strcmp($area, ''iCA1'')) && strcmp($ripmodtag, ''y'')'};
    %cellpairfilter = {'allcomb','(strcmp($area, ''CA1'') || strcmp($area, ''iCA1'')) && strcmp($ripmodtag, ''n'')',('strcmp($area, ''CA1'') || strcmp($area, ''iCA1'')) && strcmp($ripmodtag, ''n'')'};
    
    % %CA1-PFC
    % %--------
    %cellpairfilter = {'allcomb','(strcmp($area, ''CA1'') || strcmp($area, ''iCA1'')) && strcmp($ripmodtag, ''y'')','strcmp($area, ''PFC'') && strcmp($ripmodtag, ''y'')'};
    %cellpairfilter = {'allcomb','(strcmp($area, ''CA1'') || strcmp($area, ''iCA1'')) && strcmp($ripmodtag, ''n'')','strcmp($area, ''PFC'') && strcmp($ripmodtag, ''n'')'};
    
    % %CA1all-PFCmodulated
    % %------------------------
    %cellpairfilter = {'allcomb','(strcmp($area, ''CA1'') || strcmp($area, ''iCA1'')) && ($numspikes >= 100)','strcmp($area, ''PFC'') && strcmp($ripmodtag, ''y'')'};
    %cellpairfilter = {'allcomb','(strcmp($area, ''CA1'') || strcmp($area, ''iCA1'')) && ($numspikes >= 100)','strcmp($area, ''PFC'') && strcmp($ripmodtag, ''n'')'};
    
    % Time filter - none.
    % -----------
    riptetfilter = '(isequal($descrip, ''riptet''))';
    
    % Iterator
    % --------
    iterator = 'singlecellanal';  % Have defined cellpairfilter. Can also use cellpair iterator with cell defn
    
    % Filter creation
    % ----------------
    modf = createfilter('animal',animals,'days',dayfilter,'epochs',runepochfilter, 'cellpairs',...
        cellpairfilter, 'iterator', iterator);
    
    disp('Done Filter Creation');
    
    % Set analysis function
    % ----------------------
    
    % % Need only the ripplemod structure for all of this.
    % % For calculating ripplealigned resp from scratch, you will need spikes, ripples, tetinfo, and pos
    modf = setfilterfunction(modf,'DFAsj_getripresp_corrandcoactz',{'ripplemod'}); %
    %modf = setfilterfunction(modf,'DFAsj_getripresp_corrandcoactz',{'ripplemod','cellinfo','spikes', 'ripples', 'tetinfo', 'pos'}); %
    
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
gatherdata = 0; savegatherdata = 0;
 gatherdatafile = [savedir 'HP_ripmod_corrandcoactz_PFC_gather']; area = 'PFC'; kind = 'mod'; state = '';
% gatherdatafile = [savedir 'HP_ripunmod_corrandcoactz_PFC_gather']; area = 'PFC'; kind = 'unmod'; state = '';
% gatherdatafile = [savedir 'HP_ripmod_corrandcoactz_CA1_gather']; area = 'CA1'; kind = 'mod'; state = '';
% gatherdatafile = [savedir 'HP_ripunmod_corrandcoactz_CA1_gather']; area = 'CA1'; kind = 'unmod'; state = '';
% gatherdatafile = [savedir 'HP_ripmod_corrandcoactz_CA1PFC_gather']; area = 'CA1PFC'; kind = 'mod'; state = '';
% gatherdatafile = [savedir 'HP_ripunmod_corrandcoactz_CA1PFC_gather']; area = 'CA1PFC'; kind = 'unmod'; state = '';

%gatherdatafile = [savedir 'HP_ripmod_corrandcoactz_CA1allPFC_gather']; area = 'CA1allPFC'; kind = 'mod'; state = '';
%gatherdatafile = [savedir 'HP_ripunmod_corrandcoactz_CA1allPFC_gather']; area = 'CA1allPFC'; kind = 'unmod'; state = '';
%gatherdatafile = [savedir 'HP_ripmodsleep_corrandcoactz_CA1allPFC_gather']; area = 'CA1allPFC'; kind = 'mod'; state = 'sleep';
%gatherdatafile = [savedir 'HP_ripunmodsleep_corrandcoactz_CA1allPFC_gather']; area = 'CA1allPFC'; kind = 'unmod'; state = 'sleep';


% gatherdatafile = [savedir 'HP_ripmodsleep_corrandcoactz_PFC_gather']; area = 'PFC'; kind = 'mod'; state = 'sleep';
% gatherdatafile = [savedir 'HP_ripunmodsleep_corrandcoactz_PFC_gather']; area = 'PFC'; kind = 'unmod'; state = 'sleep';
% gatherdatafile = [savedir 'HP_ripmodsleep_corrandcoactz_CA1_gather']; area = 'CA1'; kind = 'mod'; state = 'sleep';
% gatherdatafile = [savedir 'HP_ripunmodsleep_corrandcoactz_CA1_gather']; area = 'CA1'; kind = 'unmod'; state = 'sleep';
% gatherdatafile = [savedir 'HP_ripmodsleep_corrandcoactz_CA1PFC_gather']; area = 'CA1PFC'; kind = 'mod'; state = 'sleep';
% gatherdatafile = [savedir 'HP_ripunmodsleep_corrandcoactz_CA1PFC_gather']; area = 'CA1PFC'; kind = 'unmod'; state = 'sleep';


if gatherdata
    
    % Parameters if any
    % -----------------
    
    % -------------------------------------------------------------
    
    cnt=0; 
    allanimindex=[]; alldataraster=[]; alldatahist = []; all_Nspk=[];
    
    for an = 1:length(modf)
        for i=1:length(modf(an).output{1})
                cnt=cnt+1;
                anim_index{an}(cnt,:) = modf(an).output{1}(i).index;
                % Only indexes
                animindex=[an modf(an).output{1}(i).index]; % Put animal index in front
                allanimindex = [allanimindex; animindex]; % Collect all Anim Day Epoch Tet Cell Index
                % Data
                allr(cnt) = modf(an).output{1}(i).r; 
                allp_shuf(cnt) = modf(an).output{1}(i).p_shuf;
                allp(cnt) = modf(an).output{1}(i).p;
                allr_rdm(cnt) = modf(an).output{1}(i).r_rdm; 
                allp_rdmshuf(cnt) = modf(an).output{1}(i).p_rdmshuf;
                allp_rdm(cnt) = modf(an).output{1}(i).p_rdm;
                
                allcoactivez(cnt) = modf(an).output{1}(i).coactivez; 
                allcoactivez_pshuf(cnt) = modf(an).output{1}(i).coactivez_pshuf;
                allcoactivez_rdm(cnt) = modf(an).output{1}(i).coactivez_rdm;
                allcoactivez_prdmshuf(cnt) = modf(an).output{1}(i).coactivez_prdmshuf;
        
        end
        
    end
    
    % DONT MAKE STRUCTURE OR CONSOLIDATE ACROSS EPOCHS FOR NOW
    % ---------------------------------------------------------
    % To consolidate across epochs, I have to just get trialResps from the functions.
    % Then here, I have to combine trialResps across epochs, and then run the corrcoeff and coactivez function
    
    % For post-sleep, I will have only 1 epoch.
    
    
    
    % Save
    % -----
    if savegatherdata == 1
        save(gatherdatafile);
    end
    
else % gatherdata=0
    
    load(gatherdatafile);
    
end % end gather data

corrRateRipples = nanmean(allp_shuf < 0.05)
corrRateRdm = nanmean(allp_rdmshuf < 0.05)

[r_realrdm p_realrdm] = ttest(allp_shuf<0.05,allp_rdmshuf<0.05)


% ------------------
% Population Figures
% ------------------

forppr = 0;
% If yes, everything set to redimscreen_figforppr1
% If not, everything set to redimscreen_figforppt1

figdir = '/data25/sjadhav/HPExpt/Figures/RippleMod/Popln/';
summdir = figdir;
set(0,'defaultaxesfontweight','normal'); set(0,'defaultaxeslinewidth',2);

if forppr==1
    set(0,'defaultaxesfontsize',16);
    tfont = 18; % title font
    xfont = 16;
    yfont = 16;
else
    set(0,'defaultaxesfontsize',40);
    tfont = 40;
    xfont = 40;
    yfont = 40;
end

if strcmp(state,'sleep'),
    statename = 'Sleep';
else
    statename = 'Run';
end


if 1
    % 1) Histogram for Corr - Real data and Rdm 
    % -----------------------------------------------------
    
    figure; hold on; redimscreen_figforppt1;
    %set(gcf, 'Position',[205 136 723 446]);
    %xaxis = min(allr):0.1:max(allr);
    xaxis = -1:0.04:1;
    h = histc(allr,xaxis); normh = h./max(h);
    plot(xaxis,normh,clr,'Linewidth',3);
    h = histc(allr_rdm,xaxis); normh = h./max(h);
    plot(xaxis,normh,[clr '--'],'Linewidth',3);
    %h = histc(allr_shuf,xaxis); normh = h./max(h);
    %plot(xaxis,normh,[clr ':'],'Linewidth',3);
    legend('Ripple ev','Random ev');

    title(sprintf('%s %s - Rip%s units: Corr Coeff Hist', area, statename, kind),'FontSize',tfont,'Fontweight','normal')
    xlabel('Corr Coeff','FontSize',xfont,'Fontweight','normal');
    ylabel(sprintf('Fraction of pairs (Total %d)', length(allr)),'FontSize',yfont,'Fontweight','normal');

    if strcmp(area,'PFC')
        set(gca,'XLim',[-0.4 0.5]);
    else
        set(gca,'XLim',[-0.2 0.25]);
    end
    text(0.07,0.7,sprintf('Corr Rate Rip: %0.2f',corrRateRipples),'FontSize',30,'Fontweight','normal');
    text(0.07,0.6,sprintf('Corr Rate Rdm: %0.2f',corrRateRdm),'FontSize',30,'Fontweight','normal');
    text(0.07,0.5,sprintf('Diff Sig: %0.3f',p_realrdm),'FontSize',30,'Fontweight','normal');
        
    figfile = [figdir,area,'_',statename,'_Ripple',kind,'_CorrCoeffHist']
    if savefig1==1,
        print('-depsc2', figfile); print('-djpeg', figfile); saveas(gcf,figfile,'fig');
    end
    
end


if 1
     % 1) Histogram for Coactivez - Real data and Rdm 
    % ------------------------------------------------------------------
    
    figure; hold on; redimscreen_figforppt1;
    %set(gcf, 'Position',[205 136 723 446]);
    xaxis = min(allcoactivez):0.4:max(allcoactivez);
    h = histc(allcoactivez,xaxis); normh = h./max(h);
    plot(xaxis,normh,clr,'Linewidth',3);
    h = histc(allcoactivez_rdm,xaxis); normh = h./max(h);
    plot(xaxis,normh,[clr '--'],'Linewidth',3);
    %h = histc(allcoactivez_shuf,xaxis); normh = h./max(h);
    %plot(xaxis,normh,[clr ':'],'Linewidth',2);
    legend('Ripple ev','Random ev');
    
   
    title(sprintf('%s %s - Rip%s units: CoactiveZ Hist', area, statename, kind),'FontSize',tfont,'Fontweight','normal')
    xlabel('Coactive Z','FontSize',xfont,'Fontweight','normal');
    ylabel('Fraction of cells','FontSize',yfont,'Fontweight','normal');

    set(gca,'XLim',[-4 5]);
    %[h,p] = ttest(allcoactivez,allcoactivez_rdm); % need a nonparametric test - like shuffle
    %text(2,0.5,sprintf('Diff Sig: %0.3f',pval),'FontSize',30,'Fontweight','normal');
        
    
    figfile = [figdir,area,'_',statename,'_Ripple',kind,'_CoactivezHist']
    if savefig1==1,
        print('-depsc2', figfile); print('-djpeg', figfile); saveas(gcf,figfile,'fig');
    end
    
    
end



% ------------------------------------------------------------------
% COMBINING PLOTS ACROSS FILES
% ------------------------------------------------------------------


keyboard;

savedir = '/data25/sjadhav/HPExpt/ProcessedData/';
figdir = '/data25/sjadhav/HPExpt/Figures/RippleMod/Popln/';

% Define area
area = 'CA1allPFC'; state ='sleep'; % state = '';, or state = 'sleep';

if strcmp(state,'sleep'),
    statename = 'Sleep';
else
    statename = 'Run';
end


% Modulated Units
load([savedir 'HP_ripmod',state,'_corrandcoactz_',area,'_gather'])
% Corr Fig
figure(1); hold on; redimscreen_figforppt1;
set(gcf, 'Position',[205 136 723 446]);
%xaxis = min(allr):0.1:max(allr);
xaxis = -1:0.05:1;
h = histc(allr,xaxis); normh = h./max(h);
plot(xaxis,normh,clr,'Linewidth',3);

% % Coactive Z Fig
% figure(2); hold on; redimscreen_figforppt1;
% set(gcf, 'Position',[205 136 723 446]);
% xaxis = min(allcoactivez):0.5:max(allcoactivez);
% h = histc(allcoactivez,xaxis); normh = h./max(h);
% plot(xaxis,normh,clr,'Linewidth',3);

allrm = allr; allpm = allp_shuf; allcozm = allcoactivez;


% UnModulated Units
load([savedir 'HP_ripunmod',state,'_corrandcoactz_',area,'_gather'])
% Corr Fig
figure(1); hold on;
xaxis = -1:0.05:1;
h = histc(allr,xaxis); normh = h./max(h);
plot(xaxis,normh,[clr '--'],'Linewidth',3);

% % Coactive Z Fig
% figure(2); hold on; redimscreen_figforppt1;
% set(gcf, 'Position',[205 136 723 446]);
% xaxis = min(allcoactivez):0.5:max(allcoactivez);
% h = histc(allcoactivez,xaxis); normh = h./max(h);
% plot(xaxis,normh,[clr '--'],'Linewidth',3);

legend('Rip Mod','Rip Unmod');
title(sprintf('%s %s - units: Corr Coeff Hist', area, statename),'FontSize',tfont,'Fontweight','normal')
xlabel('Corr Coeff','FontSize',xfont,'Fontweight','normal');
ylabel('Fraction of cells','FontSize',yfont,'Fontweight','normal');


corrRateRipMod = nanmean(allpm_shuf < 0.05), corrRateRipUnMod = nanmean(allp_shuf < 0.05),
[r_modunmod p_modunmod] = ttest2(allpm_shuf<0.05,allp_shuf<0.05)

set(gca,'XLim',[-0.2 0.25]);
text(0.07,0.7,sprintf('Corr Rate Mod: %0.2f',corrRateRipMod),'FontSize',30,'Fontweight','normal');
text(0.07,0.6,sprintf('Corr Rate Unmod: %0.2f',corrRateRipUnMod),'FontSize',30,'Fontweight','normal');
text(0.07,0.5,sprintf('Diff Sig: %0.3f',p_modunmod),'FontSize',30,'Fontweight','normal');

figfile = [figdir,area,'_',statename,'_RippleModvsUnmod_CorrCoeffHist']
if savefig1==1,
    print('-depsc2', figfile); print('-djpeg', figfile); saveas(gcf,figfile,'fig');
end











