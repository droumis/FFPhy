
% Will use getripples withing the analysis function to Nripples

clear; %close all;
runscript = 0;
savedata = 0; % save data option - only works if runscript is also on
figopt1 = 0; % Figure Options - Individual cells

savedir = '/data25/sjadhav/HPExpt/ProcessedData/';
%savefile = [savedir 'HP_riprate'];  clr = 'b'; % area = 'PFC';
savefile = [savedir 'HP_riprate_iri']; % With inte-ripple interval of 1 sec enforced
savefig1=0;


% Plot options
plotanimidx =  []; % To pick animals for plotting
plotdays = []; % If you only load data when runscript=0 and savedata=0, then this field will supplant days


% If runscript, run Datafilter and save data
if runscript == 1
    
    %Animal selection
    %-----------------------------------------------------
    animals = {'HPa','HPb'};
    
    %Filter creation
    %-----------------------------------------------------
    
    % Epoch filter
    % -------------
    dayfilter = '1:8'; % Shantanu - I am adding day filter to parse out epoch filter
    % Either Only do 1st w-track. 2 or 1 epochs per day
    % Or do Wtr1 and Wtr2, 2 epochs per day
    wtr1epochfilter = 'isequal($environment, ''wtr1'')';
    wtr2epochfilter = 'isequal($environment, ''wtr2'')';
    %runepochfilter = 'isequal($environment, ''wtr1'') || isequal($environment, ''wtr2'')';
    sleepepochfilter = 'isequal($type, ''sleep'')'; % Only pre and post sleep marked as sleep
    %postsleepepochfilter = 'isequal($environment, ''sleep'')'; % Only pre and post sleep marked as sleep
    
    % Time filter
    % -----------
    riptetfilter = '(isequal($descrip, ''riptet''))';
    
    % Iterator
    % --------
    %iterator = 'epochbehaveanal';
    iterator = 'multitetrodeanal'; % Can use this with riptetfilter, but will not use it in functions
    
    % Filter creation
    % ----------------
    wtr1f = createfilter('animal',animals,'days',dayfilter,'epochs',wtr1epochfilter,'eegtetrodes',riptetfilter,'iterator', iterator);
    wtr2f = createfilter('animal',animals,'days',dayfilter,'epochs',wtr2epochfilter,'eegtetrodes',riptetfilter,'iterator', iterator);
    sleepf = createfilter('animal',animals,'days',dayfilter,'epochs',sleepepochfilter,'eegtetrodes',riptetfilter,'iterator', iterator);
    
    disp('Done Filter Creation');
    
    % Set analysis function
    % ----------------------
    % doiri = inter ripple interval cane be limited to 1 sec
    wtr1f = setfilterfunction(wtr1f,'DFAsj_HPexpt_getriprate',{'ripples', 'tetinfo', 'pos'},'doiri',1); %
    wtr2f = setfilterfunction(wtr2f,'DFAsj_HPexpt_getriprate',{'ripples', 'tetinfo', 'pos'},'doiri',1); %
    sleepf = setfilterfunction(sleepf,'DFAsj_HPexpt_getriprate',{'ripples', 'tetinfo', 'pos'},'doiri',1); %
    % Going to call sj_getripples_tetinfo within function
    
    % Run analysis
    % ------------
    disp('wtr1');
    wtr1f = runfilter(wtr1f);
    disp('wtr2'); wtr2f = runfilter(wtr2f);
    disp('sleep'); sleepf = runfilter(sleepf);
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
%gatherdatafile = [savedir 'HP_riprate_gather'];
gatherdatafile = [savedir 'HP_riprate_iri_gather'];

if gatherdata
    % Get data
    % ----------
    % Wtr1
    % ------
    cnt=0; %cnt across animals
    allanimindex_wtr1=[]; allriprate_wtr1=[];
    for an = 1:length(wtr1f)
        for i=1:length(wtr1f(an).output{1})
            currdayep = [ wtr1f(an).output{1}(i).index(1,1),  wtr1f(an).output{1}(i).index(1,2)];
            anim_index{an}(i,:) = currdayep;
            % Only indexes
            animindex_wtr1=[an currdayep]; % Put animal index in front
            allanimindex_wtr1 = [allanimindex_wtr1; animindex_wtr1]; % Collect all Anim Day Epoch  Index
            % Data
            cnt=cnt+1;
            allriprate_wtr1(cnt) = wtr1f(an).output{1}(i).riprate;
        end
    end
    
    % Wtr2
    % ------
    cnt=0; %cnt across animals
    allanimindex_wtr2=[]; allriprate_wtr2=[];
    for an = 1:length(wtr2f)
        for i=1:length(wtr2f(an).output{1})
            currdayep = [ wtr2f(an).output{1}(i).index(1,1),  wtr2f(an).output{1}(i).index(1,2)];
            anim_index{an}(i,:) = currdayep;
            % Only indexes
            animindex_wtr2=[an currdayep]; % Put animal index in front
            allanimindex_wtr2 = [allanimindex_wtr2; animindex_wtr2]; % Collect all Anim Day Epoch  Index
            % Data
            cnt=cnt+1;
            allriprate_wtr2(cnt) = wtr2f(an).output{1}(i).riprate;
        end
    end
    
    % Sleep
    % ------
    cnt=0; %cnt across animals
    allanimindex_sleep=[]; allriprate_sleep=[];
    for an = 1:length(sleepf)
        for i=1:length(sleepf(an).output{1})
            currdayep = [ sleepf(an).output{1}(i).index(1,1),  sleepf(an).output{1}(i).index(1,2)];
            anim_index{an}(i,:) = currdayep;
            % Only indexes
            animindex_sleep=[an currdayep]; % Put animal index in front
            allanimindex_sleep = [allanimindex_sleep; animindex_sleep]; % Collect all Anim Day Epoch  Index
            % Data
            cnt=cnt+1;
            allriprate_sleep(cnt) = sleepf(an).output{1}(i).riprate;
        end
    end
    
    
    
    % Combine across animals, and for wtr1, also across epochs
    %---------------------------------------------
    
    % Wtr1
    %------
    daylist = unique(allanimindex_wtr1(:,[2]),'rows'); % Collapse across animals and epochs
    for d = 1:length(daylist)
        currday = daylist(d);
        currrate=[];
        curridxs = find(allanimindex_wtr1(:,2,:)==currday);
        for a = 1:length(curridxs)
            currrate = [currrate; allriprate_wtr1(curridxs(a))];
        end
        allriprate_wtr1c{currday} = currrate;
        allriprate_wtr1m(currday) = mean(currrate);
        allriprate_wtr1e(currday) = sem(currrate);
    end
    
    
    % Wtr2
    %------
    daylist = unique(allanimindex_wtr2(:,[2]),'rows'); % Collapse across animals and epochs
    for d = 1:size(daylist,1)
        currday = daylist(d);
        currrate=[];
        curridxs = find(allanimindex_wtr2(:,2,:)==currday);
        for a = 1:length(curridxs)
            currrate = [currrate; allriprate_wtr2(curridxs(a))];
        end
        allriprate_wtr2c{d} = currrate;
        allriprate_wtr2m(d,:) = mean(currrate);
        allriprate_wtr2e(d,:) = sem(currrate);
    end
    
    % Sleep - Keep pre and postsleep separate. Pre-sleep is always epoch1. Postsleep can be 5 or 7
    % ----------------------------------------
    dayepidxs = unique(allanimindex_sleep(:,[2 3]),'rows'); % Collapse across animals
    for ind = 1:size(dayepidxs,1)
        currrate=[];
        for a = 1:size(allanimindex_sleep,1)
            if dayepidxs(ind,:)==allanimindex_sleep(a,[2 3]) % Day-Epoch match
                currrate = [currrate; allriprate_sleep(a)];
            end
        end
        
        currday = dayepidxs(ind,1); currep = dayepidxs(ind,2);
        if currep==1
            allriprate_presleepc{currday} = currrate;
            allriprate_presleepm(currday,:) = mean(currrate);
            allriprate_presleepe(currday,:) = sem(currrate);
        else
            allriprate_postsleepc{currday} = currrate;
            allriprate_postsleepm(currday,:) = mean(currrate);
            allriprate_postsleepe(currday,:) = sem(currrate);
        end
    end
    
    % Save
    % -----
    if savegatherdata == 1
        save(gatherdatafile);
    end
    
else % gatherdata=0
    
    load(gatherdatafile);
    
end % end gather data


% ------------------
% Population Figures
% ------------------

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


if 1
    % 1) Ripple rate on Wtr1 and Wtr2
    % ----------------------------------------
    
    figure; hold on; redimscreen_figforppt1;
    set(gcf, 'Position',[205 136 723 446]);
    errorbar(allriprate_wtr1m,allriprate_wtr1e,'bo-','MarkerSize',6,'LineWidth',2);
    errorbar(allriprate_wtr2m,allriprate_wtr2e,'rs-','MarkerSize',6,'LineWidth',2);
    set(gca,'XLim',[0.5 8.5]);
    title(sprintf('SWR rate on W-tracks'));
    xlabel('Day','FontSize',xfont,'Fontweight','normal');
    ylabel('SWR rate (Hz)','FontSize',yfont,'Fontweight','normal');
    set(gca,'YLim',[0 0.9]);
    if savefig1==1,
        figfile = [figdir,'HP_RippleRate_Run'];
        print('-dpdf', figfile); print('-djpeg', figfile); saveas(gcf,figfile,'fig');
    end
    
end


if 1
    % 1) Ripple rate on Pre and Post Sleep
    % ----------------------------------------
    
    figure; hold on; redimscreen_figforppt1;
    set(gcf, 'Position',[205 136 723 446]);
    errorbar(allriprate_presleepm,allriprate_presleepe,'bo-','MarkerSize',6,'LineWidth',2);
    errorbar(allriprate_postsleepm,allriprate_postsleepe,'rs-','MarkerSize',6,'LineWidth',2);
    set(gca,'XLim',[0.5 8.5]);
    title(sprintf('SWR rate in pre (blue) and post sleep (red)'));
    xlabel('Day','FontSize',xfont,'Fontweight','normal');
    ylabel('SWR rate (Hz)','FontSize',yfont,'Fontweight','normal');
    set(gca,'YLim',[0 1]);
    if savefig1==1,
        figfile = [figdir,'HP_RippleRate_Sleep'];
        print('-dpdf', figfile); print('-djpeg', figfile); saveas(gcf,figfile,'fig');
    end
    
end


