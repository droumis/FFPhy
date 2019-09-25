
% Plotting of both groups together
% DIO is controlled for

clear; close all;
runscript = 1;
savedata = 1; % save data option - only works if runscript is also on
figopt1=0; % Figure Options for Ripple Rate

figopt2=0; % Figure Options for Ripple Size
% Sub-options under figopt2
figoptsize = 1; % For plotting ripple-size plots
figoptntet = 0;

ntet=1;

savedir = '/data25/sjadhav/RippleInterruption/ProcessedData/';

% First set - Use DFAsj_getriprate_noDIO
% savefile1 = [savedir 'RippleRate_ExpGrp_nostim'];
% savefile2 = [savedir 'RippleRate_ConGrp_nostim'];
% savefile = [savedir 'RippleRate_All_nostim'];

% Second Set - Use DFTF_getstimtimes to filter out stim times, and then use
savefile1 = [savedir 'RippleRate_ExpGrp_filtstim'];
savefile2 = [savedir 'RippleRate_ConGrp_filtstim'];
savefile = [savedir 'RippleRate_All_filtstim'];

% If runscript, run Datafilter and save data
if runscript == 1
    
    %Animal selection
    %-----------------------------------------------------
    Expanimals = {'REc','REd','REe','REf'};
    Conanimals = {'RCa','RCb','RCc','RCd'};
    
    %Filter creation
    %--------------------------------------------------------
    
    % epoch filter
    dayfilter = '1:8'; % Shantanu - I am adding day filter to parse out epoch filter
    epochfilter = 'isequal($type, ''sleep'')';
    
    % cell filter None
    % placecellfilter = '(strcmp($tag, ''PyrSR'') || strcmp($tag, ''PyrS''))';
    
    % Tet filter for ripple detection
    riptetfilter = '(isequal($descrip, ''riptet''))';
    
    % time filter - Speed <=2 for ripples in sleep box
    
    % For Set 1
    %timefilter = {{'DFTFsj_getvelpos', '(($absvel <= 2))'}};
    
    % For Set 2 - Filter out stimtimes as well
    timefilter = {{'DFTFsj_getvelpos', '(($absvel <= 2))'},...
        {'DFTFsj_getstimtimes','($nstim == 0)','tetfilter',riptetfilter}}; % Default timewin is 0.1=100ms
    
    % iterator
    iterator = 'multitetrodeanal'; % / iterator = eeganal;
    
    % filter creation
    Expripf = createfilter('animal',Expanimals,'days',dayfilter,'epochs',epochfilter,'eegtetrodes',riptetfilter,'excludetime', timefilter,'iterator', iterator);
    Conripf = createfilter('animal',Conanimals,'days',dayfilter,'epochs',epochfilter,'eegtetrodes',riptetfilter,'excludetime', timefilter,'iterator', iterator);
    
    % set analysis function
    
    % For Set 1
    %     Expripf = setfilterfunction(Expripf, 'DFAsj_getriprate_noDIO', {'ripples','DIO'}, 'numtetrodes', ntet,'minthresh',3);
    %     Conripf = setfilterfunction(Conripf, 'DFAsj_getriprate_noDIO', {'ripples','DIO'}, 'numtetrodes', ntet,'minthresh',3);
    
    % For Set 2
    Expripf = setfilterfunction(Expripf, 'DFAsj_getriprate', {'ripplesep1'}, 'numtetrodes', ntet,'minthresh',3);
    Conripf = setfilterfunction(Conripf, 'DFAsj_getriprate', {'ripplesep1'}, 'numtetrodes', ntet,'minthresh',3);
    
    % run analysis
    Expripf = runfilter(Expripf);  % Ripple rate
    Conripf = runfilter(Conripf);
    
    %--------------------- Finished Filter Function Run -------------------
    
    disp('Finished running filter');
    
    if savedata == 1
        clear figopt1 figopt2 figoptsize figoptntet runscript savedata
        save(savefile1,'Expripf','ntet');
        save(savefile2,'Conripf','ntet');
        save(savefile);
    end
    
else
    
    load(savefile);
    
end  % end runscript

if ~exist('savedata')
    return
end




%--------------------- End Run Script -------------------

% ------------------------------
% Figure and Font Sizes

forppr = 1;
% If yes, everything set to redimscreen_figforppr1
% If not, everything set to redimscreen_figforppt1

figdir = '/data25/sjadhav/RippleInterruption/Figures/01AugSep11_RippleDisFigs/Behavior/';
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

% --------Parameters-------------------------------
novel=1:3; nov=novel;
fam=4:8;


% ---------Extract Data and Get ripple rate---------

str=['Exp';'Con'];
days = unique(Expripf(1).epochs{1}(:,1));
allepochs = unique(Expripf(1).epochs{1}(:,2));


for g = 1:size(str,1)   % Do Exp and Con groups separately
    
    % Across animals: Inititalize to store ripplerate_ep, ripplepercenttime_ep, and totalstilltime_ep;
    for ep = 1:length(allepochs)
        currep = allepochs(ep);
        eval([str(g,:),'riprate_ep',num2str(currep),'=[];']);
        eval([str(g,:),'ripper_ep',num2str(currep),'=[];']);
        eval([str(g,:),'stilltime_ep',num2str(currep),'=[];']);
    end
    
    
    eval(['totanim = length(',str(g,:),'ripf);']);
    
    for an=1:totanim % Across animals
        
        % For each animal - initialize - I am not using this later
        % Not Being used
        eval([str(g,:),'riprate = [];']);
        eval([str(g,:),'ripper = [];'])
        eval([str(g,:),'stilltime = [];'])
        % Not being used
        
        for d=1:length(days)
            currday = days(d);
            for ep = 1:length(allepochs)
                currep = allepochs(ep);
                
                index = eval(['find( (',str(g,:),'ripf(an).epochs{1}(:,1)==currday) & (',str(g,:),'ripf(an).epochs{1}(:,2)==currep) );']);
                eval([str(g,:),'riprate(d,ep) =',str(g,:),'ripf(an).output{1}(index).rip(1);']);
                eval([str(g,:),'ripper(d,ep) =',str(g,:),'ripf(an).output{1}(index).rip(2);']);
                eval([str(g,:),'stilltime(d,ep) =',str(g,:),'ripf(an).output{1}(index).rip(3);']);
                
                eval([str(g,:),'riprate_ep',num2str(currep),'(an,d) =',str(g,:),'ripf(an).output{1}(index).rip(1);']);
                eval([str(g,:),'ripper_ep',num2str(currep),'(an,d) =',str(g,:),'ripf(an).output{1}(index).rip(2);']);
                eval([str(g,:),'stilltime_ep',num2str(currep),'(an,d) =',str(g,:),'ripf(an).output{1}(index).rip(3);']);
                
                % Ripple  Size
                eval([str(g,:),'ripntet_ep',num2str(currep),'{an}{d} =',str(g,:),'ripf(an).output{1}(index).ripntet;']);
                eval([str(g,:),'ripsize_ep',num2str(currep),'{an}{d} =',str(g,:),'ripf(an).output{1}(index).ripsize;']);
                
            end
        end
        
        % Not being used
        eval([str(g,:),'allriprate{an}=',str(g,:),'riprate;']);
        eval([str(g,:),'allripper{an}=',str(g,:),'ripper;']);
        eval([str(g,:),'allstilltime{an}=',str(g,:),'stilltime;']);
        % Not being used
        
    end
    
    %eval(['clear ',str(g,:),'ripf;']);
  
    % Calculate
    
    %  REc day 8 - for ntet=2? Noise?- (Implement NaNs)
    %     if ntet==2
    %         Expriprate_ep1(1,8)=NaN; Expriprate_ep1(1,4)=NaN;
    %         Expriprate_ep3(1,8)=NaN; Expriprate_ep1(1,4)=NaN;
    %         Expriprate_ep5(1,8)=NaN; Expriprate_ep1(1,4)=NaN;
    %     else
    %         Expriprate_ep1(1,8)=NaN;
    %         Expriprate_ep3(1,8)=NaN;
    %         Expriprate_ep5(1,8)=NaN;
    %     end
    
    % Sleep epochs - Get means of variables, and Fraction Increase/Dec in variables (ripple rate) as fraction of sleep 1
    sleeps = [1,3,5];
    for ep = sleeps
        
        % For Ripple Rate
        % --------------
        
        eval([str(g,:),'frriprate_ep',num2str(ep),' =',str(g,:),'riprate_ep',num2str(ep),'./',str(g,:),'riprate_ep1;']); % Fraction change: Matrix divide
        
        % Avg across all days
        eval([str(g,:),'allmean_riprateep' num2str(ep) '= nanmean(nanmean(',str(g,:),'riprate_ep',num2str(ep),'));']) % Across animals, then across days
        %eval([str(g,:),'allerr_riprateep' num2str(ep) '= nansem(nanmean(',str(g,:),'riprate_ep',num2str(ep),',2));']) % Across days - mean for each animal - then sem of that
        % Take sem across days and animals: each measure is independent
        eval([str(g,:),'allerr_riprateep' num2str(ep) '= nansem(',str(g,:),'riprate_ep',num2str(ep),'(:));']) % Across days and animals - sem of whole thing
        
        % Avg across novel days
        %eval([str(g,:),'novmean_riprateep' num2str(ep) '= nanmean(nanmean(',str(g,:),'riprate_ep',num2str(ep),'(:,nov)));'])
        %eval([str(g,:),'noverr_riprateep' num2str(ep) '= nansem(nanmean(',str(g,:),'riprate_ep',num2str(ep),'(:,nov),2));']) % Across days - mean for each animal - then sem of that
        % Change - like above, sem of all data together
        % Avg across fam days
        %eval([str(g,:),'fammean_riprateep' num2str(ep) '= nanmean(nanmean(',str(g,:),'riprate_ep',num2str(ep),'(:,fam)));'])
        %eval([str(g,:),'famerr_riprateep' num2str(ep) '= nansem(nanmean(',str(g,:),'riprate_ep',num2str(ep),'(:,fam),2));']) % Across days - mean for each animal - then sem of that
        % Change - like above, sem of all data together
        
        % Get for each day
        for d=1:length(days)
            currday=days(d);
            eval([str(g,:),'daymean_riprateep',num2str(ep),'(',num2str(currday),')= nanmean(',str(g,:),'riprate_ep',num2str(ep),'(:,d));']);
            eval([str(g,:),'dayerr_riprateep',num2str(ep),'(',num2str(currday),')= nansem(',str(g,:),'riprate_ep',num2str(ep),'(:,d));']);
        end
        
        
        % For Ripple Percent
        %-------------------
        
        %eval([str(g,:),'frripper_ep',num2str(ep),' =',str(g,:),'ripper_ep',num2str(ep),'./',str(g,:),'ripper_ep1;']); % Fraction change: Matrix divide
        
         % Avg across all days
        eval([str(g,:),'allmean_ripperep' num2str(ep) '= nanmean(nanmean(',str(g,:),'ripper_ep',num2str(ep),'));']) % Across animals, then across days
        %eval([str(g,:),'allerr_ripperep' num2str(ep) '= nansem(nanmean(',str(g,:),'ripper_ep',num2str(ep),',2));']) % Across days - mean for each animal - then sem of that
        % Take sem across days and animals: each measure is independent
        eval([str(g,:),'allerr_ripperep' num2str(ep) '= nansem(',str(g,:),'ripper_ep',num2str(ep),'(:));']) % Across days and animals - sem of whole thing
        
        % Avg across novel days
        %eval([str(g,:),'novmean_ripperep' num2str(ep) '= nanmean(nanmean(',str(g,:),'ripper_ep',num2str(ep),'(:,nov)));'])
        %eval([str(g,:),'noverr_ripperep' num2str(ep) '= nansem(nanmean(',str(g,:),'ripper_ep',num2str(ep),'(:,nov),2));']) % Across days - mean for each animal - then sem of that
        % Change - like above, sem of all data together
        % Avg across fam days
        %eval([str(g,:),'fammean_ripperep' num2str(ep) '= nanmean(nanmean(',str(g,:),'ripper_ep',num2str(ep),'(:,fam)));'])
        %eval([str(g,:),'famerr_ripperep' num2str(ep) '= nansem(nanmean(',str(g,:),'ripper_ep',num2str(ep),'(:,fam),2));']) % Across days - mean for each animal - then sem of that
        % Change - like above, sem of all data together
        
        
        % For Still Time
        %-------------------
        
        % Avg across all days
        eval([str(g,:),'allmean_stilltimeep' num2str(ep) '= nanmean(nanmean(',str(g,:),'stilltime_ep',num2str(ep),'));']) % Across animals, then across days
        %eval([str(g,:),'allerr_stilltimeep' num2str(ep) '= nansem(nanmean(',str(g,:),'stilltime_ep',num2str(ep),',2));']) % Across days - mean for each animal - then sem of that
        % Take sem across days and animals: each measure is independent
        eval([str(g,:),'allerr_stilltimeep' num2str(ep) '= nansem(',str(g,:),'stilltime_ep',num2str(ep),'(:));']) % Across days and animals - sem of whole thing
        
        % Avg across novel days
        %eval([str(g,:),'novmean_stilltimeep' num2str(ep) '= nanmean(nanmean(',str(g,:),'stilltime_ep',num2str(ep),'(:,nov)));'])
        %eval([str(g,:),'noverr_stilltimeep' num2str(ep) '= nansem(nanmean(',str(g,:),'stilltime_ep',num2str(ep),'(:,nov),2));']) % Across days - mean for each animal - then sem of that
        % Change - like above, sem of all data together
        % Avg across fam days
        %eval([str(g,:),'fammean_stilltimeep' num2str(ep) '= nanmean(nanmean(',str(g,:),'stilltime_ep',num2str(ep),'(:,fam)));'])
        %eval([str(g,:),'famerr_stilltimeep' num2str(ep) '= nansem(nanmean(',str(g,:),'stilltime_ep',num2str(ep),'(:,fam),2));']) % Across days - mean for each animal - then sem of that
        % Change - like above, sem of all data together
        
    end
    
    
    % Get means, etc of Fraction change vs days
    
    postsleeps = [3,5];
    for ep = postsleeps
        
        % For Ripple Rate
        % --------------
        
        % Per day
        eval([str(g,:),'meanfr_riprateep' num2str(ep) '= nanmean(',str(g,:),'frriprate_ep',num2str(ep),');']) % Per day - across animals
        eval([str(g,:),'errfr_riprateep' num2str(ep) '= nansem(',str(g,:),'frriprate_ep',num2str(ep),');']) % Per day - across animals
        % Assign to day... also
        eval([str(g,:),'daymeanfr_riprateep' num2str(ep),'(',num2str(currday),')= nanmean(',str(g,:),'frriprate_ep',num2str(ep),');']) % Per day - across animals
        eval([str(g,:),'dayerrfr_riprateep' num2str(ep),'(',num2str(currday),')= nansem(',str(g,:),'frriprate_ep',num2str(ep),');']) % Per day - across animals
        
        % Avg across all days
        eval([str(g,:),'allmeanfr_riprateep' num2str(ep) '= nanmean(nanmean(',str(g,:),'frriprate_ep',num2str(ep),'));']) % Across animals, then across days
        %eval([str(g,:),'allerrfr_riprateep' num2str(ep) '= nansem(nanmean(',str(g,:),'frriprate_ep',num2str(ep),',2));']) % Across days - mean for each animal - then sem of that
        % Take sem across days and animals: each measure is independent
        eval([str(g,:),'allerrfr_riprateep' num2str(ep) '= nansem(',str(g,:),'frriprate_ep',num2str(ep),'(:));']) % Across days and animals - sem of whole thing
               
%         % Avg across novel days
%         eval([str(g,:),'novmeanfr_riprateep' num2str(ep) '= nanmean(nanmean(',str(g,:),'frriprate_ep',num2str(ep),'(:,nov)));'])
%         eval([str(g,:),'noverrfr_riprateep' num2str(ep) '= nansem(nanmean(',str(g,:),'frriprate_ep',num2str(ep),'(:,nov),2));']) % Across days - mean for each animal - then sem of that
%         % Change - like above, sem of all data together
%         % Avg across fam days
%         eval([str(g,:),'fammeanfr_riprateep' num2str(ep) '= nanmean(nanmean(',str(g,:),'frriprate_ep',num2str(ep),'(:,fam)));'])
%         eval([str(g,:),'famerrfr_riprateep' num2str(ep) '= nansem(nanmean(',str(g,:),'frriprate_ep',num2str(ep),'(:,fam),2));']) % Across days - mean for each animal - then sem of that       
      

        % For Ripple Percent
        % --------------
        
        % Per day
        eval([str(g,:),'meanfr_ripperep' num2str(ep) '= nanmean(',str(g,:),'frripper_ep',num2str(ep),');']) % Per day - across animals
        eval([str(g,:),'errfr_ripperep' num2str(ep) '= nansem(',str(g,:),'frripper_ep',num2str(ep),');']) % Per day - across animals
        % Assign to day... also
        eval([str(g,:),'daymeanfr_ripperep' num2str(ep),'(',num2str(currday),')= nanmean(',str(g,:),'frripper_ep',num2str(ep),');']) % Per day - across animals
        eval([str(g,:),'dayerrfr_ripperep' num2str(ep),'(',num2str(currday),')= nansem(',str(g,:),'frripper_ep',num2str(ep),');']) % Per day - across animals
        
        % Avg across all days
        eval([str(g,:),'allmeanfr_ripperep' num2str(ep) '= nanmean(nanmean(',str(g,:),'frripper_ep',num2str(ep),'));']) % Across animals, then across days
        %eval([str(g,:),'allerrfr_ripperep' num2str(ep) '= nansem(nanmean(',str(g,:),'frripper_ep',num2str(ep),',2));']) % Across days - mean for each animal - then sem of that
        % Take sem across days and animals: each measure is independent
        eval([str(g,:),'allerrfr_riprateep' num2str(ep) '= nansem(',str(g,:),'frripper_ep',num2str(ep),'(:));']) % Across days and animals - sem of whole thing
        
%         % Avg across novel days
%         eval([str(g,:),'novmeanfr_ripperep' num2str(ep) '= nanmean(nanmean(',str(g,:),'frripper_ep',num2str(ep),'(:,nov)));'])
%         eval([str(g,:),'noverrfr_ripperep' num2str(ep) '= nansem(nanmean(',str(g,:),'frripper_ep',num2str(ep),'(:,nov),2));']) % Across days - mean for each animal - then sem of that
%         % Change - like above, sem of all data together
%         % Avg across fam days
%         eval([str(g,:),'fammeanfr_ripperep' num2str(ep) '= nanmean(nanmean(',str(g,:),'frripper_ep',num2str(ep),'(:,fam)));'])
%         eval([str(g,:),'famerrfr_ripperep' num2str(ep) '= nansem(nanmean(',str(g,:),'frripper_ep',num2str(ep),'(:,fam),2));']) % Across days - mean for each animal - then sem of that
        
        
        
        
        % Pairs of sliding days - Not necessary: abandon for now
        cnt=0;
        for nd = 1:length(days)-1
            d = days(nd);
            cnt=cnt+1;
            dayscorr = [];
            currdays = [d, d+1];
            eval([str(g,:),'meanfrslide_riprateep' num2str(ep) '(nd)= nanmean(nanmean(',str(g,:),'frriprate_ep',num2str(ep),'(:,currdays),2));']);
            eval([str(g,:),'errfrslide_riprateep' num2str(ep) '(nd)= nansem(nanmean(',str(g,:),'frriprate_ep',num2str(ep),'(:,currdays),2));']);
            startday(nd) = d;
        end
        
    end
    
    
end   % end str = Exp and Con





%****************************************
% Figures -
%*******************************************



if figopt1 == 1
    
    
    % ********************** Ripple Rate ****************************
    % ***************************************************************
    
    %-----------------------------------------------------------
    %-----------------------------------------------------------
    
    % 1) Bar plot for all, novel, familiar days, & 2) Across days  - Mean Ripple Rates
    
    set(0,'defaultaxesfontsize',20);set(0,'defaultaxesfontweight','normal');
    set(0,'defaultaxeslinewidth',2);
    
    % 1) Bar plot for all, novel, familiar days
    
    % Across all days
    
    sleeps = [1,3,5];
    for ep=1:length(sleeps)
        currep=sleeps(ep);
        eval(['allmean_riprate(',num2str(ep),',:) = [Conallmean_riprateep',num2str(currep),';  Expallmean_riprateep',num2str(currep),'];'])
    end
    
    figure; hold on; redimscreen_figforppt1;
    bar(allmean_riprate,'grouped');
    errorbar(0.8,Conallmean_riprateep1,Conallerr_riprateep1,'b');
    errorbar(1.2,Expallmean_riprateep1,Expallerr_riprateep1,'m');
    errorbar(1.8,Conallmean_riprateep3,Conallerr_riprateep3,'b');
    errorbar(2.2,Expallmean_riprateep3,Expallerr_riprateep3,'m');
    errorbar(2.8,Conallmean_riprateep5,Conallerr_riprateep5,'b');
    errorbar(3.2,Expallmean_riprateep5,Expallerr_riprateep5,'m');
    
    % ----------- Stats -------------
    % Comparing animals: each animal has one number: average across days
    
    % 1)
    % Con against Con for different epochs
    [p13Con_allriprate,h13Con_allriprate] = ranksum(nanmean(Conriprate_ep1,2),nanmean(Conriprate_ep3,2));
    [p15Con_allriprate,h15Con_allriprate] = ranksum(nanmean(Conriprate_ep1,2),nanmean(Conriprate_ep5,2));
    if h13Con_allriprate==1,
        mul = sign(Conallmean_riprateep3);
        plot(1.8, (Conallmean_riprateep3+Conallerr_riprateep3+1.2*mul*Conallerr_riprateep3), 'b*','MarkerSize',8);
    end
    if h15Con_allriprate==1,
        mul = sign(Conallmean_riprateep5);
        plot(2.8, (Conallmean_riprateep5+Conallerr_riprateep5+1.2*mul*Conallerr_riprateep5), 'b*','MarkerSize',8);
    end
    
    % 2)
    % Exp against Exp for different epochs
    [p13Exp_allriprate,h13Exp_allriprate] = ranksum(nanmean(Expriprate_ep1,2),nanmean(Expriprate_ep3,2));
    [p15Exp_allriprate,h15Exp_allriprate] = ranksum(nanmean(Expriprate_ep1,2),nanmean(Expriprate_ep5,2));
    if h13Exp_allriprate==1,
        mul = sign(Expallmean_riprateep3);
        plot(2.2, (Expallmean_riprateep3+Expallerr_riprateep3+1.2*mul*Expallerr_riprateep3), 'r*','MarkerSize',8);
    end
    if h15Exp_allriprate==1,
        mul = sign(Expallmean_riprateep5);
        plot(3.2, (Expallmean_riprateep5+Expallerr_riprateep5+1.2*mul*Expallerr_riprateep5), 'r*','MarkerSize',8);
    end
    
    % 3)
    % Con against Exp for each epoch
    [p11_allriprate,h11_allriprate] = ranksum(nanmean(Conriprate_ep1,2),nanmean(Expriprate_ep1,2));
    [p33_allriprate,h33_allriprate] = ranksum(nanmean(Conriprate_ep3,2),nanmean(Expriprate_ep3,2));
    [p55_allriprate,h55_allriprate] = ranksum(nanmean(Conriprate_ep5,2),nanmean(Expriprate_ep5,2));
    
    % Plot midway between the two graphs
    if h11_allriprate==1,
        mul = sign(Conallmean_riprateep1);
        plot(1, (Conallmean_riprateep1+Conallerr_riprateep1+1.2*mul*Conallerr_riprateep1), 'g*','MarkerSize',12);
    end
    
    if h33_allriprate==1,
        mul = sign(Conallmean_riprateep3);
        plot(2, (Conallmean_riprateep3+Conallerr_riprateep3+1.2*mul*Conallerr_riprateep3), 'g*','MarkerSize',12);
    end
    
    if h55_allriprate==1,
        mul = sign(Conallmean_riprateep5);
        plot(3, (Conallmean_riprateep5+Conallerr_riprateep5+1.2*mul*Conallerr_riprateep5), 'g*','MarkerSize',12);
    end
    
    % -------------------------
    
    title('Ripple Rate during Sleep: All Days');
    ylabel('Ripple Rate (Hz)');
    set(gca,'XTick',[1:3],'XTickLabel',{'Pre';'Post1';'Post2'},'FontSize',16,'Fontweight','normal');
    %set(gca,'XLim',[0.5 3.5]);
    %axis([0 4 0 0.5])
    
    
    %-----------------------------------------------------------
    %-----------------------------------------------------------
    
    % Across Novel days
    
    if ~isempty(intersect(days,novel))
        noveldays = intersect(days,novel);
        
        sleeps = [1,3,5];
        for ep=1:length(sleeps)
            currep=sleeps(ep);
            eval(['novmean_riprate(',num2str(ep),',:) = [Connovmean_riprateep',num2str(currep),';  Expnovmean_riprateep',num2str(currep),'];'])
        end
        
        figure; hold on; redimscreen_figforppt1;
        bar(novmean_riprate,'grouped');
        errorbar(0.8,Connovmean_riprateep1,Connoverr_riprateep1,'b');
        errorbar(1.2,Expnovmean_riprateep1,Expnoverr_riprateep1,'m');
        errorbar(1.8,Connovmean_riprateep3,Connoverr_riprateep3,'b');
        errorbar(2.2,Expnovmean_riprateep3,Expnoverr_riprateep3,'m');
        errorbar(2.8,Connovmean_riprateep5,Connoverr_riprateep5,'b');
        errorbar(3.2,Expnovmean_riprateep5,Expnoverr_riprateep5,'m');
        
        
        % ----------- Stats -------------
        % Comparing animals: each animal has one number: average across days
        
        % 1)
        % Con against Con for different epochs
        
        [p13Con_novriprate,h13Con_novriprate] = ranksum(nanmean(Conriprate_ep1(:,nov),2),nanmean(Conriprate_ep3(:,nov),2));
        [p15Con_novriprate,h15Con_novriprate] = ranksum(nanmean(Conriprate_ep1(:,nov),2),nanmean(Conriprate_ep5(:,nov),2));
        if h13Con_novriprate==1,
            mul = sign(Connovmean_riprateep3);
            plot(1.8, (Connovmean_riprateep3+Connoverr_riprateep3+1.2*mul*Connoverr_riprateep3), 'b*','MarkerSize',8);
        end
        if h15Con_novriprate==1,
            mul = sign(Connovmean_riprateep5);
            plot(2.8, (Connovmean_riprateep5+Connoverr_riprateep5+1.2*mul*Connoverr_riprateep5), 'b*','MarkerSize',8);
        end
        
        % 2)
        % Exp against Exp for different epochs
        [p13Exp_novriprate,h13Exp_novriprate] = ranksum(nanmean(Expriprate_ep1(:,nov),2),nanmean(Expriprate_ep3(:,nov),2));
        [p15Exp_novriprate,h15Exp_novriprate] = ranksum(nanmean(Expriprate_ep1(:,nov),2),nanmean(Expriprate_ep5(:,nov),2));
        if h13Exp_novriprate==1,
            mul = sign(Expnovmean_riprateep3);
            plot(2.2, (Expnovmean_riprateep3+Expnoverr_riprateep3+1.2*mul*Expnoverr_riprateep3), 'r*','MarkerSize',8);
        end
        if h15Exp_novriprate==1,
            mul = sign(Expnovmean_riprateep5);
            plot(3.2, (Expnovmean_riprateep5+Expnoverr_riprateep5+1.2*mul*Expnoverr_riprateep5), 'r*','MarkerSize',8);
        end
        
        % 3)
        % Con against Exp for each epoch
        [p11_novriprate,h11_novriprate] = ranksum(nanmean(Conriprate_ep1(:,nov),2),nanmean(Expriprate_ep1(:,nov),2));
        [p33_novriprate,h33_novriprate] = ranksum(nanmean(Conriprate_ep3(:,nov),2),nanmean(Expriprate_ep3(:,nov),2));
        [p55_novriprate,h55_novriprate] = ranksum(nanmean(Conriprate_ep5(:,nov),2),nanmean(Expriprate_ep5(:,nov),2));
        
        % Plot midway between the two graphs
        if h11_novriprate==1,
            mul = sign(Connovmean_riprateep1);
            plot(1, (Connovmean_riprateep1+Connoverr_riprateep1+1.2*mul*Connoverr_riprateep1), 'g*','MarkerSize',12);
        end
        
        if h33_novriprate==1,
            mul = sign(Connovmean_riprateep3);
            plot(2, (Connovmean_riprateep3+Connoverr_riprateep3+1.2*mul*Connoverr_riprateep3), 'g*','MarkerSize',12);
        end
        
        if h55_novriprate==1,
            mul = sign(Connovmean_riprateep5);
            plot(3, (Connovmean_riprateep5+Connoverr_riprateep5+1.2*mul*Connoverr_riprateep5), 'g*','MarkerSize',12);
        end
        
        % -------------------------
        
        
        title('Ripple Rate during Sleep: Novel Days');
        ylabel('Ripple Rate (Hz)');
        set(gca,'XTick',[1:3],'XTickLabel',{'Pre';'Post1';'Post2'},'FontSize',16,'Fontweight','normal');
    end
    
    
    
    %-----------------------------------------------------------
    %-----------------------------------------------------------
    
    % Across Familiar days
    
    if ~isempty(intersect(days,fam))
        famdays = intersect(days,fam);
        
        sleeps = [1,3,5];
        for ep=1:length(sleeps)
            currep=sleeps(ep);
            eval(['fammean_riprate(',num2str(ep),',:) = [Confammean_riprateep',num2str(currep),';  Expfammean_riprateep',num2str(currep),'];'])
        end
        
        figure; hold on; redimscreen_figforppt1;
        bar(fammean_riprate,'grouped');
        errorbar(0.8,Confammean_riprateep1,Confamerr_riprateep1,'b');
        errorbar(1.2,Expfammean_riprateep1,Expfamerr_riprateep1,'m');
        errorbar(1.8,Confammean_riprateep3,Confamerr_riprateep3,'b');
        errorbar(2.2,Expfammean_riprateep3,Expfamerr_riprateep3,'m');
        errorbar(2.8,Confammean_riprateep5,Confamerr_riprateep5,'b');
        errorbar(3.2,Expfammean_riprateep5,Expfamerr_riprateep5,'m');
        
        % ----------- Stats -------------
        % Comparing animals: each animal has one number: average across days
        
        % 1)
        % Con against Con for different epochs
        
        [p13Con_famriprate,h13Con_famriprate] = ranksum(nanmean(Conriprate_ep1(:,fam),2),nanmean(Conriprate_ep3(:,fam),2));
        [p15Con_famriprate,h15Con_famriprate] = ranksum(nanmean(Conriprate_ep1(:,fam),2),nanmean(Conriprate_ep5(:,fam),2));
        if h13Con_famriprate==1,
            mul = sign(Confammean_riprateep3);
            plot(1.8, (Confammean_riprateep3+Confamerr_riprateep3+1.2*mul*Confamerr_riprateep3), 'b*','MarkerSize',8);
        end
        if h15Con_famriprate==1,
            mul = sign(Confammean_riprateep5);
            plot(2.8, (Confammean_riprateep5+Confamerr_riprateep5+1.2*mul*Confamerr_riprateep5), 'b*','MarkerSize',8);
        end
        
        % 2)
        % Exp against Exp for different epochs
        [p13Exp_famriprate,h13Exp_famriprate] = ranksum(nanmean(Expriprate_ep1(:,fam),2),nanmean(Expriprate_ep3(:,fam),2));
        [p15Exp_famriprate,h15Exp_famriprate] = ranksum(nanmean(Expriprate_ep1(:,fam),2),nanmean(Expriprate_ep5(:,fam),2));
        if h13Exp_famriprate==1,
            mul = sign(Expfammean_riprateep3);
            plot(2.2, (Expfammean_riprateep3+Expfamerr_riprateep3+1.2*mul*Expfamerr_riprateep3), 'r*','MarkerSize',8);
        end
        if h15Exp_famriprate==1,
            mul = sign(Expfammean_riprateep5);
            plot(3.2, (Expfammean_riprateep5+Expfamerr_riprateep5+1.2*mul*Expfamerr_riprateep5), 'r*','MarkerSize',8);
        end
        
        % 3)
        % Con against Exp for each epoch
        [p11_famriprate,h11_famriprate] = ranksum(nanmean(Conriprate_ep1(:,fam),2),nanmean(Expriprate_ep1(:,fam),2));
        [p33_famriprate,h33_famriprate] = ranksum(nanmean(Conriprate_ep3(:,fam),2),nanmean(Expriprate_ep3(:,fam),2));
        [p55_famriprate,h55_famriprate] = ranksum(nanmean(Conriprate_ep5(:,fam),2),nanmean(Expriprate_ep5(:,fam),2));
        
        % Plot midway between the two graphs
        if h11_famriprate==1,
            mul = sign(Confammean_riprateep1);
            plot(1, (Confammean_riprateep1+Confamerr_riprateep1+1.2*mul*Confamerr_riprateep1), 'g*','MarkerSize',12);
        end
        
        if h33_famriprate==1,
            mul = sign(Confammean_riprateep3);
            plot(2, (Confammean_riprateep3+Confamerr_riprateep3+1.2*mul*Confamerr_riprateep3), 'g*','MarkerSize',12);
        end
        
        if h55_famriprate==1,
            mul = sign(Confammean_riprateep5);
            plot(3, (Confammean_riprateep5+Confamerr_riprateep5+1.2*mul*Confamerr_riprateep5), 'g*','MarkerSize',12);
        end
        
        % -------------------------
        
        title('Ripple Rate during Sleep: Familiar Days');
        ylabel('Ripple Rate (Hz)');
        set(gca,'XTick',[1:3],'XTickLabel',{'Pre';'Post1';'Post2'},'FontSize',16,'Fontweight','normal');
    end
    
    
    
    
    
    
    
    
    
    % ************ Fraction Change In Ripple Rate *******************
    % ***************************************************************
    
    %-----------------------------------------------------------
    %-----------------------------------------------------------
    
    % 1) Bar plot for all, novel, familiar days, & 2) Across days
    
    set(0,'defaultaxesfontsize',20);set(0,'defaultaxesfontweight','normal');
    set(0,'defaultaxeslinewidth',2);
    
    % 1) Bar plot for all, novel, familiar days
    
    % Across all days
    
    sleeps = [3,5];
    for ep=1:length(sleeps)
        currep=sleeps(ep);
        eval(['allmeanfr_riprate(',num2str(ep),',:) = [Conallmeanfr_riprateep',num2str(currep),';  Expallmeanfr_riprateep',num2str(currep),'];'])
    end
    
    figure; hold on; redimscreen_figforppt1;
    bar([2 3], allmeanfr_riprate,'grouped');
    errorbar(1.8,Conallmeanfr_riprateep3,Conallerrfr_riprateep3,'b');
    errorbar(2.2,Expallmeanfr_riprateep3,Expallerrfr_riprateep3,'m');
    errorbar(2.8,Conallmeanfr_riprateep5,Conallerrfr_riprateep5,'b');
    errorbar(3.2,Expallmeanfr_riprateep5,Expallerrfr_riprateep5,'m');
    
    % ----------- Stats -------------
    % Comparing animals: each animal has one number: average across days
    
    % 1)
    % Con against Con for different epochs
    [p35Con_allripratefr,h35Con_allripratefr] = ranksum(nanmean(Confrriprate_ep3,2),nanmean(Confrriprate_ep5,2));
    if h35Con_allripratefr==1,
        mul = sign(Conallmeanfr_riprateep5);
        plot(2.8, (Conallmeanfr_riprateep5+Conallerrfr_riprateep5+1.2*mul*Conallerrfr_riprateep5), 'b*','MarkerSize',8);
    end
    
    
    % 2)
    % Exp against Exp for different epochs
    [p35Exp_allripratefr,h35Exp_allripratefr] = ranksum(nanmean(Expfrriprate_ep3,2),nanmean(Expfrriprate_ep5,2));
    if h35Exp_allripratefr==1,
        mul = sign(Expallmeanfr_riprateep5);
        plot(3.2, (Expallmeanfr_riprateep5+Expallerrfr_riprateep5+1.2*mul*Expallerrfr_riprateep5), 'r*','MarkerSize',8);
    end
    
    % 3)
    % Con against Exp for each epoch
    [p33_allripratefr,h33_allripratefr] = ranksum(nanmean(Confrriprate_ep3,2),nanmean(Expfrriprate_ep3,2));
    [p55_allripratefr,h55_allripratefr] = ranksum(nanmean(Confrriprate_ep5,2),nanmean(Expfrriprate_ep5,2));
    
    % Plot midway between the two graphs
    
    if h33_allripratefr==1,
        mul = sign(Conallmeanfr_riprateep3);
        plot(2, (Conallmeanfr_riprateep3+Conallerrfr_riprateep3+1.2*mul*Conallerrfr_riprateep3), 'g*','MarkerSize',12);
    end
    
    if h55_allripratefr==1,
        mul = sign(Conallmeanfr_riprateep5);
        plot(3, (Conallmeanfr_riprateep5+Conallerrfr_riprateep5+1.2*mul*Conallerrfr_riprateep5), 'g*','MarkerSize',12);
    end
    
    % -------------------------
    
    title('Fraction Change in Ripple Rate during Sleep: All Days');
    ylabel('Fraction Change');
    set(gca,'XTick',[2:3],'XTickLabel',{'Post1';'Post2'},'FontSize',16,'Fontweight','normal');
    %set(gca,'XLim',[0.5 3.5]);
    %axis([0 4 0 0.5])
    
    
    
    %-----------------------------------------------------------
    %-----------------------------------------------------------
    
    % Across Novel days
    
    sleeps = [3,5];
    for ep=1:length(sleeps)
        currep=sleeps(ep);
        eval(['novmeanfr_riprate(',num2str(ep),',:) = [Connovmeanfr_riprateep',num2str(currep),';  Expnovmeanfr_riprateep',num2str(currep),'];'])
    end
    
    figure; hold on; redimscreen_figforppt1;
    bar([2 3], novmeanfr_riprate,'grouped');
    errorbar(1.8,Connovmeanfr_riprateep3,Connoverrfr_riprateep3,'b');
    errorbar(2.2,Expnovmeanfr_riprateep3,Expnoverrfr_riprateep3,'m');
    errorbar(2.8,Connovmeanfr_riprateep5,Connoverrfr_riprateep5,'b');
    errorbar(3.2,Expnovmeanfr_riprateep5,Expnoverrfr_riprateep5,'m');
    
    % ----------- Stats -------------
    % Comparing animals: each animal has one number: average across days
    
    % 1)
    % Con against Con for different epochs
    [p35Con_novripratefr,h35Con_novripratefr] = ranksum(nanmean(Confrriprate_ep3(:,nov),2),nanmean(Confrriprate_ep5(:,nov),2));
    if h35Con_novripratefr==1,
        mul = sign(Connovmeanfr_riprateep5);
        plot(2.8, (Connovmeanfr_riprateep5+Connoverrfr_riprateep5+1.2*mul*Connoverrfr_riprateep5), 'b*','MarkerSize',8);
    end
    
    
    % 2)
    % Exp against Exp for different epochs
    [p35Exp_novripratefr,h35Exp_novripratefr] = ranksum(nanmean(Expfrriprate_ep3(:,nov),2),nanmean(Expfrriprate_ep5(:,nov),2));
    if h35Exp_novripratefr==1,
        mul = sign(Expnovmeanfr_riprateep5);
        plot(3.2, (Expnovmeanfr_riprateep5+Expnoverrfr_riprateep5+1.2*mul*Expnoverrfr_riprateep5), 'r*','MarkerSize',8);
    end
    
    % 3)
    % Con against Exp for each epoch
    [p33_novripratefr,h33_novripratefr] = ranksum(nanmean(Confrriprate_ep3(:,nov),2),nanmean(Expfrriprate_ep3(:,nov),2));
    [p55_novripratefr,h55_novripratefr] = ranksum(nanmean(Confrriprate_ep5(:,nov),2),nanmean(Expfrriprate_ep5(:,nov),2));
    
    % Plot midway between the two graphs
    
    if h33_novripratefr==1,
        mul = sign(Connovmeanfr_riprateep3);
        plot(2, (Connovmeanfr_riprateep3+Connoverrfr_riprateep3+1.2*mul*Connoverrfr_riprateep3), 'g*','MarkerSize',12);
    end
    
    if h55_novripratefr==1,
        mul = sign(Connovmeanfr_riprateep5);
        plot(3, (Connovmeanfr_riprateep5+Connoverrfr_riprateep5+1.2*mul*Connoverrfr_riprateep5), 'g*','MarkerSize',12);
    end
    
    % -------------------------
    
    title('Fraction Change in Ripple Rate during Sleep: Nov Days');
    ylabel('Fraction Change');
    set(gca,'XTick',[2:3],'XTickLabel',{'Post1';'Post2'},'FontSize',16,'Fontweight','normal');
    %set(gca,'XLim',[0.5 3.5]);
    %axis([0 4 0 0.5])
    
    
    
    
    







% ***************  Plots For Ripple Size  ******************
% ***************************************************************

% Ripple  Size
%eval([str(g,:),'ripntet_ep',num2str(currep),'{an}{d} =',str(g,:),'ripf(an).output{1}(index).ripntet;']);
%eval([str(g,:),'ripsize_ep',num2str(currep),'{an}{d} =',str(g,:),'ripf(an).output{1}(index).ripsize;']);
% str=['Exp';'Con'];

if figopt2==1
    
    nov = 1:3;
    fam =4:8;
    
    dayscell{1} = days;
    dayscell{2} = nov;
    dayscell{3} = fam;
    
    strlabel{1} = 'All Days';
    strlabel{2} = 'Nov Days';
    strlabel{3} = 'Fam Days';
    
    ripsizeedges = 0:20;
    ripntetrange = 0:6;
    
    % First calculate
    
    % Either
    % 1) gather data for each animal for all days in a vector and then
    % hist, OR/AND
    % 2) hist for each day. Then take mean across days for each animal
    
 
    
    % Loop over dayscell
    for n = 1:length(dayscell)
        
        currdays = dayscell{n};
        
        
        % Initialize for both groups each time you loop over dayscell
        % Make vector for each epoch for ripsize and ripntet of both groups
        % Could Use this for statistics
        for g = 1:size(str,1)   % Do Exp and Con groups separately
            for ep=1:length(allepochs)
                currep=allepochs(ep);
                eval([str(g,:),'ripsize_ep',num2str(currep),'vec=[];']);
                eval([str(g,:),'ripntet_ep',num2str(currep),'vec=[];']);
                eval([str(g,:),'ripsizehist_ep',num2str(currep),'vec=[];']);
                eval([str(g,:),'ripntethist_ep',num2str(currep),'vec=[];']);
                eval([str(g,:),'ripsizehist_ep',num2str(currep),'=[];']);
                eval([str(g,:),'ripntethist_ep',num2str(currep),'=[];']);
            end
        end
        
        
        for g = 1:size(str,1)   % Do Exp and Con groups separately
            
            % Get nanim for current grp
            eval(['totanim = length(',str(g,:),'ripntet_ep1);']);
            
            for an=1:totanim % Loop over anim
                
                % Initialize for each animal
                for ep=1:length(allepochs)
                    currep=allepochs(ep);
                    eval(['ripsizehist_ep',num2str(ep),'= [];']);
                    eval(['ripntethist_ep',num2str(ep),'= [];']);
                end
                
                for dy = 1:length(currdays)
                    d = currdays(dy);
                    for ep=1:length(allepochs)
                        
                        currep=allepochs(ep);
                        eval(['currripsize =',str(g,:),'ripsize_ep',num2str(currep),'{an}{d};']);
                        eval(['currripntet =',str(g,:),'ripntet_ep',num2str(currep),'{an}{d};']);
                        
                        % Store all raw values in vector
                        eval([str(g,:),'ripsize_ep',num2str(currep),'vec=[',str(g,:),'ripsize_ep',num2str(currep),'vec, currripsize];']);
                        eval([str(g,:),'ripntet_ep',num2str(currep),'vec=[',str(g,:),'ripntet_ep',num2str(currep),'vec, currripntet];']);
                        
                        % Make histogram
                        h = histc(currripsize,ripsizeedges);
                        h = cumsum(h); h = h./max(h);  % Cumulative proportion for current day and epoch
                        eval(['ripsizehist_ep',num2str(currep),'(dy,:) = h;']);
                        
                        hn = histc(currripntet,ripntetrange);
                        hn = cumsum(hn); hn = hn./max(hn);
                        eval(['ripntethist_ep',num2str(currep),'(dy,:) = hn;']); % Simple hist
                        
                        % Store histogram in vector
                        eval([str(g,:),'ripsizehist_ep',num2str(currep),'vec=[',str(g,:),'ripsizehist_ep',num2str(currep),'vec, h];']);
                        eval([str(g,:),'ripntethist_ep',num2str(currep),'vec=[',str(g,:),'ripntethist_ep',num2str(currep),'vec, hn];']);
                        
                    end % end epoch
                end % end day
                
                % Unpack for each animal and a) average across days for plotting
                for ep=1:length(allepochs)
                    currep=allepochs(ep);
                    eval([str(g,:),'ripsizehist_ep',num2str(currep),'(an,:) = mean(ripsizehist_ep',num2str(currep),',1);']);
                    eval([str(g,:),'ripntethist_ep',num2str(currep),'(an,:) = mean(ripntethist_ep',num2str(currep),',1);']);
                end
            end % end anim
            
        end % Loop over str: Exp and Con
        
        
        % Plotting for current days cell
        
        
        %%%%%%%%%%%%%%% RipSize Distribution %%%%%%%%%
        
        if figoptsize ==1
            
            % 1) Control Group: All epochs
            
            figure; hold on; redimscreen_figforppt1;
            
            plot(ripsizeedges,mean(Conripsizehist_ep1),'k.-','Linewidth',2,'MarkerSize',18);
            jbfill(ripsizeedges, mean(Conripsizehist_ep1)+sem(Conripsizehist_ep1),...
                mean(Conripsizehist_ep1)-sem(Conripsizehist_ep1),'k','k',1,0.3);
            
            %         plot(ripsizeedges,mean(Conripsizehist_ep3),'c.-','Linewidth',2,'MarkerSize',18);
            %         jbfill(ripsizeedges, mean(Conripsizehist_ep3)+sem(Conripsizehist_ep3),...
            %             mean(Conripsizehist_ep3)-sem(Conripsizehist_ep3),'c','c',1,0.3);
            
            plot(ripsizeedges,mean(Conripsizehist_ep5),'b.-','Linewidth',2,'MarkerSize',18);
            jbfill(ripsizeedges, mean(Conripsizehist_ep5)+sem(Conripsizehist_ep5),...
                mean(Conripsizehist_ep5)-sem(Conripsizehist_ep5),'b','b',1,0.3);
            
            %[h15Con,p15Con] = kstest2(Conripsize_ep1vec,Conripsize_ep5vec); % All raw values
            %[h15Con,p15Con] = kstest2(Conripsizehist_ep1(:),Conripsizehist_ep5(:)); % nanim*bins; 1 distribn for all anim - avg across all days
            [h15Con,p15Con] = kstest2(Conripsizehist_ep1vec,Conripsizehist_ep5vec); % nanim*days*bins; 1 distribn for each day across animals
            if h15Con == 1
                plot(10, 0.5, 'r*','MarkerSize',12);
            end
            text(10,0.4,['p = ',num2str(p15Con)],'FontSize',16);
            
            set(gca,'XLim',[1 16]);
            set(gca,'YLim',[0 1.02]);
            
            title(['Con - Ripple Size during Sleep: ',strlabel{n}],'FontSize',24,'Fontweight','normal');
            ylabel('Cumulative Proportion','FontSize',24,'Fontweight','normal');
            xlabel('Ripple Size (stdev)','FontSize',24,'Fontweight','normal')
            
            
            % 2) Exp Group: All epochs
            
            figure; hold on; redimscreen_figforppt1;
            
            plot(ripsizeedges,mean(Expripsizehist_ep1),'k.-','Linewidth',2,'MarkerSize',18);
            jbfill(ripsizeedges, mean(Expripsizehist_ep1)+sem(Expripsizehist_ep1),...
                mean(Expripsizehist_ep1)-sem(Expripsizehist_ep1),'k','k',1,0.3);
       
            
%                     plot(ripsizeedges,mean(Expripsizehist_ep3),'m.-','Linewidth',2,'MarkerSize',18);
%                     jbfill(ripsizeedges, mean(Expripsizehist_ep3)+sem(Expripsizehist_ep3),...
%                         mean(Expripsizehist_ep3)-sem(Expripsizehist_ep3),'m','m',1,0.3);
            
            plot(ripsizeedges,mean(Expripsizehist_ep5),'r.-','Linewidth',2,'MarkerSize',18);
            jbfill(ripsizeedges, mean(Expripsizehist_ep5)+sem(Expripsizehist_ep5),...
                mean(Expripsizehist_ep5)-sem(Expripsizehist_ep5),'r','r',1,0.3);
            
            [h15Exp,p15Exp] = kstest2(Conripsizehist_ep1vec,Expripsizehist_ep5vec); h15Exp;
            if h15Exp == 1
                plot(10, 0.5, 'r*','MarkerSize',12);
            end
            text(10,0.4,['p = ',num2str(p15Exp)],'FontSize',16);
            
            set(gca,'XLim',[1 16]);
            set(gca,'YLim',[0 1.02]);
            
            title(['Exp - Ripple Size during Sleep: ',strlabel{n}],'FontSize',24,'Fontweight','normal');
            ylabel('Cumulative Proportion','FontSize',24,'Fontweight','normal');
            xlabel('Ripple Size (stdev)','FontSize',24,'Fontweight','normal')
            
            
            % 3) Con vs. Exp : Epoch 1
            
            figure; hold on; redimscreen_figforppt1;
            
            plot(ripsizeedges,mean(Conripsizehist_ep1),'b.-','Linewidth',2,'MarkerSize',18);
            jbfill(ripsizeedges, mean(Conripsizehist_ep1)+sem(Conripsizehist_ep1),...
                mean(Conripsizehist_ep1)-sem(Conripsizehist_ep1),'b','b',1,0.3);
            
            plot(ripsizeedges,mean(Expripsizehist_ep1),'r.-','Linewidth',2,'MarkerSize',18);
            jbfill(ripsizeedges, mean(Expripsizehist_ep1)+sem(Expripsizehist_ep1),...
                mean(Expripsizehist_ep1)-sem(Expripsizehist_ep1),'r','r',1,0.3);
            
            set(gca,'XLim',[1 16]);
            set(gca,'YLim',[0 1.02]);
            
            [h11,p11] = kstest2(Conripsizehist_ep1vec,Expripsizehist_ep1vec); h11;
            if h11 == 1
                plot(10, 0.5, 'r*','MarkerSize',12);
            end
            text(10,0.4,['p = ',num2str(p11)],'FontSize',16);
            
            title(['Con vs Exp - Ripple Size during Pre-Sleep: ',strlabel{n}],'FontSize',24,'Fontweight','normal');
            ylabel('Cumulative Proportion','FontSize',24,'Fontweight','normal');
            xlabel('Ripple Size (stdev)','FontSize',24,'Fontweight','normal')
            
            
            % 4) Con vs. Exp : Epoch 5
            
            figure; hold on; redimscreen_figforppt1;
            
            plot(ripsizeedges,mean(Conripsizehist_ep5),'b.-','Linewidth',2,'MarkerSize',18);
            jbfill(ripsizeedges, mean(Conripsizehist_ep5)+sem(Conripsizehist_ep5),...
                mean(Conripsizehist_ep5)-sem(Conripsizehist_ep5),'b','b',1,0.3);
            
            plot(ripsizeedges,mean(Expripsizehist_ep5),'r.-','Linewidth',2,'MarkerSize',18);
            jbfill(ripsizeedges, mean(Expripsizehist_ep5)+sem(Expripsizehist_ep5),...
                mean(Expripsizehist_ep5)-sem(Expripsizehist_ep5),'r','r',1,0.3);
            
            [h55,p55] = kstest2(Conripsizehist_ep5vec,Expripsizehist_ep5vec); h55;
            if h55 == 1
                plot(10, 0.5, 'r*','MarkerSize',12);
            end
            text(10,0.4,['p = ',num2str(p55)],'FontSize',16);
            
            set(gca,'XLim',[1 16]);
            set(gca,'YLim',[0 1.02]);
            
            title(['Con vs Exp - Ripple Size during Post-Sleep: ',strlabel{n}],'FontSize',24,'Fontweight','normal');
            ylabel('Cumulative Proportion','FontSize',24,'Fontweight','normal');
            xlabel('Ripple Size (stdev)','FontSize',24,'Fontweight','normal')
            
        end
        
        
        
        %%%%%%%%%%%%%%% RipNtet Distribution %%%%%%%%%
        
        %%% Need to change Stat Distribution for this
        
        % 1) Control Group: All epochs
        
        if figoptntet == 1
            
            figure; hold on; redimscreen_figforppt1;
            
            plot(ripntetrange,mean(Conripntethist_ep1),'k.-','Linewidth',2,'MarkerSize',24);
            jbfill(ripntetrange, mean(Conripntethist_ep1)+sem(Conripntethist_ep1),...
                mean(Conripntethist_ep1)-sem(Conripntethist_ep1),'k','k',1,0.3);
            
            %         plot(ripntetrange,mean(Conripntethist_ep3),'c.-','Linewidth',2,'MarkerSize',24);
            %         jbfill(ripntetrange, mean(Conripntethist_ep3)+sem(Conripntethist_ep3),...
            %             mean(Conripntethist_ep3)-sem(Conripntethist_ep3),'c','c',1,0.3);
            
            plot(ripntetrange,mean(Conripntethist_ep5),'b.-','Linewidth',2,'MarkerSize',24);
            jbfill(ripntetrange, mean(Conripntethist_ep5)+sem(Conripntethist_ep5),...
                mean(Conripntethist_ep5)-sem(Conripntethist_ep5),'b','b',1,0.3);
            
            [h15Con,p15Con] = kstest2(Conripntethist_ep1vec,Conripntethist_ep5vec); h15Con;
            if h15Con == 1
                plot(4, 0.4, 'r*','MarkerSize',12);
            end
            text(4, 0.3,['p = ',num2str(p15Con)],'FontSize',16);
            
            
            %         set(gca,'XLim',[1 16]);
            %         set(gca,'YLim',[0 1.02]);
            
            title(['Con - Ripple Ntet during Sleep: ',strlabel{n}],'FontSize',24,'Fontweight','normal');
            ylabel('Proportion','FontSize',24,'Fontweight','normal');
            xlabel('Ripple Size (stdev)','FontSize',24,'Fontweight','normal')
            
            
            % 2) Exp Group: All epochs
            
            figure; hold on; redimscreen_figforppt1;
            
            plot(ripntetrange,mean(Expripntethist_ep1),'k.-','Linewidth',2,'MarkerSize',24);
            jbfill(ripntetrange, mean(Expripntethist_ep1)+sem(Expripntethist_ep1),...
                mean(Expripntethist_ep1)-sem(Expripntethist_ep1),'k','k',1,0.3);
            
            %         plot(ripntetrange,mean(Expripntethist_ep3),'m.-','Linewidth',2,'MarkerSize',24);
            %         jbfill(ripntetrange, mean(Expripntethist_ep3)+sem(Expripntethist_ep3),...
            %             mean(Expripntethist_ep3)-sem(Expripntethist_ep3),'m','m',1,0.3);
            
            plot(ripntetrange,mean(Expripntethist_ep5),'r.-','Linewidth',2,'MarkerSize',24);
            jbfill(ripntetrange, mean(Expripntethist_ep5)+sem(Expripntethist_ep5),...
                mean(Expripntethist_ep5)-sem(Expripntethist_ep5),'r','r',1,0.3);
            
            [h15Exp,p15Exp] = kstest2(Expripntethist_ep1vec,Expripntethist_ep5vec); h15Exp;
            if h15Exp == 1
                plot(4, 0.4, 'r*','MarkerSize',12);
            end
            text(4, 0.3,['p = ',num2str(p15Exp)],'FontSize',16);
            
            
            title(['Exp - Ripple Ntet during Sleep: ',strlabel{n}],'FontSize',24,'Fontweight','normal');
            ylabel('Proportion','FontSize',24,'Fontweight','normal');
            xlabel('Ripple Size (stdev)','FontSize',24,'Fontweight','normal')
            
            
            % 3) Con vs. Exp : Epoch 1
            
            figure; hold on; redimscreen_figforppt1;
            
            plot(ripntetrange,mean(Conripntethist_ep1),'b.-','Linewidth',2,'MarkerSize',24);
            jbfill(ripntetrange, mean(Conripntethist_ep1)+sem(Conripntethist_ep1),...
                mean(Conripntethist_ep1)-sem(Conripntethist_ep1),'b','b',1,0.3);
            
            plot(ripntetrange,mean(Expripntethist_ep1),'r.-','Linewidth',2,'MarkerSize',24);
            jbfill(ripntetrange, mean(Expripntethist_ep1)+sem(Expripntethist_ep1),...
                mean(Expripntethist_ep1)-sem(Expripntethist_ep1),'r','r',1,0.3);
            
            
            [h11,p11] = kstest2(Conripntethist_ep1vec,Expripntethist_ep1vec); h11;
            if h11 == 1
                plot(4, 0.4, 'r*','MarkerSize',12);
            end
            text(4, 0.3,['p = ',num2str(p11)],'FontSize',16);
            
            
            title(['Con vs Exp - Ripple Ntet during Pre-Sleep: ',strlabel{n}],'FontSize',24,'Fontweight','normal');
            ylabel('Proportion','FontSize',24,'Fontweight','normal');
            xlabel('Ripple Size (stdev)','FontSize',24,'Fontweight','normal')
            
            
            % 4) Con vs. Exp : Epoch 5
            
            figure; hold on; redimscreen_figforppt1;
            
            plot(ripntetrange,mean(Conripntethist_ep5),'b.-','Linewidth',2,'MarkerSize',24);
            jbfill(ripntetrange, mean(Conripntethist_ep5)+sem(Conripntethist_ep5),...
                mean(Conripntethist_ep5)-sem(Conripntethist_ep5),'b','b',1,0.3);
            
            plot(ripntetrange,mean(Expripntethist_ep5),'r.-','Linewidth',2,'MarkerSize',24);
            jbfill(ripntetrange, mean(Expripntethist_ep5)+sem(Expripntethist_ep5),...
                mean(Expripntethist_ep5)-sem(Expripntethist_ep5),'r','r',1,0.3);
            
            [h55,p55] = kstest2(Conripntethist_ep5vec,Expripntethist_ep5vec); h55;
            if h55 == 1
                plot(4, 0.4, 'r*','MarkerSize',12);
            end
            text(4, 0.3,['p = ',num2str(p55)],'FontSize',16);
            
            
            title(['Con vs Exp - Ripple Ntet during Post-Sleep: ',strlabel{n}],'FontSize',24,'Fontweight','normal');
            ylabel('Proportion','FontSize',24,'Fontweight','normal');
            xlabel('Ripple Size (stdev)','FontSize',24,'Fontweight','normal')
            
        end % figoptntet
        
        
    end % Loop over dayscell
    
end % if figopt2











    
    
    
    
    
    
    %-----------------------------------------------------------
    %-----------------------------------------------------------
    
    % Across Familiar days
    
    
    
    sleeps = [3,5];
    for ep=1:length(sleeps)
        currep=sleeps(ep);
        eval(['fammeanfr_riprate(',num2str(ep),',:) = [Confammeanfr_riprateep',num2str(currep),';  Expfammeanfr_riprateep',num2str(currep),'];'])
    end
    
    figure; hold on; redimscreen_figforppt1;
    bar([2 3], fammeanfr_riprate,'grouped');
    errorbar(1.8,Confammeanfr_riprateep3,Confamerrfr_riprateep3,'b');
    errorbar(2.2,Expfammeanfr_riprateep3,Expfamerrfr_riprateep3,'m');
    errorbar(2.8,Confammeanfr_riprateep5,Confamerrfr_riprateep5,'b');
    errorbar(3.2,Expfammeanfr_riprateep5,Expfamerrfr_riprateep5,'m');
    
    % ----------- Stats -------------
    % Comparing animals: each animal has one number: average across days
    
    % 1)
    % Con against Con for different epochs
    [p35Con_famripratefr,h35Con_famripratefr] = ranksum(nanmean(Confrriprate_ep3(:,fam),2),nanmean(Confrriprate_ep5(:,fam),2));
    if h35Con_famripratefr==1,
        mul = sign(Confammeanfr_riprateep5);
        plot(2.8, (Confammeanfr_riprateep5+Confamerrfr_riprateep5+1.2*mul*Confamerrfr_riprateep5), 'b*','MarkerSize',8);
    end
    
    
    % 2)
    % Exp against Exp for different epochs
    [p35Exp_famripratefr,h35Exp_famripratefr] = ranksum(nanmean(Expfrriprate_ep3(:,fam),2),nanmean(Expfrriprate_ep5(:,fam),2));
    if h35Exp_famripratefr==1,
        mul = sign(Expfammeanfr_riprateep5);
        plot(3.2, (Expfammeanfr_riprateep5+Expfamerrfr_riprateep5+1.2*mul*Expfamerrfr_riprateep5), 'r*','MarkerSize',8);
    end
    
    % 3)
    % Con against Exp for each epoch
    [p33_famripratefr,h33_famripratefr] = ranksum(nanmean(Confrriprate_ep3(:,fam),2),nanmean(Expfrriprate_ep3(:,fam),2));
    [p55_famripratefr,h55_famripratefr] = ranksum(nanmean(Confrriprate_ep5(:,fam),2),nanmean(Expfrriprate_ep5(:,fam),2));
    
    % Plot midway between the two graphs
    
    if h33_famripratefr==1,
        mul = sign(Confammeanfr_riprateep3);
        plot(2, (Confammeanfr_riprateep3+Confamerrfr_riprateep3+1.2*mul*Confamerrfr_riprateep3), 'g*','MarkerSize',12);
    end
    
    if h55_famripratefr==1,
        mul = sign(Confammeanfr_riprateep5);
        plot(3, (Confammeanfr_riprateep5+Confamerrfr_riprateep5+1.2*mul*Confamerrfr_riprateep5), 'g*','MarkerSize',12);
    end
    
    % -------------------------
    
    title('Fraction Change in Ripple Rate during Sleep: Fam Days');
    ylabel('Fraction Change');
    set(gca,'XTick',[2:3],'XTickLabel',{'Post1';'Post2'},'FontSize',16,'Fontweight','normal');
    %set(gca,'XLim',[0.5 3.5]);
    %axis([0 4 0 0.5])
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    % ********** Repeat Mean Plots for StillTime   *****************
    
    % ********************** Still Time ****************************
    % ***************************************************************
    
    
    %-----------------------------------------------------------
    %-----------------------------------------------------------
    
    % 1) Bar plot for all, novel, familiar days, - Still Time
    
    set(0,'defaultaxesfontsize',20);set(0,'defaultaxesfontweight','normal');
    set(0,'defaultaxeslinewidth',2);
    
    % 1) Bar plot for all, novel, familiar days
    
    % Across all days
    
    sleeps = [1,3,5];
    for ep=1:length(sleeps)
        currep=sleeps(ep);
        eval(['allmean_stilltime(',num2str(ep),',:) = [Conallmean_stilltimeep',num2str(currep),';  Expallmean_stilltimeep',num2str(currep),'];'])
    end
    
    figure; hold on; redimscreen_figforppt1;
    bar(allmean_stilltime,'grouped');
    errorbar(0.8,Conallmean_stilltimeep1,Conallerr_stilltimeep1,'b');
    errorbar(1.2,Expallmean_stilltimeep1,Expallerr_stilltimeep1,'m');
    errorbar(1.8,Conallmean_stilltimeep3,Conallerr_stilltimeep3,'b');
    errorbar(2.2,Expallmean_stilltimeep3,Expallerr_stilltimeep3,'m');
    errorbar(2.8,Conallmean_stilltimeep5,Conallerr_stilltimeep5,'b');
    errorbar(3.2,Expallmean_stilltimeep5,Expallerr_stilltimeep5,'m');
    
    % ----------- Stats -------------
    % Comparing animals: each animal has one number: average across days
    
    % 1)
    % Con against Con for different epochs
    [p13Con_allstilltime,h13Con_allstilltime] = ranksum(nanmean(Constilltime_ep1,2),nanmean(Constilltime_ep3,2));
    [p15Con_allstilltime,h15Con_allstilltime] = ranksum(nanmean(Constilltime_ep1,2),nanmean(Constilltime_ep5,2));
    if h13Con_allstilltime==1,
        mul = sign(Conallmean_stilltimeep3);
        plot(1.8, (Conallmean_stilltimeep3+Conallerr_stilltimeep3+1.2*mul*Conallerr_stilltimeep3), 'b*','MarkerSize',8);
    end
    if h15Con_allstilltime==1,
        mul = sign(Conallmean_stilltimeep5);
        plot(2.8, (Conallmean_stilltimeep5+Conallerr_stilltimeep5+1.2*mul*Conallerr_stilltimeep5), 'b*','MarkerSize',8);
    end
    
    % 2)
    % Exp against Exp for different epochs
    [p13Exp_allstilltime,h13Exp_allstilltime] = ranksum(nanmean(Expstilltime_ep1,2),nanmean(Expstilltime_ep3,2));
    [p15Exp_allstilltime,h15Exp_allstilltime] = ranksum(nanmean(Expstilltime_ep1,2),nanmean(Expstilltime_ep5,2));
    if h13Exp_allstilltime==1,
        mul = sign(Expallmean_stilltimeep3);
        plot(2.2, (Expallmean_stilltimeep3+Expallerr_stilltimeep3+1.2*mul*Expallerr_stilltimeep3), 'r*','MarkerSize',8);
    end
    if h15Exp_allstilltime==1,
        mul = sign(Expallmean_stilltimeep5);
        plot(3.2, (Expallmean_stilltimeep5+Expallerr_stilltimeep5+1.2*mul*Expallerr_stilltimeep5), 'r*','MarkerSize',8);
    end
    
    % 3)
    % Con against Exp for each epoch
    [p11_allstilltime,h11_allstilltime] = ranksum(nanmean(Constilltime_ep1,2),nanmean(Expstilltime_ep1,2));
    [p33_allstilltime,h33_allstilltime] = ranksum(nanmean(Constilltime_ep3,2),nanmean(Expstilltime_ep3,2));
    [p55_allstilltime,h55_allstilltime] = ranksum(nanmean(Constilltime_ep5,2),nanmean(Expstilltime_ep5,2));
    
    % Plot midway between the two graphs
    if h11_allstilltime==1,
        mul = sign(Conallmean_stilltimeep1);
        plot(1, (Conallmean_stilltimeep1+Conallerr_stilltimeep1+1.2*mul*Conallerr_stilltimeep1), 'g*','MarkerSize',12);
    end
    
    if h33_allstilltime==1,
        mul = sign(Conallmean_stilltimeep3);
        plot(2, (Conallmean_stilltimeep3+Conallerr_stilltimeep3+1.2*mul*Conallerr_stilltimeep3), 'g*','MarkerSize',12);
    end
    
    if h55_allstilltime==1,
        mul = sign(Conallmean_stilltimeep5);
        plot(3, (Conallmean_stilltimeep5+Conallerr_stilltimeep5+1.2*mul*Conallerr_stilltimeep5), 'g*','MarkerSize',12);
    end
    
    % -------------------------
    
    title('Still Time during Sleep: All Days');
    ylabel('Still Time (s)');
    set(gca,'XTick',[1:3],'XTickLabel',{'Pre';'Post1';'Post2'},'FontSize',16,'Fontweight','normal');
    %set(gca,'XLim',[0.5 3.5]);
    %axis([0 4 0 0.5])
    
    
    %-----------------------------------------------------------
    %-----------------------------------------------------------
    
    % Across Novel days
    
    if ~isempty(intersect(days,novel))
        noveldays = intersect(days,novel);
        
        sleeps = [1,3,5];
        for ep=1:length(sleeps)
            currep=sleeps(ep);
            eval(['novmean_stilltime(',num2str(ep),',:) = [Connovmean_stilltimeep',num2str(currep),';  Expnovmean_stilltimeep',num2str(currep),'];'])
        end
        
        figure; hold on; redimscreen_figforppt1;
        bar(novmean_stilltime,'grouped');
        errorbar(0.8,Connovmean_stilltimeep1,Connoverr_stilltimeep1,'b');
        errorbar(1.2,Expnovmean_stilltimeep1,Expnoverr_stilltimeep1,'m');
        errorbar(1.8,Connovmean_stilltimeep3,Connoverr_stilltimeep3,'b');
        errorbar(2.2,Expnovmean_stilltimeep3,Expnoverr_stilltimeep3,'m');
        errorbar(2.8,Connovmean_stilltimeep5,Connoverr_stilltimeep5,'b');
        errorbar(3.2,Expnovmean_stilltimeep5,Expnoverr_stilltimeep5,'m');
        
        
        % ----------- Stats -------------
        % Comparing animals: each animal has one number: average across days
        
        % 1)
        % Con against Con for different epochs
        
        [p13Con_novstilltime,h13Con_novstilltime] = ranksum(nanmean(Constilltime_ep1(:,nov),2),nanmean(Constilltime_ep3(:,nov),2));
        [p15Con_novstilltime,h15Con_novstilltime] = ranksum(nanmean(Constilltime_ep1(:,nov),2),nanmean(Constilltime_ep5(:,nov),2));
        if h13Con_novstilltime==1,
            mul = sign(Connovmean_stilltimeep3);
            plot(1.8, (Connovmean_stilltimeep3+Connoverr_stilltimeep3+1.2*mul*Connoverr_stilltimeep3), 'b*','MarkerSize',8);
        end
        if h15Con_novstilltime==1,
            mul = sign(Connovmean_stilltimeep5);
            plot(2.8, (Connovmean_stilltimeep5+Connoverr_stilltimeep5+1.2*mul*Connoverr_stilltimeep5), 'b*','MarkerSize',8);
        end
        
        % 2)
        % Exp against Exp for different epochs
        [p13Exp_novstilltime,h13Exp_novstilltime] = ranksum(nanmean(Expstilltime_ep1(:,nov),2),nanmean(Expstilltime_ep3(:,nov),2));
        [p15Exp_novstilltime,h15Exp_novstilltime] = ranksum(nanmean(Expstilltime_ep1(:,nov),2),nanmean(Expstilltime_ep5(:,nov),2));
        if h13Exp_novstilltime==1,
            mul = sign(Expnovmean_stilltimeep3);
            plot(2.2, (Expnovmean_stilltimeep3+Expnoverr_stilltimeep3+1.2*mul*Expnoverr_stilltimeep3), 'r*','MarkerSize',8);
        end
        if h15Exp_novstilltime==1,
            mul = sign(Expnovmean_stilltimeep5);
            plot(3.2, (Expnovmean_stilltimeep5+Expnoverr_stilltimeep5+1.2*mul*Expnoverr_stilltimeep5), 'r*','MarkerSize',8);
        end
        
        % 3)
        % Con against Exp for each epoch
        [p11_novstilltime,h11_novstilltime] = ranksum(nanmean(Constilltime_ep1(:,nov),2),nanmean(Expstilltime_ep1(:,nov),2));
        [p33_novstilltime,h33_novstilltime] = ranksum(nanmean(Constilltime_ep3(:,nov),2),nanmean(Expstilltime_ep3(:,nov),2));
        [p55_novstilltime,h55_novstilltime] = ranksum(nanmean(Constilltime_ep5(:,nov),2),nanmean(Expstilltime_ep5(:,nov),2));
        
        % Plot midway between the two graphs
        if h11_novstilltime==1,
            mul = sign(Connovmean_stilltimeep1);
            plot(1, (Connovmean_stilltimeep1+Connoverr_stilltimeep1+1.2*mul*Connoverr_stilltimeep1), 'g*','MarkerSize',12);
        end
        
        if h33_novstilltime==1,
            mul = sign(Connovmean_stilltimeep3);
            plot(2, (Connovmean_stilltimeep3+Connoverr_stilltimeep3+1.2*mul*Connoverr_stilltimeep3), 'g*','MarkerSize',12);
        end
        
        if h55_novstilltime==1,
            mul = sign(Connovmean_stilltimeep5);
            plot(3, (Connovmean_stilltimeep5+Connoverr_stilltimeep5+1.2*mul*Connoverr_stilltimeep5), 'g*','MarkerSize',12);
        end
        
        % -------------------------
        
        
        title('Still Time during Sleep: Novel Days');
        ylabel('Still Time (s)');
        set(gca,'XTick',[1:3],'XTickLabel',{'Pre';'Post1';'Post2'},'FontSize',16,'Fontweight','normal');
    end
    
    
    
    %-----------------------------------------------------------
    %-----------------------------------------------------------
    
    % Across Familiar days
    
    if ~isempty(intersect(days,fam))
        famdays = intersect(days,fam);
        
        sleeps = [1,3,5];
        for ep=1:length(sleeps)
            currep=sleeps(ep);
            eval(['fammean_stilltime(',num2str(ep),',:) = [Confammean_stilltimeep',num2str(currep),';  Expfammean_stilltimeep',num2str(currep),'];'])
        end
        
        figure; hold on; redimscreen_figforppt1;
        bar(fammean_stilltime,'grouped');
        errorbar(0.8,Confammean_stilltimeep1,Confamerr_stilltimeep1,'b');
        errorbar(1.2,Expfammean_stilltimeep1,Expfamerr_stilltimeep1,'m');
        errorbar(1.8,Confammean_stilltimeep3,Confamerr_stilltimeep3,'b');
        errorbar(2.2,Expfammean_stilltimeep3,Expfamerr_stilltimeep3,'m');
        errorbar(2.8,Confammean_stilltimeep5,Confamerr_stilltimeep5,'b');
        errorbar(3.2,Expfammean_stilltimeep5,Expfamerr_stilltimeep5,'m');
        
        % ----------- Stats -------------
        % Comparing animals: each animal has one number: average across days
        
        % 1)
        % Con against Con for different epochs
        
        [p13Con_famstilltime,h13Con_famstilltime] = ranksum(nanmean(Constilltime_ep1(:,fam),2),nanmean(Constilltime_ep3(:,fam),2));
        [p15Con_famstilltime,h15Con_famstilltime] = ranksum(nanmean(Constilltime_ep1(:,fam),2),nanmean(Constilltime_ep5(:,fam),2));
        if h13Con_famstilltime==1,
            mul = sign(Confammean_stilltimeep3);
            plot(1.8, (Confammean_stilltimeep3+Confamerr_stilltimeep3+1.2*mul*Confamerr_stilltimeep3), 'b*','MarkerSize',8);
        end
        if h15Con_famstilltime==1,
            mul = sign(Confammean_stilltimeep5);
            plot(2.8, (Confammean_stilltimeep5+Confamerr_stilltimeep5+1.2*mul*Confamerr_stilltimeep5), 'b*','MarkerSize',8);
        end
        
        % 2)
        % Exp against Exp for different epochs
        [p13Exp_famstilltime,h13Exp_famstilltime] = ranksum(nanmean(Expstilltime_ep1(:,fam),2),nanmean(Expstilltime_ep3(:,fam),2));
        [p15Exp_famstilltime,h15Exp_famstilltime] = ranksum(nanmean(Expstilltime_ep1(:,fam),2),nanmean(Expstilltime_ep5(:,fam),2));
        if h13Exp_famstilltime==1,
            mul = sign(Expfammean_stilltimeep3);
            plot(2.2, (Expfammean_stilltimeep3+Expfamerr_stilltimeep3+1.2*mul*Expfamerr_stilltimeep3), 'r*','MarkerSize',8);
        end
        if h15Exp_famstilltime==1,
            mul = sign(Expfammean_stilltimeep5);
            plot(3.2, (Expfammean_stilltimeep5+Expfamerr_stilltimeep5+1.2*mul*Expfamerr_stilltimeep5), 'r*','MarkerSize',8);
        end
        
        % 3)
        % Con against Exp for each epoch
        [p11_famstilltime,h11_famstilltime] = ranksum(nanmean(Constilltime_ep1(:,fam),2),nanmean(Expstilltime_ep1(:,fam),2));
        [p33_famstilltime,h33_famstilltime] = ranksum(nanmean(Constilltime_ep3(:,fam),2),nanmean(Expstilltime_ep3(:,fam),2));
        [p55_famstilltime,h55_famstilltime] = ranksum(nanmean(Constilltime_ep5(:,fam),2),nanmean(Expstilltime_ep5(:,fam),2));
        
        % Plot midway between the two graphs
        if h11_famstilltime==1,
            mul = sign(Confammean_stilltimeep1);
            plot(1, (Confammean_stilltimeep1+Confamerr_stilltimeep1+1.2*mul*Confamerr_stilltimeep1), 'g*','MarkerSize',12);
        end
        
        if h33_famstilltime==1,
            mul = sign(Confammean_stilltimeep3);
            plot(2, (Confammean_stilltimeep3+Confamerr_stilltimeep3+1.2*mul*Confamerr_stilltimeep3), 'g*','MarkerSize',12);
        end
        
        if h55_famstilltime==1,
            mul = sign(Confammean_stilltimeep5);
            plot(3, (Confammean_stilltimeep5+Confamerr_stilltimeep5+1.2*mul*Confamerr_stilltimeep5), 'g*','MarkerSize',12);
        end
        
        % -------------------------
        
        title('Still Time during Sleep: Familiar Days');
        ylabel('Still Time (s)');
        set(gca,'XTick',[1:3],'XTickLabel',{'Pre';'Post1';'Post2'},'FontSize',16,'Fontweight','normal');
    end
    
    
    
    
    
    
    
    
    
    
    % *********  Vs Days: Fraction Change In Ripple Rate ************
    % ***************************************************************
    
    
    %-----------------------------------------------------------
    %-----------------------------------------------------------
    
    % 1) Vs. Each Day
    
    figure; hold on; redimscreen_figforppt1;
    errorbar(Conmeanfr_riprateep3,Conerrfr_riprateep3,'b--','LineWidth',1);
    errorbar(Conmeanfr_riprateep5,Conerrfr_riprateep5,'b','LineWidth',2);
    
    errorbar(Expmeanfr_riprateep3,Experrfr_riprateep3,'r--','LineWidth',1);
    errorbar(Expmeanfr_riprateep5,Experrfr_riprateep5,'r','LineWidth',2);
    
    % Need to find difference from 1 using distribution around 1
    %     for d = 2:length(days)
    %          [p13_frriprate(d),h13_frriprate(d)] = ranksum(frriprate_ep1(:,d), frriprate_ep3(:,d) );
    %     end
    
    if ntet==2
        set(gca,'YLim',[0 4.0]);
    else
        set(gca,'YLim',[0 2.2]);
    end
    set(gca,'XLim',[0.5 max(days)+0.5]);
    title('Increase in Ripple Rate during Post-Beh Sleep');
    ylabel('Fraction Change in Ripple Rate'); xlabel('Days');
    
    %  Vs. Pairs of Sliding Days
    
    figure; hold on; redimscreen_figforppt1;
    errorbar(Conmeanfrslide_riprateep3,Conerrfrslide_riprateep3,'b--','LineWidth',1);
    errorbar(Conmeanfrslide_riprateep5,Conerrfrslide_riprateep5,'b','LineWidth',2);
    
    errorbar(Expmeanfrslide_riprateep3,Experrfrslide_riprateep3,'r--','LineWidth',1);
    errorbar(Expmeanfrslide_riprateep5,Experrfrslide_riprateep5,'r','LineWidth',2);
    
    if ntet==2
        set(gca,'YLim',[0 4.0]);
    else
        set(gca,'YLim',[0 2.2]);
    end
    set(gca,'XLim',[0.5 max(days)-0.5]);
    title('Increase in Ripple Rate during Post-Beh Sleep');
    ylabel('Fraction Change in Ripple Rate'); xlabel('Pairs of Sliding Days');
    
end % end figopt1






%     % Repeat All Plots for Percent Time during Ripples (Eqvt to Rip Rate*Rip Lth)
%
%     % ********************** Ripple Percent ****************************
%     % ***************************************************************
%
%
%
%     %-----------------------------------------------------------
%     %-----------------------------------------------------------
%
%     % 1) Bar plot for all, novel, familiar days, & 2) Across days  - Mean Ripple Per Time
%
%     set(0,'defaultaxesfontsize',20);set(0,'defaultaxesfontweight','normal');
%     set(0,'defaultaxeslinewidth',2);
%
%     % 1) Bar plot for all, novel, familiar days
%
%     % Across all days
%
%     sleeps = [1,3,5];
%     for ep=1:length(sleeps)
%         currep=sleeps(ep);
%         eval(['allmean_ripper(',num2str(ep),',:) = [Conallmean_ripperep',num2str(currep),';  Expallmean_ripperep',num2str(currep),'];'])
%     end
%
%     figure; hold on; redimscreen_figforppt1;
%     bar(allmean_ripper,'grouped');
%     errorbar(0.8,Conallmean_ripperep1,Conallerr_ripperep1,'b');
%     errorbar(1.2,Expallmean_ripperep1,Expallerr_ripperep1,'m');
%     errorbar(1.8,Conallmean_ripperep3,Conallerr_ripperep3,'b');
%     errorbar(2.2,Expallmean_ripperep3,Expallerr_ripperep3,'m');
%     errorbar(2.8,Conallmean_ripperep5,Conallerr_ripperep5,'b');
%     errorbar(3.2,Expallmean_ripperep5,Expallerr_ripperep5,'m');
%
%     % ----------- Stats -------------
%     % Comparing animals: each animal has one number: average across days
%
%     % 1)
%     % Con against Con for different epochs
%     [p13Con_allripper,h13Con_allripper] = ranksum(nanmean(Conripper_ep1,2),nanmean(Conripper_ep3,2));
%     [p15Con_allripper,h15Con_allripper] = ranksum(nanmean(Conripper_ep1,2),nanmean(Conripper_ep5,2));
%     if h13Con_allripper==1,
%         mul = sign(Conallmean_ripperep3);
%         plot(1.8, (Conallmean_ripperep3+Conallerr_ripperep3+1.2*mul*Conallerr_ripperep3), 'b*','MarkerSize',8);
%     end
%     if h15Con_allripper==1,
%         mul = sign(Conallmean_ripperep5);
%         plot(2.8, (Conallmean_ripperep5+Conallerr_ripperep5+1.2*mul*Conallerr_ripperep5), 'b*','MarkerSize',8);
%     end
%
%     % 2)
%     % Exp against Exp for different epochs
%     [p13Exp_allripper,h13Exp_allripper] = ranksum(nanmean(Expripper_ep1,2),nanmean(Expripper_ep3,2));
%     [p15Exp_allripper,h15Exp_allripper] = ranksum(nanmean(Expripper_ep1,2),nanmean(Expripper_ep5,2));
%     if h13Exp_allripper==1,
%         mul = sign(Expallmean_ripperep3);
%         plot(2.2, (Expallmean_ripperep3+Expallerr_ripperep3+1.2*mul*Expallerr_ripperep3), 'r*','MarkerSize',8);
%     end
%     if h15Exp_allripper==1,
%         mul = sign(Expallmean_ripperep5);
%         plot(3.2, (Expallmean_ripperep5+Expallerr_ripperep5+1.2*mul*Expallerr_ripperep5), 'r*','MarkerSize',8);
%     end
%
%     % 3)
%     % Con against Exp for each epoch
%     [p11_allripper,h11_allripper] = ranksum(nanmean(Conripper_ep1,2),nanmean(Expripper_ep1,2));
%     [p33_allripper,h33_allripper] = ranksum(nanmean(Conripper_ep3,2),nanmean(Expripper_ep3,2));
%     [p55_allripper,h55_allripper] = ranksum(nanmean(Conripper_ep5,2),nanmean(Expripper_ep5,2));
%
%     % Plot midway between the two graphs
%     if h11_allripper==1,
%         mul = sign(Conallmean_ripperep1);
%         plot(1, (Conallmean_ripperep1+Conallerr_ripperep1+1.2*mul*Conallerr_ripperep1), 'g*','MarkerSize',12);
%     end
%
%     if h33_allripper==1,
%         mul = sign(Conallmean_ripperep3);
%         plot(2, (Conallmean_ripperep3+Conallerr_ripperep3+1.2*mul*Conallerr_ripperep3), 'g*','MarkerSize',12);
%     end
%
%     if h55_allripper==1,
%         mul = sign(Conallmean_ripperep5);
%         plot(3, (Conallmean_ripperep5+Conallerr_ripperep5+1.2*mul*Conallerr_ripperep5), 'g*','MarkerSize',12);
%     end
%
%     % -------------------------
%
%     title('Percent Ripple Time during Sleep: All Days');
%     ylabel('Percent Ripple Time');
%     set(gca,'XTick',[1:3],'XTickLabel',{'Pre';'Post1';'Post2'},'FontSize',16,'Fontweight','normal');
%     %set(gca,'XLim',[0.5 3.5]);
%     %axis([0 4 0 0.5])
%
%
%     %-----------------------------------------------------------
%     %-----------------------------------------------------------
%
%     % Across Novel days
%
%     if ~isempty(intersect(days,novel))
%         noveldays = intersect(days,novel);
%
%         sleeps = [1,3,5];
%         for ep=1:length(sleeps)
%             currep=sleeps(ep);
%             eval(['novmean_ripper(',num2str(ep),',:) = [Connovmean_ripperep',num2str(currep),';  Expnovmean_ripperep',num2str(currep),'];'])
%         end
%
%         figure; hold on; redimscreen_figforppt1;
%         bar(novmean_ripper,'grouped');
%         errorbar(0.8,Connovmean_ripperep1,Connoverr_ripperep1,'b');
%         errorbar(1.2,Expnovmean_ripperep1,Expnoverr_ripperep1,'m');
%         errorbar(1.8,Connovmean_ripperep3,Connoverr_ripperep3,'b');
%         errorbar(2.2,Expnovmean_ripperep3,Expnoverr_ripperep3,'m');
%         errorbar(2.8,Connovmean_ripperep5,Connoverr_ripperep5,'b');
%         errorbar(3.2,Expnovmean_ripperep5,Expnoverr_ripperep5,'m');
%
%
%         % ----------- Stats -------------
%         % Comparing animals: each animal has one number: average across days
%
%         % 1)
%         % Con against Con for different epochs
%
%         [p13Con_novripper,h13Con_novripper] = ranksum(nanmean(Conripper_ep1(:,nov),2),nanmean(Conripper_ep3(:,nov),2));
%         [p15Con_novripper,h15Con_novripper] = ranksum(nanmean(Conripper_ep1(:,nov),2),nanmean(Conripper_ep5(:,nov),2));
%         if h13Con_novripper==1,
%             mul = sign(Connovmean_ripperep3);
%             plot(1.8, (Connovmean_ripperep3+Connoverr_ripperep3+1.2*mul*Connoverr_ripperep3), 'b*','MarkerSize',8);
%         end
%         if h15Con_novripper==1,
%             mul = sign(Connovmean_ripperep5);
%             plot(2.8, (Connovmean_ripperep5+Connoverr_ripperep5+1.2*mul*Connoverr_ripperep5), 'b*','MarkerSize',8);
%         end
%
%         % 2)
%         % Exp against Exp for different epochs
%         [p13Exp_novripper,h13Exp_novripper] = ranksum(nanmean(Expripper_ep1(:,nov),2),nanmean(Expripper_ep3(:,nov),2));
%         [p15Exp_novripper,h15Exp_novripper] = ranksum(nanmean(Expripper_ep1(:,nov),2),nanmean(Expripper_ep5(:,nov),2));
%         if h13Exp_novripper==1,
%             mul = sign(Expnovmean_ripperep3);
%             plot(2.2, (Expnovmean_ripperep3+Expnoverr_ripperep3+1.2*mul*Expnoverr_ripperep3), 'r*','MarkerSize',8);
%         end
%         if h15Exp_novripper==1,
%             mul = sign(Expnovmean_ripperep5);
%             plot(3.2, (Expnovmean_ripperep5+Expnoverr_ripperep5+1.2*mul*Expnoverr_ripperep5), 'r*','MarkerSize',8);
%         end
%
%         % 3)
%         % Con against Exp for each epoch
%         [p11_novripper,h11_novripper] = ranksum(nanmean(Conripper_ep1(:,nov),2),nanmean(Expripper_ep1(:,nov),2));
%         [p33_novripper,h33_novripper] = ranksum(nanmean(Conripper_ep3(:,nov),2),nanmean(Expripper_ep3(:,nov),2));
%         [p55_novripper,h55_novripper] = ranksum(nanmean(Conripper_ep5(:,nov),2),nanmean(Expripper_ep5(:,nov),2));
%
%         % Plot midway between the two graphs
%         if h11_novripper==1,
%             mul = sign(Connovmean_ripperep1);
%             plot(1, (Connovmean_ripperep1+Connoverr_ripperep1+1.2*mul*Connoverr_ripperep1), 'g*','MarkerSize',12);
%         end
%
%         if h33_novripper==1,
%             mul = sign(Connovmean_ripperep3);
%             plot(2, (Connovmean_ripperep3+Connoverr_ripperep3+1.2*mul*Connoverr_ripperep3), 'g*','MarkerSize',12);
%         end
%
%         if h55_novripper==1,
%             mul = sign(Connovmean_ripperep5);
%             plot(3, (Connovmean_ripperep5+Connoverr_ripperep5+1.2*mul*Connoverr_ripperep5), 'g*','MarkerSize',12);
%         end
%
%         % -------------------------
%
%
%         title('Percent Ripple Time during Sleep: Novel Days');
%         ylabel('Percent Ripple Time');
%         set(gca,'XTick',[1:3],'XTickLabel',{'Pre';'Post1';'Post2'},'FontSize',16,'Fontweight','normal');
%     end
%
%
%
%     %-----------------------------------------------------------
%     %-----------------------------------------------------------
%
%     % Across Familiar days
%
%     if ~isempty(intersect(days,fam))
%         famdays = intersect(days,fam);
%
%         sleeps = [1,3,5];
%         for ep=1:length(sleeps)
%             currep=sleeps(ep);
%             eval(['fammean_ripper(',num2str(ep),',:) = [Confammean_ripperep',num2str(currep),';  Expfammean_ripperep',num2str(currep),'];'])
%         end
%
%         figure; hold on; redimscreen_figforppt1;
%         bar(fammean_ripper,'grouped');
%         errorbar(0.8,Confammean_ripperep1,Confamerr_ripperep1,'b');
%         errorbar(1.2,Expfammean_ripperep1,Expfamerr_ripperep1,'m');
%         errorbar(1.8,Confammean_ripperep3,Confamerr_ripperep3,'b');
%         errorbar(2.2,Expfammean_ripperep3,Expfamerr_ripperep3,'m');
%         errorbar(2.8,Confammean_ripperep5,Confamerr_ripperep5,'b');
%         errorbar(3.2,Expfammean_ripperep5,Expfamerr_ripperep5,'m');
%
%         % ----------- Stats -------------
%         % Comparing animals: each animal has one number: average across days
%
%         % 1)
%         % Con against Con for different epochs
%
%         [p13Con_famripper,h13Con_famripper] = ranksum(nanmean(Conripper_ep1(:,fam),2),nanmean(Conripper_ep3(:,fam),2));
%         [p15Con_famripper,h15Con_famripper] = ranksum(nanmean(Conripper_ep1(:,fam),2),nanmean(Conripper_ep5(:,fam),2));
%         if h13Con_famripper==1,
%             mul = sign(Confammean_ripperep3);
%             plot(1.8, (Confammean_ripperep3+Confamerr_ripperep3+1.2*mul*Confamerr_ripperep3), 'b*','MarkerSize',8);
%         end
%         if h15Con_famripper==1,
%             mul = sign(Confammean_ripperep5);
%             plot(2.8, (Confammean_ripperep5+Confamerr_ripperep5+1.2*mul*Confamerr_ripperep5), 'b*','MarkerSize',8);
%         end
%
%         % 2)
%         % Exp against Exp for different epochs
%         [p13Exp_famripper,h13Exp_famripper] = ranksum(nanmean(Expripper_ep1(:,fam),2),nanmean(Expripper_ep3(:,fam),2));
%         [p15Exp_famripper,h15Exp_famripper] = ranksum(nanmean(Expripper_ep1(:,fam),2),nanmean(Expripper_ep5(:,fam),2));
%         if h13Exp_famripper==1,
%             mul = sign(Expfammean_ripperep3);
%             plot(2.2, (Expfammean_ripperep3+Expfamerr_ripperep3+1.2*mul*Expfamerr_ripperep3), 'r*','MarkerSize',8);
%         end
%         if h15Exp_famripper==1,
%             mul = sign(Expfammean_ripperep5);
%             plot(3.2, (Expfammean_ripperep5+Expfamerr_ripperep5+1.2*mul*Expfamerr_ripperep5), 'r*','MarkerSize',8);
%         end
%
%         % 3)
%         % Con against Exp for each epoch
%         [p11_famripper,h11_famripper] = ranksum(nanmean(Conripper_ep1(:,fam),2),nanmean(Expripper_ep1(:,fam),2));
%         [p33_famripper,h33_famripper] = ranksum(nanmean(Conripper_ep3(:,fam),2),nanmean(Expripper_ep3(:,fam),2));
%         [p55_famripper,h55_famripper] = ranksum(nanmean(Conripper_ep5(:,fam),2),nanmean(Expripper_ep5(:,fam),2));
%
%         % Plot midway between the two graphs
%         if h11_famripper==1,
%             mul = sign(Confammean_ripperep1);
%             plot(1, (Confammean_ripperep1+Confamerr_ripperep1+1.2*mul*Confamerr_ripperep1), 'g*','MarkerSize',12);
%         end
%
%         if h33_famripper==1,
%             mul = sign(Confammean_ripperep3);
%             plot(2, (Confammean_ripperep3+Confamerr_ripperep3+1.2*mul*Confamerr_ripperep3), 'g*','MarkerSize',12);
%         end
%
%         if h55_famripper==1,
%             mul = sign(Confammean_ripperep5);
%             plot(3, (Confammean_ripperep5+Confamerr_ripperep5+1.2*mul*Confamerr_ripperep5), 'g*','MarkerSize',12);
%         end
%
%         % -------------------------
%
%         title('Percent Ripple Time during Sleep: Familiar Days');
%         ylabel('Percent Ripple Time (Hz)');
%         set(gca,'XTick',[1:3],'XTickLabel',{'Pre';'Post1';'Post2'},'FontSize',16,'Fontweight','normal');
%     end
%
%
%
%
%
%
%
%     % ********** Fraction Change In Percent Ripple Time *************
%     % ***************************************************************
%
%     %-----------------------------------------------------------
%     %-----------------------------------------------------------
%
%     % 1) Bar plot for all, novel, familiar days, & 2) Across days
%
%     set(0,'defaultaxesfontsize',20);set(0,'defaultaxesfontweight','normal');
%     set(0,'defaultaxeslinewidth',2);
%
%     % 1) Bar plot for all, novel, familiar days
%
%     % Across all days
%
%     sleeps = [3,5];
%     for ep=1:length(sleeps)
%         currep=sleeps(ep);
%         eval(['allmeanfr_ripper(',num2str(ep),',:) = [Conallmeanfr_ripperep',num2str(currep),';  Expallmeanfr_ripperep',num2str(currep),'];'])
%     end
%
%     figure; hold on; redimscreen_figforppt1;
%     bar([2 3], allmeanfr_ripper,'grouped');
%     errorbar(1.8,Conallmeanfr_ripperep3,Conallerrfr_ripperep3,'b');
%     errorbar(2.2,Expallmeanfr_ripperep3,Expallerrfr_ripperep3,'m');
%     errorbar(2.8,Conallmeanfr_ripperep5,Conallerrfr_ripperep5,'b');
%     errorbar(3.2,Expallmeanfr_ripperep5,Expallerrfr_ripperep5,'m');
%
%     % ----------- Stats -------------
%     % Comparing animals: each animal has one number: average across days
%
%     % 1)
%     % Con against Con for different epochs
%     [p35Con_allripperfr,h35Con_allripperfr] = ranksum(nanmean(Confrripper_ep3,2),nanmean(Confrripper_ep5,2));
%     if h35Con_allripperfr==1,
%         mul = sign(Conallmeanfr_ripperep5);
%         plot(2.8, (Conallmeanfr_ripperep5+Conallerrfr_ripperep5+1.2*mul*Conallerrfr_ripperep5), 'b*','MarkerSize',8);
%     end
%
%
%     % 2)
%     % Exp against Exp for different epochs
%     [p35Exp_allripperfr,h35Exp_allripperfr] = ranksum(nanmean(Expfrripper_ep3,2),nanmean(Expfrripper_ep5,2));
%     if h35Exp_allripperfr==1,
%         mul = sign(Expallmeanfr_ripperep5);
%         plot(3.2, (Expallmeanfr_ripperep5+Expallerrfr_ripperep5+1.2*mul*Expallerrfr_ripperep5), 'r*','MarkerSize',8);
%     end
%
%     % 3)
%     % Con against Exp for each epoch
%     [p33_allripperfr,h33_allripperfr] = ranksum(nanmean(Confrripper_ep3,2),nanmean(Expfrripper_ep3,2));
%     [p55_allripperfr,h55_allripperfr] = ranksum(nanmean(Confrripper_ep5,2),nanmean(Expfrripper_ep5,2));
%
%     % Plot midway between the two graphs
%
%     if h33_allripperfr==1,
%         mul = sign(Conallmeanfr_ripperep3);
%         plot(2, (Conallmeanfr_ripperep3+Conallerrfr_ripperep3+1.2*mul*Conallerrfr_ripperep3), 'g*','MarkerSize',12);
%     end
%
%     if h55_allripperfr==1,
%         mul = sign(Conallmeanfr_ripperep5);
%         plot(3, (Conallmeanfr_ripperep5+Conallerrfr_ripperep5+1.2*mul*Conallerrfr_ripperep5), 'g*','MarkerSize',12);
%     end
%
%     % -------------------------
%
%     title('Fraction Change in Percent Ripple Time during Sleep: All Days');
%     ylabel('Fraction Change');
%     set(gca,'XTick',[2:3],'XTickLabel',{'Post1';'Post2'},'FontSize',16,'Fontweight','normal');
%     %set(gca,'XLim',[0.5 3.5]);
%     %axis([0 4 0 0.5])
%
%
%
%     %-----------------------------------------------------------
%     %-----------------------------------------------------------
%
%     % Across Novel days
%
%     sleeps = [3,5];
%     for ep=1:length(sleeps)
%         currep=sleeps(ep);
%         eval(['novmeanfr_ripper(',num2str(ep),',:) = [Connovmeanfr_ripperep',num2str(currep),';  Expnovmeanfr_ripperep',num2str(currep),'];'])
%     end
%
%     figure; hold on; redimscreen_figforppt1;
%     bar([2 3], novmeanfr_ripper,'grouped');
%     errorbar(1.8,Connovmeanfr_ripperep3,Connoverrfr_ripperep3,'b');
%     errorbar(2.2,Expnovmeanfr_ripperep3,Expnoverrfr_ripperep3,'m');
%     errorbar(2.8,Connovmeanfr_ripperep5,Connoverrfr_ripperep5,'b');
%     errorbar(3.2,Expnovmeanfr_ripperep5,Expnoverrfr_ripperep5,'m');
%
%     % ----------- Stats -------------
%     % Comparing animals: each animal has one number: average across days
%
%     % 1)
%     % Con against Con for different epochs
%     [p35Con_novripperfr,h35Con_novripperfr] = ranksum(nanmean(Confrripper_ep3(:,nov),2),nanmean(Confrripper_ep5(:,nov),2));
%     if h35Con_novripperfr==1,
%         mul = sign(Connovmeanfr_ripperep5);
%         plot(2.8, (Connovmeanfr_ripperep5+Connoverrfr_ripperep5+1.2*mul*Connoverrfr_ripperep5), 'b*','MarkerSize',8);
%     end
%
%
%     % 2)
%     % Exp against Exp for different epochs
%     [p35Exp_novripperfr,h35Exp_novripperfr] = ranksum(nanmean(Expfrripper_ep3(:,nov),2),nanmean(Expfrripper_ep5(:,nov),2));
%     if h35Exp_novripperfr==1,
%         mul = sign(Expnovmeanfr_ripperep5);
%         plot(3.2, (Expnovmeanfr_ripperep5+Expnoverrfr_ripperep5+1.2*mul*Expnoverrfr_ripperep5), 'r*','MarkerSize',8);
%     end
%
%     % 3)
%     % Con against Exp for each epoch
%     [p33_novripperfr,h33_novripperfr] = ranksum(nanmean(Confrripper_ep3(:,nov),2),nanmean(Expfrripper_ep3(:,nov),2));
%     [p55_novripperfr,h55_novripperfr] = ranksum(nanmean(Confrripper_ep5(:,nov),2),nanmean(Expfrripper_ep5(:,nov),2));
%
%     % Plot midway between the two graphs
%
%     if h33_novripperfr==1,
%         mul = sign(Connovmeanfr_ripperep3);
%         plot(2, (Connovmeanfr_ripperep3+Connoverrfr_ripperep3+1.2*mul*Connoverrfr_ripperep3), 'g*','MarkerSize',12);
%     end
%
%     if h55_novripperfr==1,
%         mul = sign(Connovmeanfr_ripperep5);
%         plot(3, (Connovmeanfr_ripperep5+Connoverrfr_ripperep5+1.2*mul*Connoverrfr_ripperep5), 'g*','MarkerSize',12);
%     end
%
%     % -------------------------
%
%     title('Fraction Change in Percent Ripple Time during Sleep: Nov Days');
%     ylabel('Fraction Change');
%     set(gca,'XTick',[2:3],'XTickLabel',{'Post1';'Post2'},'FontSize',16,'Fontweight','normal');
%     %set(gca,'XLim',[0.5 3.5]);
%     %axis([0 4 0 0.5])
%
%
%
%
%
%
%     %-----------------------------------------------------------
%     %-----------------------------------------------------------
%
%     % Across Familiar days
%
%
%
%     sleeps = [3,5];
%     for ep=1:length(sleeps)
%         currep=sleeps(ep);
%         eval(['fammeanfr_ripper(',num2str(ep),',:) = [Confammeanfr_ripperep',num2str(currep),';  Expfammeanfr_ripperep',num2str(currep),'];'])
%     end
%
%     figure; hold on; redimscreen_figforppt1;
%     bar([2 3], fammeanfr_ripper,'grouped');
%     errorbar(1.8,Confammeanfr_ripperep3,Confamerrfr_ripperep3,'b');
%     errorbar(2.2,Expfammeanfr_ripperep3,Expfamerrfr_ripperep3,'m');
%     errorbar(2.8,Confammeanfr_ripperep5,Confamerrfr_ripperep5,'b');
%     errorbar(3.2,Expfammeanfr_ripperep5,Expfamerrfr_ripperep5,'m');
%
%     % ----------- Stats -------------
%     % Comparing animals: each animal has one number: average across days
%
%     % 1)
%     % Con against Con for different epochs
%     [p35Con_famripperfr,h35Con_famripperfr] = ranksum(nanmean(Confrripper_ep3(:,fam),2),nanmean(Confrripper_ep5(:,fam),2));
%     if h35Con_famripperfr==1,
%         mul = sign(Confammeanfr_ripperep5);
%         plot(2.8, (Confammeanfr_ripperep5+Confamerrfr_ripperep5+1.2*mul*Confamerrfr_ripperep5), 'b*','MarkerSize',8);
%     end
%
%
%     % 2)
%     % Exp against Exp for different epochs
%     [p35Exp_famripperfr,h35Exp_famripperfr] = ranksum(nanmean(Expfrripper_ep3(:,fam),2),nanmean(Expfrripper_ep5(:,fam),2));
%     if h35Exp_famripperfr==1,
%         mul = sign(Expfammeanfr_ripperep5);
%         plot(3.2, (Expfammeanfr_ripperep5+Expfamerrfr_ripperep5+1.2*mul*Expfamerrfr_ripperep5), 'r*','MarkerSize',8);
%     end
%
%     % 3)
%     % Con against Exp for each epoch
%     [p33_famripperfr,h33_famripperfr] = ranksum(nanmean(Confrripper_ep3(:,fam),2),nanmean(Expfrripper_ep3(:,fam),2));
%     [p55_famripperfr,h55_famripperfr] = ranksum(nanmean(Confrripper_ep5(:,fam),2),nanmean(Expfrripper_ep5(:,fam),2));
%
%     % Plot midway between the two graphs
%
%     if h33_famripperfr==1,
%         mul = sign(Confammeanfr_ripperep3);
%         plot(2, (Confammeanfr_ripperep3+Confamerrfr_ripperep3+1.2*mul*Confamerrfr_ripperep3), 'g*','MarkerSize',12);
%     end
%
%     if h55_famripperfr==1,
%         mul = sign(Confammeanfr_ripperep5);
%         plot(3, (Confammeanfr_ripperep5+Confamerrfr_ripperep5+1.2*mul*Confamerrfr_ripperep5), 'g*','MarkerSize',12);
%     end
%
%     % -------------------------
%
%     title('Fraction Change in Percent Ripple Time during Sleep: Fam Days');
%     ylabel('Fraction Change');
%     set(gca,'XTick',[2:3],'XTickLabel',{'Post1';'Post2'},'FontSize',16,'Fontweight','normal');
%     %set(gca,'XLim',[0.5 3.5]);
%     %axis([0 4 0 0.5])




















