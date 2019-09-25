
% Changing plotting from DFSsj_getriprate by making it similar to DFSsj_getriprate_run
% Also, can use ripplesep1 now instead of ripples - referenced to epoch1
% Getting rid of novel and familiar plots. If you want them again, see DFSsj_getriprate
% Also, changing ripple size calc and plots to match DFSsj_getriprate_run

% Plotting of both groups together
% DIO is controlled for

clear; close all;
runscript = 0;
savedata = 0; % save data option - only works if runscript is also on
figopt1=1; % Figure Options for Ripple Rate

figopt2=0; % Figure Options for Ripple Size
% Sub-options under figopt2
figoptsize = 1; % For plotting ripple-size plots
figoptntet = 0;
dodayripsize=0; % Per day- turned off by default
ntet=1;

savedir = '/data25/sjadhav/RippleInterruption/ProcessedData/';

% First set - Use DFAsj_getriprate_noDIO
% savefile1 = [savedir 'RippleRate_ExpGrp_nostim'];
% savefile2 = [savedir 'RippleRate_ConGrp_nostim'];
% savefile = [savedir 'RippleRate_All_nostim'];

% Second Set - Use DFTF_getstimtimes to filter out stim times, and then use
% savefile1 = [savedir 'RippleRate_ExpGrp_filtstim'];
% savefile2 = [savedir 'RippleRate_ConGrp_filtstim'];
% savefile = [savedir 'RippleRate_All_filtstim'];

% Revision 2 - More animals in Expt Grp
savefile1 = [savedir 'RippleRate_ExpGrp_filtstim_rev2'];
savefile2 = [savedir 'RippleRate_ConGrp_filtstim_rev2'];
savefile = [savedir 'RippleRate_All_filtstim_rev2'];


% If runscript, run Datafilter and save data
if runscript == 1
    
    %Animal selection
    %-----------------------------------------------------
    Expanimals = {'REc','REd','REe','REf','REg','REh'};
    %Expanimals = {'REg'};
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

figdir = '/data25/sjadhav/RippleInterruption/Figures/01AugSep11_RippleDisFigs/RippleRate/Rest/';
summdir = figdir;
set(0,'defaultaxesfontweight','normal'); set(0,'defaultaxeslinewidth',2);

if forppr==1
    set(0,'defaultaxesfontsize',16);
    tfont = 14; % title font
    xfont = 14;
    yfont = 14;
else
    set(0,'defaultaxesfontsize',24);
    tfont = 28;
    xfont = 20;
    yfont = 20;
end

clr = {'b','r','g','c','m','y','k','r'};
savefig1=0;
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
                eval([str(g,:),'ripbaseline_ep',num2str(currep),'{an}{d} =',str(g,:),'ripf(an).output{1}(index).rip_baseline;']); % Single number
                eval([str(g,:),'ripstd_ep',num2str(currep),'{an}{d} =',str(g,:),'ripf(an).output{1}(index).rip_std;']); % Single number
                eval([str(g,:),'ripthresh_ep',num2str(currep),'{an}{d} =',str(g,:),'ripf(an).output{1}(index).rip_thresh;']); % Single number
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
        
        % For Ripple Percent
        %-------------------
        
        eval([str(g,:),'frripper_ep',num2str(ep),' =',str(g,:),'ripper_ep',num2str(ep),'./',str(g,:),'ripper_ep1;']); % Fraction change: Matrix divide
        
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
        
        % Get means for each day
        for d=1:length(days)
            currday=days(d);
            eval([str(g,:),'daymean_ripperep',num2str(ep),'(',num2str(currday),')= nanmean(',str(g,:),'ripper_ep',num2str(ep),'(:,d));']);
            eval([str(g,:),'dayerr_ripperep',num2str(ep),'(',num2str(currday),')= nansem(',str(g,:),'ripper_ep',num2str(ep),'(:,d));']);
        end
        
        % For Still Time
        %-------------------
        
        % Avg across all days
        eval([str(g,:),'allmean_stilltimeep' num2str(ep) '= nanmean(nanmean(',str(g,:),'stilltime_ep',num2str(ep),'));']) % Across animals, then across days
        %eval([str(g,:),'allerr_stilltimeep' num2str(ep) '= nansem(nanmean(',str(g,:),'stilltime_ep',num2str(ep),',2));']) % Across days - mean for each animal - then sem of that
        % Take sem across days and animals: each measure is independent
        eval([str(g,:),'allerr_stilltimeep' num2str(ep) '= nanstd(',str(g,:),'stilltime_ep',num2str(ep),'(:));']) % Across days and animals - sem of whole thing
        
        % Avg across novel days
        %eval([str(g,:),'novmean_stilltimeep' num2str(ep) '= nanmean(nanmean(',str(g,:),'stilltime_ep',num2str(ep),'(:,nov)));'])
        %eval([str(g,:),'noverr_stilltimeep' num2str(ep) '= nansem(nanmean(',str(g,:),'stilltime_ep',num2str(ep),'(:,nov),2));']) % Across days - mean for each animal - then sem of that
        % Change - like above, sem of all data together
        % Avg across fam days
        %eval([str(g,:),'fammean_stilltimeep' num2str(ep) '= nanmean(nanmean(',str(g,:),'stilltime_ep',num2str(ep),'(:,fam)));'])
        %eval([str(g,:),'famerr_stilltimeep' num2str(ep) '= nansem(nanmean(',str(g,:),'stilltime_ep',num2str(ep),'(:,fam),2));']) % Across days - mean for each animal - then sem of that
        % Change - like above, sem of all data together
        
        % Get means for each day
        for d=1:length(days)
            currday=days(d);
            eval([str(g,:),'daymean_stilltimeep',num2str(ep),'(',num2str(currday),')= nanmean(',str(g,:),'stilltime_ep',num2str(ep),'(:,d));']);
            eval([str(g,:),'dayerr_stilltimeep',num2str(ep),'(',num2str(currday),')= nansem(',str(g,:),'stilltime_ep',num2str(ep),'(:,d));']);
        end
        
    end
    
    
    % Get means, etc of Fraction change vs days
    
    postsleeps = [3,5];
    for ep = postsleeps
        
        % For Ripple Rate
        % --------------
        
        % Per day
        eval([str(g,:),'meanfr_riprateep' num2str(ep) '= nanmean(',str(g,:),'frriprate_ep',num2str(ep),');']) % Per day - across animals
        eval([str(g,:),'errfr_riprateep' num2str(ep) '= nansem(',str(g,:),'frriprate_ep',num2str(ep),');']) % Per day - across animals
        
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
        
        % Get means for each day
        for d=1:length(days)
            currday=days(d);
            eval([str(g,:),'daymeanfr_riprateep',num2str(ep),'(',num2str(currday),')= nanmean(',str(g,:),'frriprate_ep',num2str(ep),'(:,d));']);
            eval([str(g,:),'dayerrfr_riprateep',num2str(ep),'(',num2str(currday),')= nansem(',str(g,:),'frriprate_ep',num2str(ep),'(:,d));']);
        end
        
        
        % For Ripple Percent
        % --------------
        
        % Per day
        eval([str(g,:),'meanfr_ripperep' num2str(ep) '= nanmean(',str(g,:),'frripper_ep',num2str(ep),');']) % Per day - across animals
        eval([str(g,:),'errfr_ripperep' num2str(ep) '= nansem(',str(g,:),'frripper_ep',num2str(ep),');']) % Per day - across animals
        
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
        
        % Get means for each day
        for d=1:length(days)
            currday=days(d);
            eval([str(g,:),'daymeanfr_ripperep',num2str(ep),'(',num2str(currday),')= nanmean(',str(g,:),'frripper_ep',num2str(ep),'(:,d));']);
            eval([str(g,:),'dayerrfr_ripperep',num2str(ep),'(',num2str(currday),')= nansem(',str(g,:),'frripper_ep',num2str(ep),'(:,d));']);
        end
        
        
        
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
    
    % ***************   Ripple Size  ******************
    % *************************************************
    
    % Histogram Edges
    ripsizeedges = 3:0.5:12;
    ripntetrange = 0:6;
    % Get nanim for current grp
    eval(['totanim = length(',str(g,:),'ripsize_ep3);']);
    % Make vector to combine across days and animals: ripsize. Keep epochs separate for now. Easy to combine later
    for ep=1:length(allepochs)
        currep=allepochs(ep);
        eval([str(g,:),'ripsize_ep',num2str(currep),'_alldays=[];']); % combine across days. Keep epochs separate.
        eval([str(g,:),'ripntet_ep',num2str(currep),'_alldays=[];']); % combine across days. Keep epochs separate.
    end
    
    % Loop over each day and get separately
    for n = 1:length(days)
        
        currday = days(n);
        for ep=1:length(allepochs) % loop over the two run epochs
            currep=allepochs(ep);
            % Make vector to combine across animals for current day - keep epochs separate for now - Stats.
            eval([str(g,:),'ripsize_ep',num2str(currep),'_day{',num2str(currday),'}=[];']); % combine across animals for current day.
            eval([str(g,:),'ripntet_ep',num2str(currep),'_day{',num2str(currday),'}=[];']); % combine across animals for current day.
            
            % Get from each animal and average for day
            for an=1:totanim % Loop over anim
                
                % Get current data
                eval(['currripsize =',str(g,:),'ripsize_ep',num2str(currep),'{an}{currday};']);
                eval(['currripntet =',str(g,:),'ripntet_ep',num2str(currep),'{an}{currday};']);
                % ripple parameters
                eval(['curr_ripbaseline =',str(g,:),'ripbaseline_ep',num2str(currep),'{an}{currday};']);
                eval(['curr_ripstd =',str(g,:),'ripstd_ep',num2str(currep),'{an}{currday};']);
                
                % Make histogram for current day and epoch - ripsize
                h = histc(currripsize,ripsizeedges);
                h = cumsum(h); h = h./max(h);  % Cumulative proportion for current day and epoch
                % Store. You will take mean and sem across animals for plotting for day
                if length(find(isnan(h)))>0 % you get allnans when no ripp.
                    h(1)=0; h(2:length(ripsizeedges))=1;
                end
                eval([str(g,:),'ripsizehist_dayanim_ep',num2str(currep),'(currday,an,:) = h;']);
                % ripntet
                hn = histc(currripntet,ripntetrange);
                hn = cumsum(hn); hn = hn./max(hn);
                if length(find(isnan(h)))>0 % you get allnans.
                    hn(1)=0; hn(2:length(ripsizeedges))=1;
                end
                % Store. You will take mean and sem across animals for plotting for day
                eval([str(g,:),'ripntethist_dayanim_ep',num2str(currep),'(currday,an,:) = hn;']);
                
                % Store vectors
                % Put raw values in vector for current day from all animals - Stats for each day based on raw stds
                eval([str(g,:),'ripsize_ep',num2str(currep),'_day{',num2str(currday),'}=[',str(g,:),'ripsize_ep',num2str(currep),'_day{',num2str(currday),'},currripsize];']);
                eval([str(g,:),'ripntet_ep',num2str(currep),'_day{',num2str(currday),'}=[',str(g,:),'ripntet_ep',num2str(currep),'_day{',num2str(currday),'},currripsize];']);
                
                % Put raw values in vector for all days and all animals
                eval([str(g,:),'ripsize_ep',num2str(currep),'_alldays=[',str(g,:),'ripsize_ep',num2str(currep),'_alldays, currripsize];']);
                eval([str(g,:),'ripntet_ep',num2str(currep),'_alldays=[',str(g,:),'ripntet_ep',num2str(currep),'_alldays, currripntet];']);
                
            end % end loop over all animals for curr day and epoch - keeping epochs separate for now
        end % end loop over run epochs
    end % end loop over days for ripple size
    
    
    
    
end   % end str = Exp and Con





%****************************************
% Figures -
%*******************************************


% Can combine epoch3 and epoch5 data (= post-sleep) or plot them separate
epcomb=1;


if figopt1 == 1
    
    % ********************** Ripple Rate ****************************
    % ***************************************************************
    
    %-----------------------------------------------------------
    %-----------------------------------------------------------
    %  1) Bar Plot - Mean Ripple Rates 2) vs days  - Mean Ripple Rates
    
    
    % 1) Bar plot for all days
    figure; hold on;
    if forppr==1
        redimscreen_figforppr1;
    else
        redimscreen_figforppt1;
    end
    
    % Epoch 2 and 4 combined
    %-----------------------
    if epcomb==1
        % Plot as bars
        %         bar(0.5,mean(Conriprate_ep1(:)),'b');
        %         bar(1.5,mean(Expriprate_ep1(:)),'r');
        %         bar(2.5,mean([Conriprate_ep3(:); Conriprate_ep5(:)]),'b');
        %         bar(3.5,mean([Expriprate_ep3(:); Expriprate_ep5(:)]),'r');
        %         errorbar(0.5,mean(Conriprate_ep1(:)),sem(Conriprate_ep1(:)),'k','LineWidth',2);
        %         errorbar(1.5,mean(Expriprate_ep1(:)),sem(Expriprate_ep1(:)),'k','LineWidth',2);
        %         errorbar(2.5,mean([Conriprate_ep3(:); Conriprate_ep5(:)]),sem([Conriprate_ep3(:); Conriprate_ep5(:)]),'k','LineWidth',2);
        %         errorbar(3.5,mean([Expriprate_ep3(:); Expriprate_ep5(:)]),sem([Expriprate_ep3(:); Expriprate_ep5(:)]),'k','LineWidth',2);
        % Plot as groups
        errorbar([1,2], [mean(Conriprate_ep1(:));mean([Conriprate_ep3(:); Conriprate_ep5(:)])],...
            [sem(Conriprate_ep1(:));sem([Conriprate_ep3(:); Conriprate_ep5(:)])],'bd-','MarkerSize',12,'LineWidth',3);
        errorbar([1,2], [mean(Expriprate_ep1(:));mean([Expriprate_ep3(:); Expriprate_ep5(:)])],...
            [sem(Expriprate_ep1(:));sem([Expriprate_ep3(:); Expriprate_ep5(:)])],'ro-','MarkerSize',12,'LineWidth',3);
        % ----------- Stats -------------
        %Nanimals*Ndays  pts in each group: each day is a measure
        % Con against Exp for epoch1 and epochs3and5 combined
        % ---------------
        [h11_allriprate,p11_allriprate] = ttest2(Conriprate_ep1(:),Expriprate_ep1(:)); % p=0.3779
        [h3535_allriprate,p3535_allriprate] = ttest2(mean([Conriprate_ep3(:), Conriprate_ep5(:)],2),mean([Expriprate_ep3(:), Expriprate_ep5(:)],2)); %p=0.6894,32X32
        [h135_allriprateCon,p135_allriprateCon] = ttest2(Conriprate_ep1(:),mean([Conriprate_ep3(:), Conriprate_ep5(:)],2)); %p=0.2005,32X32
        [h135_allriprateExp,p135_allriprateExp] = ttest2(Expriprate_ep1(:),mean([Expriprate_ep3(:), Expriprate_ep5(:)],2)); % p=0.6535,32X32
        % Use anova_rm/anova2/anovan/rm_anova2
        % anova_rm: pgrp=0.9605; ptime=0; pint=0.1689
        [panovarm_allriprate] = anova_rm({[Conriprate_ep1(:),mean([Conriprate_ep3(:),Conriprate_ep5(:)],2)] [Expriprate_ep1(:),mean([Expriprate_ep3(:),Expriprate_ep5(:)],2)]},'off');
        
%         % anova2 with epochs3and5 meaned. each day is a measure. factors = Grp and Pre-Post condition. 4*8=32 repeats
%         % anova2: pgrp=0.3553,prows(pre-post)=0.2642,pint=0.6981
%         [panova2_allriprate, table, stats] = anova2(  [[Conriprate_ep1(:);mean([Conriprate_ep3(:),Conriprate_ep5(:)],2)],[Expriprate_ep1(:);mean([Expriprate_ep3(:),Expriprate_ep5(:)],2)]],32,'on');
%         % multcompare on the anova2
%         
%         % Using tukey-kramer = hsd by default. Can use 'ctype','bonferroni'
%         % 1a) - % compare grps: Grps same
%         [c,m,h,nms] = multcompare(stats,'estimate','column','display','off');
%         % 1b) - - % compare Pre Post. Same
%         [c,m,h,nms] = multcompare(stats,'estimate','row','display','off');
        
        % rm_anova2 will treat both day and time as factors for different subjects. Wont treat day as independent factor
        
        % Plot midway between the two bars
        if h11_allriprate==1,
            mul = sign(mean(Conriprate_ep1(:)));
            plot(1.0, (mean(Conriprate_ep1(:))) + sem(Conriprate_ep1(:)) + ...
                1.1*mul*(sem(Conriprate_ep1(:))), 'r*','MarkerSize',12,'LineWidth',2);
        end
        % Plot midway between the two bars
        if h3535_allriprate==1,
            mul = sign(mean([Conriprate_ep3(:); Conriprate_ep5(:)]));
            plot(3.0, (mean([Conriprate_ep3(:); Conriprate_ep5(:)]) + sem([Conriprate_ep3(:); Conriprate_ep5(:)]) + ...
                1.1*mul*(sem([Conriprate_ep3(:); Conriprate_ep5(:)]))), 'r*','MarkerSize',12,'LineWidth',2);
        end
        % ConvsCon and ExpvsExp
        % -------------------------
        [h135Con_allriprate,p135Con_allriprate] = ttest2(Conriprate_ep1(:),[Conriprate_ep3(:); Conriprate_ep5(:)]); % p=0.1340
        [h135Exp_allriprate,p135Exp_allriprate] = ttest2(Expriprate_ep1(:),[Expriprate_ep3(:); Expriprate_ep5(:)]); % p=0.5980
        % Plot over Post
        if h135Con_allriprate==1,
            mul = sign(mean([Conriprate_ep3(:); Conriprate_ep5(:)]));
            plot(2.5, (mean([Conriprate_ep3(:); Conriprate_ep5(:)]) + sem([Conriprate_ep3(:); Conriprate_ep5(:)]) + ...
                1.2*mul*sem([Conriprate_ep3(:); Conriprate_ep5(:)])), 'r+','MarkerSize',12,'LineWidth',2);
        end
        if h135Exp_allriprate==1,
            mul = sign(mean([Expriprate_ep3(:); Expriprate_ep5(:)]));
            plot(3.5, (mean([Expriprate_ep3(:); Expriprate_ep5(:)]) + sem([Expriprate_ep3(:); Expriprate_ep5(:)]) + ...
                1.2*mul*sem([Expriprate_ep3(:); Expriprate_ep5(:)])), 'r+','MarkerSize',12,'LineWidth',2);
        end
        title('Ripple Rate during Rest: All Days');
        ylabel('Ripple Rate (Hz)');
        set(gca,'XTick',[1,2],'XTickLabel',{'Pre';'Post'},'FontSize',xfont,'Fontweight','normal');
        set(gca,'YLim',[0 1.2]);
        %set(gca,'XLim',[0.5 3.5]);
        %axis([0 4 0 0.5])
        if savefig1==1,
            figfile = [figdir,'0RippleRateRest_PostComb'];
            print('-dpdf', figfile);
            print('-djpeg', figfile);
            saveas(gcf,figfile,'fig');
        end
        
        
    else
        for ep=1:length(sleeps)
            currep=sleeps(ep);
            eval(['allmean_riprate(',num2str(ep),',:) = [Conallmean_riprateep',num2str(currep),';  Expallmean_riprateep',num2str(currep),'];'])
        end
        
        % Plot as bars
        %         bar(allmean_riprate,'grouped');
        %         errorbar(0.8,Conallmean_riprateep1,Conallerr_riprateep1,'k');
        %         errorbar(1.2,Expallmean_riprateep1,Expallerr_riprateep1,'k');
        %         errorbar(1.8,Conallmean_riprateep3,Conallerr_riprateep3,'k');
        %         errorbar(2.2,Expallmean_riprateep3,Expallerr_riprateep3,'k');
        %         errorbar(2.8,Conallmean_riprateep5,Conallerr_riprateep5,'k');
        %         errorbar(3.2,Expallmean_riprateep5,Expallerr_riprateep5,'k');
        % Plot as groups
        errorbar([1:3], [mean(Conriprate_ep1(:));mean(Conriprate_ep3(:)); mean(Conriprate_ep5(:))],...
            [sem(Conriprate_ep1(:));sem(Conriprate_ep3(:)); sem(Conriprate_ep5(:))],'bd-','MarkerSize',12,'LineWidth',3);
        errorbar([1:3], [mean(Expriprate_ep1(:));mean(Expriprate_ep3(:)); mean(Expriprate_ep5(:))],...
            [sem(Expriprate_ep1(:));sem(Expriprate_ep3(:)); sem(Expriprate_ep5(:))],'ro-','MarkerSize',12,'LineWidth',3);
        % ----------- Stats -------------
        % Comparing animals: each animal has one number: average across days
        % Con against Exp for each epoch
        [h11_allriprate,p11_allriprate] = ttest2(Conriprate_ep1(:),Expriprate_ep1(:));
        [h33_allriprate,p33_allriprate] = ttest2(Conriprate_ep3(:),Expriprate_ep3(:));
        [h55_allriprate,p55_allriprate] = ttest2(Conriprate_ep5(:),Expriprate_ep5(:));
        % Use anova_rm/anova2/anovan/rm_anova2
        % anova_rm: pgrp=0.5644; ptime(pre-post1-post2)=0; pint=0.1092
        [panova_allriprate] = anova_rm({[Conriprate_ep1(:),Conriprate_ep3(:),Conriprate_ep5(:)] [Expriprate_ep1(:),Expriprate_ep3(:),Expriprate_ep5(:)]},'off');
        % anova2. each day is a measure. factors = Grp and Pre-Post condition. 4*8=32 repeats
        % anova2: pgrp=0.325,prows(pre-post1-post2)=0.3964,pint=0.8935
        [panova2_allriprate, table, stats, terms] = anova2([[Conriprate_ep1(:);Conriprate_ep3(:);Conriprate_ep5(:)],[Expriprate_ep1(:);Expriprate_ep3(:);Expriprate_ep5(:)]],32,'off');
        
        % rm_anova2 will treat both day and time as factors for different subjects. Wont treat day as independent factor
        
        
        
        % Plot midway between the two graphs
        if h11_allriprate==1,
            mul = sign(mean(Conriprate_ep1(:)));
            plot(1, (mean(Conriprate_ep1(:))+sem(Conriprate_ep1(:))+1.2*mul*sem(Conriprate_ep1(:))), 'r*','MarkerSize',12,'LineWidth',2);
        end
        if h33_allriprate==1,
            mul = sign(mean(Conriprate_ep3(:)));
            plot(2, (mean(Conriprate_ep3(:))+sem(Conriprate_ep1(:))+1.2*mul*sem(Conriprate_ep3(:))), 'r*','MarkerSize',12,'LineWidth',2);
        end
        if h55_allriprate==1,
            mul = sign(Conallmean_riprateep5);
            plot(3, (mean(Conriprate_ep5(:))+sem(Conriprate_ep5(:))+1.2*mul*sem(Conriprate_ep5(:))), 'r*','MarkerSize',12,'LineWidth',2);
        end
        % Con against Con for different epochs
        [h13Con_allriprate,p13Con_allriprate] = ttest2(Conriprate_ep1(:),Conriprate_ep3(:));
        [h15Con_allriprate,p15Con_allriprate] = ttest2(Conriprate_ep1(:),Conriprate_ep3(:));
        [h35Con_allriprate,p35Con_allriprate] = ttest2(Conriprate_ep1(:),Conriprate_ep3(:));
        if h13Con_allriprate==1,
            mul = sign(mean(Conriprate_ep3(:)));
            plot(1.8, (mean(Conriprate_ep3(:))+sem(Conriprate_ep3(:))+1.2*mul*sem(Conriprate_ep3(:))), 'r+','MarkerSize',12,'LineWidth',2);
        end
        if h15Con_allriprate==1,
            mul = sign(mean(Conriprate_ep5(:)));
            plot(2.8, (mean(Conriprate_ep5(:))+sem(Conriprate_ep3(:))+1.2*mul*sem(Conriprate_ep5(:))), 'r+','MarkerSize',12,'LineWidth',2);
        end
        if h35Con_allriprate==1,
            mul = sign(mean(Conriprate_ep5(:)));
            plot(2.8, (mean(Conriprate_ep5(:))+sem(Conriprate_ep3(:))+1.4*mul*sem(Conriprate_ep5(:))), 'm+','MarkerSize',12,'LineWidth',2);
        end
        % Exp against Exp for different epochs
        [h13Exp_allriprate,p13Exp_allriprate] = ttest2(Expriprate_ep1(:),Expriprate_ep3(:));
        [h15Exp_allriprate,p15Exp_allriprate] = ttest2(Expriprate_ep1(:),Expriprate_ep3(:));
        [h35Exp_allriprate,p35Exp_allriprate] = ttest2(Expriprate_ep1(:),Expriprate_ep3(:));
        if h13Exp_allriprate==1,
            mul = sign(mean(Expriprate_ep3(:)));
            plot(2.2, (mean(Expriprate_ep3(:))+sem(Expriprate_ep3(:))+1.2*mul*sem(Expriprate_ep3(:))), 'r+','MarkerSize',12,'LineWidth',2);
        end
        if h15Exp_allriprate==1,
            mul = sign(mean(Expriprate_ep5(:)));
            plot(3.2, (mean(Expriprate_ep5(:))+sem(Expriprate_ep3(:))+1.2*mul*sem(Expriprate_ep5(:))), 'r+','MarkerSize',12,'LineWidth',2);
        end
        if h35Exp_allriprate==1,
            mul = sign(mean(Expriprate_ep5(:)));
            plot(3.2, (mean(Expriprate_ep5(:))+sem(Expriprate_ep3(:))+1.4*mul*sem(Expriprate_ep5(:))), 'm+','MarkerSize',12,'LineWidth',2);
        end
        title('Ripple Rate during Sleep: All Days');
        ylabel('Ripple Rate (Hz)');
        set(gca,'XTick',[1:3],'XTickLabel',{'Pre';'Post1';'Post2'},'FontSize',16,'Fontweight','normal');
        set(gca,'YLim',[0 1.2]);
        %set(gca,'XLim',[0.5 3.5]);
        %axis([0 4 0 0.5])
        if savefig1==1,
            figfile = [figdir,'RippleRateRest_PostSep'];
            print('-dpdf', figfile);
            print('-djpeg', figfile);
            saveas(gcf,figfile,'fig');
        end
        
    end % epcomb
    
    %-----------------------------------------------------------
    %-----------------------------------------------------------
    
    % 2) Mean Ripple Rates vs days
    
    if epcomb==1
        % Epoch 3 and 5 combined
        %-----------------------
        figure; hold on;
        redimscreen_2horsubplots;
        % First combine Ep3and5 data for both Con and Exp and do stats
        for d=1:length(days)
            currday = days(d);
            Conday_riprate_post(:,currday) = [Conriprate_ep3(:,currday); Conriprate_ep5(:,currday)];
            Expday_riprate_post(:,currday) = [Expriprate_ep3(:,currday); Expriprate_ep5(:,currday)];
            Conday_riprate_pre(:,currday) = Conriprate_ep1(:,currday);
            Expday_riprate_pre(:,currday) = Expriprate_ep1(:,currday);
            [p_dayriprate_post(currday),h_dayriprate_post(currday)] = ranksum(Conday_riprate_post(:,currday),Expday_riprate_post(:,currday));
            [p_dayriprate_pre(currday),h_dayriprate_pre(currday)] = ranksum(Conday_riprate_pre(:,currday),Expday_riprate_pre(:,currday));
        end
        % Now plot
        subplot(1,2,1); hold on;
        errorbar(1:length(Conday_riprate_pre),mean(Conday_riprate_pre,1),sem(Conday_riprate_pre,1),'bd-','MarkerSize',10,'LineWidth',2);
        errorbar(1:length(Expday_riprate_pre),mean(Expday_riprate_pre,1),sem(Expday_riprate_pre,1),'ro-','MarkerSize',10,'LineWidth',2);
        % Anova
        % repeated measures anova. Grp X Day Compare separately, the pre curves for con and exp; and the post curves for con and exp
        % pgrp=0.7422,pday=0.0349,pint=0.6111
        [panova_allripratevsday0] = anova_rm({Conday_riprate_pre Expday_riprate_pre},'off');
        title('Ripple Rate during Rest-Pre');
        ylabel('Ripple Rate (Hz)'); xlabel('Days');
        % ----------- Stats -------------
        % Comparing groups for each day:
        for d=1:length(days)
            currday = days(d);
            % Plot * if significant
            if h_dayriprate_pre(currday)==1,
                mul = sign(mean(Conday_riprate_pre(:,currday)));
                plot(currday, mean(Conday_riprate_pre(:,currday)) + sem(Conday_riprate_pre(:,currday))...
                    + 1.1*mul*sem(Conday_riprate_pre(:,currday)), 'r*','MarkerSize',12);
            end
        end
        % -------------------------
        subplot(1,2,2); hold on;
        errorbar(1:length(Conday_riprate_post),mean(Conday_riprate_post,1),sem(Conday_riprate_post,1),'bd-','MarkerSize',10,'LineWidth',2);
        errorbar(1:size(Expday_riprate_post,2),mean(Expday_riprate_post,1),sem(Expday_riprate_post,1),'ro-','MarkerSize',10,'LineWidth',2);
        % Anova
        % repeated measures anova. Grp X Day Compare separately, the pre curves for con and exp; and the post curves for con and exp
        % pgrp=0.8085,pday=0,pint=0.0.0711
        [panova_allripratevsday1] = anova_rm({Conday_riprate_post Expday_riprate_post},'off');
        title('Ripple Rate during Rest-Post');
        ylabel('Ripple Rate (Hz)'); xlabel('Days')
        %set(gca,'XLim',[0.5 3.5]);
        %axis([0 4 0 0.5])
        % ----------- Stats -------------
        % Comparing groups for each day:
        for d=1:length(days)
            currday = days(d);
            % Plot * if significant
            if h_dayriprate_post(currday)==1,
                mul = sign(mean(Conday_riprate_post(:,currday)));
                plot(currday, mean(Conday_riprate_post(:,currday)) + sem(Conday_riprate_post(:,currday))...
                    + 1.1*mul*sem(Conday_riprate_post(:,currday)), 'r*','MarkerSize',12);
            end
        end
        % -------------------------
        if savefig1==1,
            figfile = [figdir,'RippleRateRestVsDays_PostComb'];
            print('-dpdf', figfile);
            print('-djpeg', figfile);
            saveas(gcf,figfile,'fig');
        end
    else
        % Epoch 3 and 5 separate
        %-----------------------
        figure; hold on;
        redimscreen_2horsubplots;
        % First do stats
        for d=1:length(days)
            currday = days(d);
            Conday_riprate_pre(:,currday) = Conriprate_ep1(:,currday); Expday_riprate_pre(:,currday) = Expriprate_ep1(:,currday);
            Conday_riprate_post1(:,currday) = Conriprate_ep3(:,currday); Expday_riprate_post1(:,currday) = Expriprate_ep3(:,currday);
            Conday_riprate_post2(:,currday) = Conriprate_ep5(:,currday); Expday_riprate_post2(:,currday) = Expriprate_ep5(:,currday);
            [p_dayriprate_pre(currday),h_dayriprate_pre(currday)] = ranksum(Conriprate_ep1(:,currday),Expriprate_ep1(:,currday));
            [p_dayriprate_post1(currday),h_dayriprate_post1(currday)] = ranksum(Conriprate_ep3(:,currday),Expriprate_ep3(:,currday));
            [p_dayriprate_post2(currday),h_dayriprate_post2(currday)] = ranksum(Conriprate_ep5(:,currday),Expriprate_ep5(:,currday));
        end
        % Now plot
        subplot(1,3,1); hold on;
        errorbar(1:length(Conday_riprate_pre),mean(Conday_riprate_pre,1),sem(Conday_riprate_pre,1),'bd-','MarkerSize',10,'LineWidth',2);
        errorbar(1:length(Expday_riprate_pre),mean(Expday_riprate_pre,1),sem(Expday_riprate_pre,1),'ro-','MarkerSize',10,'LineWidth',2);
        % Anova
        % repeated measures anova. Grp X Day Compare separately, the pre curves for con and exp; and the post1,post2 curves for con and exp
        % pgrp=0.7422,pday=0.0349,pint=0.6111
        [panova_allripratevsday0] = anova_rm({Conday_riprate_pre Expday_riprate_pre},'off');
        title('Ripple Rate during Rest-Pre');
        ylabel('Ripple Rate (Hz)'); xlabel('Days');
        set(gca,'YLim',[0 2.4]);
        % ----------- Stats -------------
        % Comparing groups for each day:
        for d=1:length(days)
            currday = days(d);
            % Plot * if significant
            if h_dayriprate_pre(currday)==1,
                mul = sign(mean(Conday_riprate_pre(:,currday)));
                plot(currday, mean(Conday_riprate_pre(:,currday)) + sem(Conday_riprate_pre(:,currday))...
                    + 1.1*mul*sem(Conday_riprate_pre(:,currday)), 'r*','MarkerSize',12);
            end
        end
        % -------------------------
        subplot(1,3,2); hold on;
        errorbar(1:length(Conday_riprate_post1),mean(Conday_riprate_post1,1),sem(Conday_riprate_post1,1),'bd-','MarkerSize',10,'LineWidth',2);
        errorbar(1:length(Expday_riprate_post1),mean(Expday_riprate_post1,1),sem(Expday_riprate_post1,1),'ro-','MarkerSize',10,'LineWidth',2);
        % Anova
        % repeated measures anova. Grp X Day Compare separately, the pre curves for con and exp; and the post1,post2 curves for con and exp
        % pgrp=0.9125,pday=0.0462,pint=0
        [panova_allripratevsday1] = anova_rm({Conday_riprate_post1 Expday_riprate_post1},'off');
        title('Ripple Rate during Rest-Post1');
        ylabel('Ripple Rate (Hz)'); xlabel('Days')
        %set(gca,'XLim',[0.5 3.5]);
        %axis([0 4 0 0.5])
        % ----------- Stats -------------
        % Comparing groups for each day:
        for d=1:length(days)
            currday = days(d);
            % Plot * if significant
            if h_dayriprate_post2(currday)==1,
                mul = sign(mean(Conday_riprate_post2(:,currday)));
                plot(currday, mean(Conday_riprate_post2(:,currday)) + sem(Conday_riprate_post2(:,currday))...
                    + 1.1*mul*sem(Conday_riprate_post2(:,currday)), 'r*','MarkerSize',12);
            end
        end
        % -------------------------
        subplot(1,3,3); hold on;
        errorbar(1:length(Conday_riprate_post2),mean(Conday_riprate_post2,1),sem(Conday_riprate_post2,1),'bd-','MarkerSize',10,'LineWidth',2);
        errorbar(1:length(Expday_riprate_post2),mean(Expday_riprate_post2,1),sem(Expday_riprate_post2,1),'ro-','MarkerSize',10,'LineWidth',2);
        % Anova
        % repeated measures anova. Grp X Day Compare separately, the pre curves for con and exp; and the post1,post2 curves for con and exp
        % pgrp=0.8375,pday=0.0363,pint=0
        [panova_allripratevsday2] = anova_rm({Conday_riprate_post2 Expday_riprate_post2},'off');
        title('Ripple Rate during Rest-Post2');
        ylabel('Ripple Rate (Hz)'); xlabel('Days')
        set(gca,'YLim',[0 2.4]);
        %set(gca,'XLim',[0.5 3.5]);
        %axis([0 4 0 0.5])
        % ----------- Stats -------------
        % Comparing groups for each day:
        for d=1:length(days)
            currday = days(d);
            % Plot * if significant
            if h_dayriprate_post2(currday)==1,
                mul = sign(mean(Conday_riprate_post2(:,currday)));
                plot(currday, mean(Conday_riprate_post2(:,currday)) + sem(Conday_riprate_post2(:,currday))...
                    + 1.1*mul*sem(Conday_riprate_post2(:,currday)), 'r*','MarkerSize',12);
            end
        end
        % -------------------------
        if savefig1==1,
            figfile = [figdir,'RippleRateRestVsDays_PostSep'];
            print('-dpdf', figfile);
            print('-djpeg', figfile);
            saveas(gcf,figfile,'fig');
        end
    end % end epcomb
    
    
    
    % ************ Fraction Change In Ripple Rate *******************
    % ***************************************************************
    
    %-----------------------------------------------------------
    %-----------------------------------------------------------
    
    % 1) Bar plot for all days, & 2) vs days
    
    % For all days
    figure; hold on;
    if forppr==1
        redimscreen_figforppr1;
    else
        redimscreen_figforppt1;
    end
    
    % Epoch 3 and 5 combined
    %-----------------------
    if epcomb==1
        % Plot bars
        %         bar(1,mean([Confrriprate_ep3(:); Confrriprate_ep5(:)]),'b');
        %         bar(2,mean([Expfrriprate_ep3(:); Expfrriprate_ep5(:)]),'r');
        %         errorbar(1,mean([Confrriprate_ep3(:); Confrriprate_ep5(:)]),sem([Confrriprate_ep3(:); Confrriprate_ep5(:)]),'k','LineWidth',2);
        %         errorbar(2,mean([Expfrriprate_ep3(:); Expfrriprate_ep5(:)]),sem([Expfrriprate_ep3(:); Expfrriprate_ep5(:)]),'k','LineWidth',2);
        % Plot as groups
        errorbar(0.5, mean([Confrriprate_ep3(:); Confrriprate_ep5(:)]),...
            sem([Confrriprate_ep3(:); Confrriprate_ep5(:)]),'bd','MarkerSize',12,'LineWidth',3);
        errorbar(1.5, mean([Expfrriprate_ep3(:); Expfrriprate_ep5(:)]),...
            sem([Expfrriprate_ep3(:); Expfrriprate_ep5(:)]),'ro','MarkerSize',12,'LineWidth',3);
        % ----------- Stats -------------
        %Nanimals*Ndays  pts in each group: each day is a measure
        % Con against Exp for epoch1 and epochs3and5 combined
        % ---------------
        [h3535_allripratefr,p3535_allripratefr] = ttest2([Confrriprate_ep3(:); Confrriprate_ep5(:)],[Expfrriprate_ep3(:); Expfrriprate_ep5(:)]);
        % Use anova
        % anova1: pgrp=0.0339 (same value as anova_rm); F=4.71
%        [panova_allripratefr] = anova1([mean([Confrriprate_ep3(:),Confrriprate_ep5(:)],2),mean([Expfrriprate_ep3(:),Expfrriprate_ep5(:)],2)],[],'off');
        % anova_rm will compare curves for the two groups across days (prevspost fraction change vs days)
        % Do below for VsDay - ConFrPost/Pre vs ExpFrPost/Pre
        % Plot midway between the two bars
        if h3535_allripratefr==1,
            mul = sign(mean([Confrriprate_ep3(:); Confrriprate_ep5(:)]));
            plot(1, (mean([Confrriprate_ep3(:); Confrriprate_ep5(:)]) + sem([Confrriprate_ep3(:); Confrriprate_ep5(:)]) + ...
                1.1*mul*(sem([Confrriprate_ep3(:); Confrriprate_ep5(:)]))), 'r*','MarkerSize',12,'Linewidth',2);
        end
        title('Fraction change in Ripple Rate during Rest: All Days');
        ylabel('Fraction Change (Post/Pre)');
        set(gca,'XTick',[0.5,1.5],'XTickLabel',{'Con';'Exp'},'FontSize',xfont,'Fontweight','normal');
        set(gca,'YLim',[0 2.2]);
        %set(gca,'XLim',[0.5 3.5]);
        %axis([0 4 0 0.5])
        if savefig1==1,
            figfile = [figdir,'RippleRateRest_FracChg_PostComb'];
            print('-dpdf', figfile);
            print('-djpeg', figfile);
            saveas(gcf,figfile,'fig');
        end
    else
        
        % Plot Bars
        %         bar(0.5,mean(Confrriprate_ep3(:)),'b');
        %         bar(1.5,mean(Expfrriprate_ep3(:)),'r');
        %         bar(2.5,mean(Confrriprate_ep5(:)),'b');
        %         bar(3.5,mean(Expfrriprate_ep5(:)),'r');
        %         errorbar(0.5,mean(Confrriprate_ep3(:)),sem(Confrriprate_ep3(:)),'k','LineWidth',2);
        %         errorbar(1.5,mean(Expfrriprate_ep3(:)),sem(Expfrriprate_ep3(:)),'k','LineWidth',2);
        %         errorbar(2.5,mean(Confrriprate_ep5(:)),sem(Confrriprate_ep5(:)),'k','LineWidth',2);
        %         errorbar(3.5,mean(Expfrriprate_ep5(:)),sem(Expfrriprate_ep5(:)),'k','LineWidth',2);
        % Plot as groups
        errorbar([1:2], [mean(Confrriprate_ep3(:)); mean(Confrriprate_ep5(:))],...
            [sem(Confrriprate_ep3(:)); sem(Confrriprate_ep5(:))],'bd-','MarkerSize',12,'LineWidth',3);
        errorbar([1:2], [mean(Expfrriprate_ep3(:)); mean(Expfrriprate_ep5(:))],...
            [sem(Expfrriprate_ep3(:)); sem(Expfrriprate_ep5(:))],'ro-','MarkerSize',12,'LineWidth',3);
        % ----------- Stats -------------
        % Con vs Exp for each epoch
        [h33_allripratefr,p33_allripratefr] = ttest2(Confrriprate_ep3(:),Expfrriprate_ep3(:));
        [h55_allripratefr,p55_allripratefr] = ttest2(Confrriprate_ep5(:),Expfrriprate_ep5(:));
        % Use anova_rm/anova2/anovan/rm_anova2
        % anova_rm: pgrp=0.0339; ptime(post1/pre-post2/pre)=0.0103; pint=0.3562
        [panova_allripratefr] = anova_rm({[Confrriprate_ep3(:),Confrriprate_ep5(:)] [Expfrriprate_ep3(:),Expfrriprate_ep5(:)]},'off');
        % anova2. each day is a measure. factors = Grp and Pre-Post condition. 4*8=32 repeats
        % anova2: pcolumns(grp)=0.0058,prows(post1-post2)=0.1310,pint=0.5952
        [panova2_allripratefr] = anova2([[Confrriprate_ep3(:);Confrriprate_ep5(:)],[Expfrriprate_ep3(:);Expfrriprate_ep5(:)]],32,'off');
        % rm_anova2 will treat both day and time as factors for different subjects. Wont treat day as independent factor
        % Plot midway
        if h33_allripratefr==1,
            mul = sign(mean(Confrriprate_ep3(:)));
            plot(1, (mean(Confrriprate_ep3(:))+sem(Confrriprate_ep3(:))+1.2*mul*sem(Confrriprate_ep3(:))), 'r*','MarkerSize',12,'LineWidth',2);
        end
        if h55_allripratefr==1,
            mul = sign(mean(Confrriprate_ep5(:)));
            plot(2, (mean(Confrriprate_ep5(:))+sem(Confrriprate_ep5(:))+1.2*mul*sem(Confrriprate_ep5(:))), 'r*','MarkerSize',12,'LineWidth',2);
        end
        % ConvsCon and ExpvsExp for post epochs
        [h35Con_allripratefr,p35Con_allripratefr] = ttest2(Confrriprate_ep3(:),Confrriprate_ep5(:));
        [h35Exp_allripratefr,p35Exp_allripratefr] = ttest2(Expfrriprate_ep3(:),Expfrriprate_ep5(:));
        %Plot above Post2
        if h35Con_allripratefr==1,
            mul = sign(mean(Confrriprate_ep5(:)));
            plot(2.1, (mean(Confrriprate_ep5(:))+sem(Confrriprate_ep5(:))), 'b+','MarkerSize',12,'LineWidth',2);
        end
        if h35Exp_allripratefr==1,
            mul = sign(mean(Expfrriprate_ep5(:)));
            plot(2.1, (mean(Expfrriprate_ep5(:))+sem(Expfrriprate_ep5(:))), 'r+','MarkerSize',12,'LineWidth',2);
        end
        
        % -------------------------
        title('Fraction Change in Ripple Rate during Sleep: All Days');
        ylabel('Fraction Change (Post/Pre)');
        set(gca,'XTick',[1,2],'XTickLabel',{'Post1';'Post2'},'FontSize',16,'Fontweight','normal');
        set(gca,'YLim',[0 2.4]);
        %set(gca,'XLim',[0.5 3.5]);
        %axis([0 4 0 0.5])
        if savefig1==1,
            figfile = [figdir,'RippleRateRest_FracChg_PostSep'];
            print('-dpdf', figfile);
            print('-djpeg', figfile);
            saveas(gcf,figfile,'fig');
        end
    end % end epcomb - Frac Change In ripple rate - All Days
    
    
    
    
    %-----------------------------------------------------------
    %-----------------------------------------------------------
    %    Vs Days: Fraction Change In Ripple Rate
    
    
    % Epoch 3 and 5 combined
    %-----------------------
    if epcomb==1
        figure; hold on;
        if forppr==1
            redimscreen_figforppr1;
        else
            redimscreen_figforppt1;
        end
        % First combine Ep3and5 data for both Con and Exp and do stats
        for d=1:length(days)
            currday = days(d);
            Confrday_riprate_post(:,currday) = [Confrriprate_ep3(:,currday); Conriprate_ep5(:,currday)];
            Expfrday_riprate_post(:,currday) = [Expfrriprate_ep3(:,currday); Expfrriprate_ep5(:,currday)];
            [p_dayripratefr_post(currday),h_dayripratefr_post(currday)] = ranksum(Confrday_riprate_post(:,currday),Expfrday_riprate_post(:,currday));
            % mean of post1 and post2 - for anova
            Confrday_riprate_posta(:,currday) = mean([Confrriprate_ep3(:,currday), Conriprate_ep5(:,currday)],2);
            Expfrday_riprate_posta(:,currday) = mean([Expfrriprate_ep3(:,currday), Expriprate_ep5(:,currday)],2);
        end
        % Now plot
        errorbar(1:length(Confrday_riprate_post),mean(Confrday_riprate_post,1),sem(Confrday_riprate_post,1),'bd-','MarkerSize',10,'LineWidth',3);
        errorbar(1:size(Expfrday_riprate_post,2),mean(Expfrday_riprate_post,1),sem(Expfrday_riprate_post,1),'ro-','MarkerSize',10,'LineWidth',3);
        % Anova
        % repeated measures anova. Grp X Day Compare the post/pre curves for con and exp
        % Pgrp=0.2829, ptime=0.3801, pint=0.583,
        [panova_allripratevsdayfr] = anova_rm({Confrday_riprate_posta Expfrday_riprate_posta},'off');
        % ----------- Stats -------------
        % Comparing groups for each day:
        for d=1:length(days)
            currday = days(d);
            % Plot * if significant
            if h_dayriprate_post(currday)==1,
                mul = sign(mean(Confrday_riprate_pre(:,currday)));
                plot(currday, mean(Confrday_riprate_pre(:,currday)) + sem(Confrday_riprate_pre(:,currday))...
                    + 1.1*mul*sem(Confrday_riprate_pre(:,currday)), 'r*','MarkerSize',12);
            end
        end
        % -------------------------
        title('Fraction change in Ripple Rate during Rest');
        ylabel('Fraction Change in Ripple Rate (Post/Pre)'); xlabel('Days');
        set(gca,'YLim',[0 2.1]);
        if savefig1==1,
            figfile = [figdir,'0RippleRateRest_FracChgVsDay_PostComb'];
            print('-dpdf', figfile);
            print('-djpeg', figfile);
            saveas(gcf,figfile,'fig');
        end
    else
        
        figure; hold on;
        redimscreen_2horsubplots;
        % First do stats
        for d=1:length(days)
            currday = days(d);
            Confrday_riprate_post1(:,currday) = Confrriprate_ep3(:,currday); Expfrday_riprate_post1(:,currday) = Expfrriprate_ep3(:,currday);
            Confrday_riprate_post2(:,currday) = Confrriprate_ep5(:,currday); Expfrday_riprate_post2(:,currday) = Expfrriprate_ep5(:,currday);
            [p_dayripratefr_post1(currday),h_dayripratefr_post1(currday)] = ranksum(Confrriprate_ep3(:,currday),Expfrriprate_ep3(:,currday));
            [p_dayripratefr_post2(currday),h_dayripratefr_post2(currday)] = ranksum(Confrriprate_ep5(:,currday),Expfrriprate_ep5(:,currday));
        end
        % Anova
        % repeated measures anova. Grp X Day Compare separately, the post1/pre curve for con and exp; and the post2/pre curve for con and exp
        [panova_allripratevsdayfr1] = anova_rm({Confrday_riprate_post1 Expfrday_riprate_post1},'off');
        % Pgrp=0.219, ptime=0.3359, pint=0.6543,
        [panova_allripratevsdayfr2] = anova_rm({Confrday_riprate_post2 Expfrday_riprate_post2},'off');
        % Pgrp=0.4723, ptime=0.219, pint=0.5863,
        % Now plot
        subplot(1,2,1); hold on;
        errorbar(1:length(Confrday_riprate_post1),mean(Confrday_riprate_post1,1),sem(Confrday_riprate_post1,1),'bd-','MarkerSize',10,'LineWidth',2);
        errorbar(1:length(Expfrday_riprate_post1),mean(Expfrday_riprate_post1,1),sem(Expfrday_riprate_post1,1),'ro-','MarkerSize',10,'LineWidth',2);
        title('Fraction change in Ripple Rate during Post1');
        ylabel('Fraction Change in Ripple Rate (Post1/Pre)'); xlabel('Days');
        set(gca,'YLim',[0 3.1]);
        % ----------- Stats -------------
        % Comparing groups for each day:
        for d=1:length(days)
            currday = days(d);
            % Plot * if significant
            if h_dayripratefr_post1(currday)==1,
                mul = sign(mean(Confrday_riprate_post1(:,currday)));
                plot(currday, mean(Confrday_riprate_post1(:,currday)) + sem(Confrday_riprate_post1(:,currday))...
                    + 1.1*mul*sem(Confrday_riprate_post1(:,currday)), 'r*','MarkerSize',12);
            end
        end
        % -------------------------
        subplot(1,2,2); hold on;
        errorbar(1:length(Confrday_riprate_post2),mean(Confrday_riprate_post2,1),sem(Confrday_riprate_post2,1),'bd-','MarkerSize',10,'LineWidth',2);
        errorbar(1:length(Expfrday_riprate_post2),mean(Expfrday_riprate_post2,1),sem(Expfrday_riprate_post2,1),'ro-','MarkerSize',10,'LineWidth',2);
        title('Fraction change in Ripple Rate during Post2');
        ylabel('Fraction Change in Ripple Rate (Post2/Pre)'); xlabel('Days');
        set(gca,'YLim',[0 3.1]);
        % ----------- Stats -------------
        % Comparing groups for each day:
        for d=1:length(days)
            currday = days(d);
            % Plot * if significant
            if h_dayripratefr_post2(currday)==1,
                mul = sign(mean(Confrday_riprate_post2(:,currday)));
                plot(currday, mean(Confrday_riprate_post2(:,currday)) + sem(Confrday_riprate_post2(:,currday))...
                    + 1.1*mul*sem(Confrday_riprate_post2(:,currday)), 'r*','MarkerSize',12);
            end
        end
        if savefig1==1,
            figfile = [figdir,'RippleRateRest_FracChgVsDay_PostSep'];
            print('-dpdf', figfile);
            print('-djpeg', figfile);
            saveas(gcf,figfile,'fig');
        end
        
    end % end epcomb- Frac change vs days
    
    %     %  Vs. Pairs of Sliding Days
    %     figure; hold on; redimscreen_figforppt1;
    %     errorbar(Conmeanfrslide_riprateep3,Conerrfrslide_riprateep3,'b--','LineWidth',1);
    %     errorbar(Conmeanfrslide_riprateep5,Conerrfrslide_riprateep5,'b','LineWidth',2);
    %     errorbar(Expmeanfrslide_riprateep3,Experrfrslide_riprateep3,'r--','LineWidth',1);
    %     errorbar(Expmeanfrslide_riprateep5,Experrfrslide_riprateep5,'r','LineWidth',2);
    %     if ntet==2
    %         set(gca,'YLim',[0 4.0]);
    %     else
    %         set(gca,'YLim',[0 2.2]);
    %     end
    %     set(gca,'XLim',[0.5 max(days)-0.5]);
    %     title('Increase in Ripple Rate during Post-Beh Sleep');
    %     ylabel('Fraction Change in Ripple Rate'); xlabel('Pairs of Sliding Days');
    
    
    
    
    % ********************** Still Time ****************************
    % ***************************************************************
    
    %-----------------------------------------------------------
    %-----------------------------------------------------------
    % 1) Bar plot for all days
    figure; hold on;
    if forppr==1
        redimscreen_figforppr1;
    else
        redimscreen_figforppt1;
    end
    
    % Epoch 2 and 4 combined
    %-----------------------
    if epcomb==1
        % Plot as bars
        bar(0.5,mean(Constilltime_ep1(:)),'b');
        bar(1.5,mean(Expstilltime_ep1(:)),'r');
        bar(3.5,mean([Constilltime_ep3(:); Constilltime_ep5(:)]),'b');
        bar(4.5,mean([Expstilltime_ep3(:); Expstilltime_ep5(:)]),'r');
        errorbar(0.5,mean(Constilltime_ep1(:)),std(Constilltime_ep1(:)),'k','LineWidth',2);
        errorbar(1.5,mean(Expstilltime_ep1(:)),std(Expstilltime_ep1(:)),'k','LineWidth',2);
        errorbar(3.5,mean([Constilltime_ep3(:); Constilltime_ep5(:)]),std([Constilltime_ep3(:); Constilltime_ep5(:)]),'k','LineWidth',2);
        errorbar(4.5,mean([Expstilltime_ep3(:); Expstilltime_ep5(:)]),std([Expstilltime_ep3(:); Expstilltime_ep5(:)]),'k','LineWidth',2);
        % Plot as groups
        %         errorbar([1,2], [mean(Constilltime_ep1(:));mean([Constilltime_ep3(:); Constilltime_ep5(:)])],...
        %             [std(Constilltime_ep1(:));sem([Constilltime_ep3(:); Constilltime_ep5(:)])],'bd-','MarkerSize',12,'LineWidth',3);
        %         errorbar([1,2], [mean(Expstilltime_ep1(:));mean([Expstilltime_ep3(:); Expstilltime_ep5(:)])],...
        %             [std(Expstilltime_ep1(:));sem([Expstilltime_ep3(:); Expstilltime_ep5(:)])],'ro-','MarkerSize',12,'LineWidth',3);
        % ----------- Stats -------------
        %Nanimals*Ndays  pts in each group: each day is a measure
        % Con against Exp for epoch1 and epochs3and5 combined
        % ---------------
        [h11_allstilltime,p11_allstilltime] = ttest2(Constilltime_ep1(:),Expstilltime_ep1(:));
        [h3535_allstilltime,p3535_allstilltime] = ttest2([Constilltime_ep3(:); Constilltime_ep5(:)],[Expstilltime_ep3(:); Expstilltime_ep5(:)]);
        % Use anova_rm/anova2/anovan/rm_anova2
        % anova_rm: pgrp=; ptime=0.1971; pint=0.0658
        [panovarm_allstilltime] = anova_rm({[Constilltime_ep1(:),mean([Constilltime_ep3(:),Constilltime_ep5(:)],2)] [Expstilltime_ep1(:),mean([Expstilltime_ep3(:),Expstilltime_ep5(:)],2)]},'off');
        % anova2 with epochs3and5 meaned. each day is a measure. factors = Grp and Pre-Post condition. 4*8=32 repeats
        % anova2: pgrp=,prows(pre-post)=0.3018,pint=0.139
%        [panova2_allstilltime] = anova2(  [[Constilltime_ep1(:);mean([Constilltime_ep3(:),Constilltime_ep5(:)],2)],[Expstilltime_ep1(:);mean([Expstilltime_ep3(:),Expstilltime_ep5(:)],2)]],32,'off');
        % rm_anova2 will treat both day and time as factors for different subjects. Wont treat day as independent factor
        title('StillTime during Rest: All Days');
        ylabel('StillTime (s)');
        set(gca,'XTick',[1,4],'XTickLabel',{'Pre';'Post'},'FontSize',xfont,'Fontweight','normal');
        set(gca,'YLim',[0 900]);
        %set(gca,'XLim',[0.5 3.5]);
        %axis([0 4 0 0.5])
        if savefig1==1,
            figfile = [figdir,'StillTimeRest_PostComb'];
            print('-dpdf', figfile);
            print('-djpeg', figfile);
            saveas(gcf,figfile,'fig');
        end
        
        
    else
        for ep=1:length(sleeps)
            currep=sleeps(ep);
            eval(['allmean_stilltime(',num2str(ep),',:) = [Conallmean_stilltimeep',num2str(currep),';  Expallmean_stilltimeep',num2str(currep),'];'])
        end
        
        % Plot as bars
        bar(allmean_stilltime,'grouped');
        errorbar(0.8,Conallmean_stilltimeep1,Conallerr_stilltimeep1,'k');
        errorbar(1.2,Expallmean_stilltimeep1,Expallerr_stilltimeep1,'k');
        errorbar(1.8,Conallmean_stilltimeep3,Conallerr_stilltimeep3,'k');
        errorbar(2.2,Expallmean_stilltimeep3,Expallerr_stilltimeep3,'k');
        errorbar(2.8,Conallmean_stilltimeep5,Conallerr_stilltimeep5,'k');
        errorbar(3.2,Expallmean_stilltimeep5,Expallerr_stilltimeep5,'k');
        % Plot as groups
        %         errorbar([1:3], [mean(Constilltime_ep1(:));mean(Constilltime_ep3(:)); mean(Constilltime_ep5(:))],...
        %             [sem(Constilltime_ep1(:));sem(Constilltime_ep3(:)); sem(Constilltime_ep5(:))],'bd-','MarkerSize',12,'LineWidth',3);
        %         errorbar([1:3], [mean(Expstilltime_ep1(:));mean(Expstilltime_ep3(:)); mean(Expstilltime_ep5(:))],...
        %             [sem(Expstilltime_ep1(:));sem(Expstilltime_ep3(:)); sem(Expstilltime_ep5(:))],'ro-','MarkerSize',12,'LineWidth',3);
        % ----------- Stats -------------
        % Comparing animals: each animal has one number: average across days
        % Con against Exp for each epoch
        [h11_allstilltime,p11_allstilltime] = ttest2(Constilltime_ep1(:),Expstilltime_ep1(:));
        [h33_allstilltime,p33_allstilltime] = ttest2(Constilltime_ep3(:),Expstilltime_ep3(:));
        [h55_allstilltime,p55_allstilltime] = ttest2(Constilltime_ep5(:),Expstilltime_ep5(:));
        % Use anova_rm/anova2/anovan/rm_anova2
        % anova_rm: pgrp=; ptime(pre-post1-post2)=0.2853; pint=0.1
        [panova_allstilltime] = anova_rm({[Constilltime_ep1(:),Constilltime_ep3(:),Constilltime_ep5(:)] [Expstilltime_ep1(:),Expstilltime_ep3(:),Expstilltime_ep5(:)]},'on');
        % anova2. each day is a measure. factors = Grp and Pre-Post condition. 4*8=32 repeats
        % anova2: pgrp=,prows(pre-post1-post2)=0.4517,pint=0.2307
        [panova2_allstilltime] = anova2([[Constilltime_ep1(:);Constilltime_ep3(:);Constilltime_ep5(:)],[Expstilltime_ep1(:);Expstilltime_ep3(:);Expstilltime_ep5(:)]],32,'on');
        % rm_anova2 will treat both day and time as factors for different subjects. Wont treat day as independent factor
        title('Still Time during Rest: All Days');
        ylabel('Still Time (s)');
        set(gca,'XTick',[1:3],'XTickLabel',{'Pre';'Post1';'Post2'},'FontSize',16,'Fontweight','normal');
        set(gca,'YLim',[0 900]);
        %set(gca,'XLim',[0.5 3.5]);
        %axis([0 4 0 0.5])
        if savefig1==1,
            figfile = [figdir,'StillTimeRest_PostSep'];
            print('-dpdf', figfile);
            print('-djpeg', figfile);
            saveas(gcf,figfile,'fig');
        end
        
    end % epcomb
    
    
    
    
    
    
end % end figopt1







% ***********************************************************************************
% ***********************************************************************************




% ***************  Plots For Ripple Size  ******************
% ***************************************************************

if figopt2==1
    
    ripsizeedgesn =[3,ripsizeedges];
    figoptsize = 1;
    
    % Get values, including combining ep3 and ep5
    %-------------------------
    % Ep1
    Conmeanhist1 = mean(squeeze(mean(Conripsizehist_dayanim_ep1,2))); % Mean along animals, then days
    Conerrhist1 = sem(squeeze(sem(Conripsizehist_dayanim_ep1,2)));
    Expmeanhist1 = mean(squeeze(mean(Expripsizehist_dayanim_ep1,2)));
    Experrhist1 = sem(squeeze(sem(Expripsizehist_dayanim_ep1,2)));
    % Ep3
    Conmeanhist3 = mean(squeeze(mean(Conripsizehist_dayanim_ep3,2))); % Mean along animals, then days
    Conerrhist3 = sem(squeeze(sem(Conripsizehist_dayanim_ep3,2)));
    Expmeanhist3 = mean(squeeze(mean(Expripsizehist_dayanim_ep3,2)));
    Experrhist3 = sem(squeeze(sem(Expripsizehist_dayanim_ep3,2)));
    % Ep5
    Conmeanhist5 = mean(squeeze(mean(Conripsizehist_dayanim_ep5,2)));
    Conerrhist5 = sem(squeeze(sem(Conripsizehist_dayanim_ep5,2)));
    Expmeanhist5 = mean(squeeze(mean(Expripsizehist_dayanim_ep5,2)));
    Experrhist5 = sem(squeeze(sem(Expripsizehist_dayanim_ep5,2)));
    % Combine post across epochs
    Conmeanhist = mean([Conmeanhist3;Conmeanhist5]);
    Conerrhist = sem([Conerrhist3;Conerrhist5]);
    Expmeanhist = mean([Expmeanhist3;Expmeanhist5]);
    Experrhist = sem([Experrhist3;Experrhist5]);
    % Vectors of histogram for current day - for stats
    Conhistvec_ep1_alldays = Conripsizehist_dayanim_ep1(:);
    Exphistvec_ep1_alldays = Expripsizehist_dayanim_ep1(:);
    % Vectors of histogram for current day - for stats
    Conhistvec_ep3_alldays = Conripsizehist_dayanim_ep3(:);
    Exphistvec_ep3_alldays = Expripsizehist_dayanim_ep3(:);
    Conhistvec_ep5_alldays = Conripsizehist_dayanim_ep5(:);
    Exphistvec_ep5_alldays = Expripsizehist_dayanim_ep5(:);
    
    
    % *********************
    % Con vs Exp
    % *********************
    if epcomb==1
        % Epoch 3 and 5 combined
        %-----------------------
        figure; hold on;
        redimscreen_2horsubplots;
        subplot(1,2,1); hold on;
        plot(ripsizeedgesn,[0,Conmeanhist1],'b.-','Linewidth',2,'MarkerSize',18);
        jbfill(ripsizeedges, [Conmeanhist1+Conerrhist1],...
            [Conmeanhist1-Conerrhist1],'b','b',1,0.3);
        plot(ripsizeedgesn,[0,Expmeanhist1],'r.-','Linewidth',2,'MarkerSize',18);
        jbfill(ripsizeedges, [Expmeanhist1+Experrhist1],...
            [Expmeanhist1-Experrhist1],'r','r',1,0.3);
        set(gca,'YLim',[0 1]); set(gca,'XLim',[0 max(ripsizeedges)+1]);
        title(['Pre: Con vs Exp - Ripple Size during Rest. All Days'],'FontSize',tfont,'Fontweight','normal');
        ylabel('Cumulative Proportion','FontSize',yfont,'Fontweight','normal');
        xlabel('Ripple Size (stdev)','FontSize',xfont,'Fontweight','normal');
        % ----------- Stats -------------
        % Raw values of ripsize in std combined across animals for all days
        [h_ripsize0,p_ripsize0] = ttest2(Conripsize_ep1_alldays,Expripsize_ep1_alldays);
        % ripsizehist in std combined across animals for current day
        [p_ripsizehist0,h_ripsizehist0] = ranksum(Conhistvec_ep1_alldays,Exphistvec_ep1_alldays);
        %         if h_ripsize0 == 1
        %             plot(7, 0.5, 'r*','MarkerSize',12);
        %         end
        %text(7,0.4,['p = ',num2str(p_ripsize0)],'FontSize',xfont);
        %---------------------------------
        subplot(1,2,2); hold on;
        plot(ripsizeedgesn,[0,Conmeanhist],'b.-','Linewidth',2,'MarkerSize',18);
        jbfill(ripsizeedges, [Conmeanhist+Conerrhist],...
            [Conmeanhist-Conerrhist],'b','b',1,0.3);
        plot(ripsizeedgesn,[0,Expmeanhist],'r.-','Linewidth',2,'MarkerSize',18);
        jbfill(ripsizeedges, [Expmeanhist+Experrhist],...
            [Expmeanhist-Experrhist],'r','r',1,0.3);
        set(gca,'YLim',[0 1]); set(gca,'XLim',[0 max(ripsizeedges)+1]);
        title(['Post: Con vs Exp - Ripple Size during Rest. All Days'],'FontSize',tfont,'Fontweight','normal');
        ylabel('Cumulative Proportion','FontSize',yfont,'Fontweight','normal');
        xlabel('Ripple Size (stdev)','FontSize',xfont,'Fontweight','normal');
        % ----------- Stats -------------
        % Raw values of ripsize in std combined across animals for all days
        [h_ripsize,p_ripsize] = ttest2([Conripsize_ep3_alldays, Conripsize_ep5_alldays],[Expripsize_ep3_alldays, Expripsize_ep5_alldays]);
        % ripsizehist in std combined across animals for current day
        [p_ripsizehist,h_ripsizehist] = ranksum([Conhistvec_ep3_alldays;Conhistvec_ep5_alldays],[Exphistvec_ep3_alldays;Exphistvec_ep5_alldays]);
        %         if h_ripsize == 1
        %             plot(7, 0.5, 'r*','MarkerSize',12);
        %         end
        %text(7,0.4,['p = ',num2str(p_ripsize)],'FontSize',xfont);
        if savefig1==1,
            figfile = [figdir,'RippleSizeRest_ConvsExp_PostComb'];
            print('-dpdf', figfile);
            print('-djpeg', figfile);
            saveas(gcf,figfile,'fig');
        end
        %---------------------------------
    else
        % Epoch 3 and 5 separate
        %-----------------------
        figure; hold on;
        redimscreen_2horsubplots;
        subplot(1,3,1); hold on;
        plot(ripsizeedgesn,[0,Conmeanhist1],'b.-','Linewidth',2,'MarkerSize',18);
        jbfill(ripsizeedges, [Conmeanhist1+Conerrhist1],...
            [Conmeanhist1-Conerrhist1],'b','b',1,0.3);
        plot(ripsizeedgesn,[0,Expmeanhist1],'r.-','Linewidth',2,'MarkerSize',18);
        jbfill(ripsizeedges, [Expmeanhist1+Experrhist1],...
            [Expmeanhist1-Experrhist1],'r','r',1,0.3);
        set(gca,'YLim',[0 1]); set(gca,'XLim',[0 max(ripsizeedges)+1]);
        title(['Pre: Con vs Exp - Ripple Size during Rest. All Days'],'FontSize',tfont,'Fontweight','normal');
        ylabel('Cumulative Proportion','FontSize',yfont,'Fontweight','normal');
        xlabel('Ripple Size (stdev)','FontSize',xfont,'Fontweight','normal');
        % ----------- Stats -------------
        % Raw values of ripsize in std combined across animals for all days
        [h_ripsize0,p_ripsize0] = ttest2(Conripsize_ep1_alldays,Expripsize_ep1_alldays);
        % ripsizehist in std combined across animals for current day
        [p_ripsizehist0,h_ripsizehist0] = ranksum(Conhistvec_ep1_alldays,Exphistvec_ep1_alldays);
        %         if h_ripsize0 == 1
        %             plot(7, 0.5, 'r*','MarkerSize',12);
        %         end
        %text(7,0.4,['p = ',num2str(p_ripsize0)],'FontSize',xfont);
        %---------------------------------
        subplot(1,3,2); hold on;
        plot(ripsizeedgesn,[0,Conmeanhist3],'b.-','Linewidth',2,'MarkerSize',18);
        jbfill(ripsizeedges, [Conmeanhist3+Conerrhist3],...
            [Conmeanhist3-Conerrhist3],'b','b',1,0.3);
        plot(ripsizeedgesn,[0,Expmeanhist3],'r.-','Linewidth',2,'MarkerSize',18);
        jbfill(ripsizeedges, [Expmeanhist3+Experrhist3],...
            [Expmeanhist3-Experrhist3],'r','r',1,0.3);
        set(gca,'YLim',[0 1]); set(gca,'XLim',[0 max(ripsizeedges)+1]);
        title(['Post1: Con vs Exp - Ripple Size during Rest. All Days'],'FontSize',tfont,'Fontweight','normal');
        ylabel('Cumulative Proportion','FontSize',yfont,'Fontweight','normal');
        xlabel('Ripple Size (stdev)','FontSize',xfont,'Fontweight','normal');
        % ----------- Stats -------------
        % Raw values of ripsize in std combined across animals for all days
        [h_ripsize1,p_ripsize1] = ttest2(Conripsize_ep3_alldays,Expripsize_ep3_alldays);
        % ripsizehist in std combined across animals for current day
        [p_ripsizehist1,h_ripsizehist1] = ranksum(Conhistvec_ep3_alldays,Exphistvec_ep3_alldays);
        %         if h_ripsize1 == 1
        %             plot(7, 0.5, 'r*','MarkerSize',12);
        %         end
        %text(7,0.4,['p = ',num2str(p_ripsize1)],'FontSize',xfont);
        %---------------------------------
        subplot(1,3,3); hold on;
        plot(ripsizeedgesn,[0,Conmeanhist5],'b.-','Linewidth',2,'MarkerSize',18);
        jbfill(ripsizeedges, [Conmeanhist5+Conerrhist5],...
            [Conmeanhist5-Conerrhist5],'b','b',1,0.3);
        plot(ripsizeedgesn,[0,Expmeanhist5],'r.-','Linewidth',2,'MarkerSize',18);
        jbfill(ripsizeedges, [Expmeanhist5+Experrhist5],...
            [Expmeanhist5-Experrhist5],'r','r',1,0.3);
        set(gca,'YLim',[0 1]); set(gca,'XLim',[0 max(ripsizeedges)+1]);
        title(['Post2: Con vs Exp - Ripple Size during Rest. All Days'],'FontSize',tfont,'Fontweight','normal');
        ylabel('Cumulative Proportion','FontSize',yfont,'Fontweight','normal');
        xlabel('Ripple Size (stdev)','FontSize',xfont,'Fontweight','normal');
        % ----------- Stats -------------
        % Raw values of ripsize in std combined across animals for all days
        [h_ripsize2,p_ripsize2] = ttest2(Conripsize_ep5_alldays,Expripsize_ep5_alldays);
        % ripsizehist in std combined across animals for current day
        [p_ripsizehist2,h_ripsizehist2] = ranksum(Conhistvec_ep5_alldays,Exphistvec_ep5_alldays);
        %         if h_ripsize2 == 1
        %             plot(7, 0.5, 'r*','MarkerSize',12);
        %         end
        %text(7,0.4,['p = ',num2str(p_ripsize2)],'FontSize',xfont);
        %---------------------------------
        if savefig1==1,
            figfile = [figdir,'RippleSizeRest_ConvsExp_PostSep'];
            print('-dpdf', figfile);
            print('-djpeg', figfile);
            saveas(gcf,figfile,'fig');
        end
    end % end epcomb
    
    
    % *****************************
    % Con vs Con   and   Exp vs Exp
    % *****************************
    if epcomb==1
        % Epoch 3 and 5 combined
        %-----------------------
        figure; hold on;
        redimscreen_2horsubplots;
        subplot(1,2,1); hold on;
        plot(ripsizeedgesn,[0,Conmeanhist1],'k.-','Linewidth',2,'MarkerSize',18);
        jbfill(ripsizeedges, [Conmeanhist1+Conerrhist1],...
            [Conmeanhist1-Conerrhist1],'k','k',1,0.3);
        plot(ripsizeedgesn,[0,Conmeanhist],'b.-','Linewidth',2,'MarkerSize',18);
        jbfill(ripsizeedges, [Conmeanhist+Conerrhist],...
            [Conmeanhist-Conerrhist],'b','b',1,0.3);
        set(gca,'YLim',[0 1]); set(gca,'XLim',[0 max(ripsizeedges)+1]);
        title(['Pre vs Post: Con - Ripple Size during Rest. All Days'],'FontSize',tfont,'Fontweight','normal');
        ylabel('Cumulative Proportion','FontSize',yfont,'Fontweight','normal');
        xlabel('Ripple Size (stdev)','FontSize',xfont,'Fontweight','normal');
        % ----------- Stats -------------
        % Raw values of ripsize in std combined across animals for all days
        [h_ripsize0,p_ripsize0] = ttest2(Conripsize_ep1_alldays,[Conripsize_ep3_alldays, Conripsize_ep5_alldays]);
        % ripsizehist in std combined across animals for current day
        [p_ripsizehist0,h_ripsizehist0] = ranksum(Conhistvec_ep1_alldays,[Conhistvec_ep3_alldays; Conhistvec_ep5_alldays]);
        %         if h_ripsize0 == 1
        %             plot(7, 0.5, 'r*','MarkerSize',12);
        %         end
        %text(7,0.4,['p = ',num2str(p_ripsize0)],'FontSize',xfont);
        %---------------------------------
        subplot(1,2,2); hold on;
        plot(ripsizeedgesn,[0,Expmeanhist1],'k.-','Linewidth',2,'MarkerSize',18);
        jbfill(ripsizeedges, [Expmeanhist1+Experrhist1],...
            [Expmeanhist1-Experrhist1],'k','k',1,0.3);
        plot(ripsizeedgesn,[0,Expmeanhist],'r.-','Linewidth',2,'MarkerSize',18);
        jbfill(ripsizeedges, [Expmeanhist+Experrhist],...
            [Expmeanhist-Experrhist],'r','r',1,0.3);
        set(gca,'YLim',[0 1]); set(gca,'XLim',[0 max(ripsizeedges)+1]);
        title(['Pre vs Post: Exp - Ripple Size during Rest. All Days'],'FontSize',tfont,'Fontweight','normal');
        ylabel('Cumulative Proportion','FontSize',yfont,'Fontweight','normal');
        xlabel('Ripple Size (stdev)','FontSize',xfont,'Fontweight','normal');
        % ----------- Stats -------------
        % Raw values of ripsize in std combined across animals for all days
        [h_ripsize,p_ripsize] = ttest2(Expripsize_ep1_alldays,[Expripsize_ep3_alldays, Expripsize_ep5_alldays]);
        % ripsizehist in std combined across animals for current day
        [p_ripsizehist,h_ripsizehist] = ranksum(Exphistvec_ep1_alldays,[Exphistvec_ep3_alldays;Exphistvec_ep5_alldays]);
        %         if h_ripsize == 1
        %             plot(7, 0.5, 'r*','MarkerSize',12);
        %         end
        %text(7,0.4,['p = ',num2str(p_ripsize)],'FontSize',xfont);
        %---------------------------------
        if savefig1==1,
            figfile = [figdir,'RippleSizeRest_PrevsPost_PostComb'];
            print('-dpdf', figfile);
            print('-djpeg', figfile);
            saveas(gcf,figfile,'fig');
        end
    else
        % Epoch 3 and 5 separate
        %-----------------------
        figure; hold on;
        redimscreen_2horsubplots;
        subplot(2,2,1); hold on;
        plot(ripsizeedgesn,[0,Conmeanhist1],'k.-','Linewidth',2,'MarkerSize',18);
        jbfill(ripsizeedges, [Conmeanhist1+Conerrhist1],...
            [Conmeanhist1-Conerrhist1],'k','k',1,0.3);
        plot(ripsizeedgesn,[0,Conmeanhist3],'b.-','Linewidth',2,'MarkerSize',18);
        jbfill(ripsizeedges, [Conmeanhist3+Conerrhist3],...
            [Conmeanhist3-Conerrhist3],'b','b',1,0.3);
        set(gca,'YLim',[0 1]); set(gca,'XLim',[0 max(ripsizeedges)+1]);
        title(['Pre vs Post1: Con  - Ripple Size during Rest. All Days'],'FontSize',tfont,'Fontweight','normal');
        ylabel('Cumulative Proportion','FontSize',yfont,'Fontweight','normal');
        xlabel('Ripple Size (stdev)','FontSize',xfont,'Fontweight','normal');
        % ----------- Stats -------------
        % Raw values of ripsize in std combined across animals for all days
        [h_ripsize0,p_ripsize0] = ttest2(Conripsize_ep1_alldays,Conripsize_ep3_alldays);
        % ripsizehist in std combined across animals for current day
        [p_ripsizehist0,h_ripsizehist0] = ranksum(Conhistvec_ep1_alldays,Conhistvec_ep3_alldays);
        %         if h_ripsize0 == 1
        %             plot(7, 0.5, 'r*','MarkerSize',12);
        %         end
        %text(7,0.4,['p = ',num2str(p_ripsize0)],'FontSize',xfont);
        %---------------------------------
        subplot(2,2,2); hold on;
        plot(ripsizeedgesn,[0,Conmeanhist1],'k.-','Linewidth',2,'MarkerSize',18);
        jbfill(ripsizeedges, [Conmeanhist1+Conerrhist1],...
            [Conmeanhist1-Conerrhist1],'k','k',1,0.3);
        plot(ripsizeedgesn,[0,Conmeanhist5],'b.-','Linewidth',2,'MarkerSize',18);
        jbfill(ripsizeedges, [Conmeanhist5+Conerrhist5],...
            [Conmeanhist5-Conerrhist5],'b','b',1,0.3);
        set(gca,'YLim',[0 1]); set(gca,'XLim',[0 max(ripsizeedges)+1]);
        title(['Pre vs Post2: Con - Ripple Size during Rest. All Days'],'FontSize',tfont,'Fontweight','normal');
        ylabel('Cumulative Proportion','FontSize',yfont,'Fontweight','normal');
        xlabel('Ripple Size (stdev)','FontSize',xfont,'Fontweight','normal');
        % ----------- Stats -------------
        % Raw values of ripsize in std combined across animals for all days
        [h_ripsize1,p_ripsize1] = ttest2(Conripsize_ep1_alldays,Conripsize_ep5_alldays);
        % ripsizehist in std combined across animals for current day
        [p_ripsizehist1,h_ripsizehist1] = ranksum(Conhistvec_ep1_alldays,Conhistvec_ep5_alldays);
        %         if h_ripsize1 == 1
        %             plot(7, 0.5, 'r*','MarkerSize',12);
        %         end
        %text(7,0.4,['p = ',num2str(p_ripsize1)],'FontSize',xfont);
        %---------------------------------
        subplot(2,2,3); hold on;
        plot(ripsizeedgesn,[0,Expmeanhist1],'k.-','Linewidth',2,'MarkerSize',18);
        jbfill(ripsizeedges, [Expmeanhist1+Experrhist1],...
            [Expmeanhist1-Experrhist1],'k','k',1,0.3);
        plot(ripsizeedgesn,[0,Expmeanhist3],'r.-','Linewidth',2,'MarkerSize',18);
        jbfill(ripsizeedges, [Expmeanhist3+Experrhist3],...
            [Expmeanhist3-Experrhist3],'r','r',1,0.3);
        set(gca,'YLim',[0 1]); set(gca,'XLim',[0 max(ripsizeedges)+1]);
        title(['Pre vs Post1: Exp - Ripple Size during Rest. All Days'],'FontSize',tfont,'Fontweight','normal');
        ylabel('Cumulative Proportion','FontSize',yfont,'Fontweight','normal');
        xlabel('Ripple Size (stdev)','FontSize',xfont,'Fontweight','normal');
        % ----------- Stats -------------
        % Raw values of ripsize in std combined across animals for all days
        [h_ripsize2,p_ripsize2] = ttest2(Expripsize_ep1_alldays,Expripsize_ep3_alldays);
        % ripsizehist in std combined across animals for current day
        [p_ripsizehist2,h_ripsizehist2] = ranksum(Exphistvec_ep1_alldays,Exphistvec_ep3_alldays);
        %         if h_ripsize2 == 1
        %             plot(7, 0.5, 'r*','MarkerSize',12);
        %         end
        %         text(7,0.4,['p = ',num2str(p_ripsize2)],'FontSize',xfont);
        %---------------------------------
        subplot(2,2,4); hold on;
        plot(ripsizeedgesn,[0,Expmeanhist1],'k.-','Linewidth',2,'MarkerSize',18);
        jbfill(ripsizeedges, [Expmeanhist1+Experrhist1],...
            [Expmeanhist1-Experrhist1],'k','k',1,0.3);
        plot(ripsizeedgesn,[0,Expmeanhist5],'r.-','Linewidth',2,'MarkerSize',18);
        jbfill(ripsizeedges, [Expmeanhist5+Experrhist5],...
            [Expmeanhist5-Experrhist5],'r','r',1,0.3);
        set(gca,'YLim',[0 1]); set(gca,'XLim',[0 max(ripsizeedges)+1]);
        title(['Pre vs Post2: Exp - Ripple Size during Rest. All Days'],'FontSize',tfont,'Fontweight','normal');
        ylabel('Cumulative Proportion','FontSize',yfont,'Fontweight','normal');
        xlabel('Ripple Size (stdev)','FontSize',xfont,'Fontweight','normal');
        % ----------- Stats -------------
        % Raw values of ripsize in std combined across animals for all days
        [h_ripsize3,p_ripsize3] = ttest2(Expripsize_ep1_alldays,Expripsize_ep5_alldays);
        % ripsizehist in std combined across animals for current day
        [p_ripsizehist3,h_ripsizehist3] = ranksum(Exphistvec_ep1_alldays,Exphistvec_ep5_alldays);
        %         if h_ripsize3 == 1
        %             plot(7, 0.5, 'r*','MarkerSize',12);
        %         end
        %text(7,0.4,['p = ',num2str(p_ripsize2)],'FontSize',xfont);
        %---------------------------------
        if savefig1==1,
            figfile = [figdir,'RippleSizeRest_PrevsPost_PostSep'];
            print('-dpdf', figfile);
            print('-djpeg', figfile);
            saveas(gcf,figfile,'fig');
        end
    end % end epcomb
    
    
    % **************************
    % For each day separately
    % **************************
    if dodayripsize==1
        %-------------------------
        days=1:size(Conripsizehist_dayanim_ep3,1);
        %days=[1,2,3,5,8];
        for n = 1:length(days)
            currday = days(n);
            % Ep1
            Conmeanhist1 = squeeze(mean(Conripsizehist_dayanim_ep1(currday,:,:),2))'; % Mean along 2nd dimension of animals
            Conerrhist1 = squeeze(sem(Conripsizehist_dayanim_ep1(currday,:,:),2))'; % sem along 2nd dimension of animals
            Expmeanhist1 = squeeze(mean(Expripsizehist_dayanim_ep1(currday,:,:),2))';
            Experrhist1 = squeeze(sem(Expripsizehist_dayanim_ep1(currday,:,:),2))';
            % Ep3
            Conmeanhist3 = squeeze(mean(Conripsizehist_dayanim_ep3(currday,:,:),2))'; % Mean along 2nd dimension of animals
            Conerrhist3 = squeeze(sem(Conripsizehist_dayanim_ep3(currday,:,:),2))'; % sem along 2nd dimension of animals
            Expmeanhist3 = squeeze(mean(Expripsizehist_dayanim_ep3(currday,:,:),2))';
            Experrhist3 = squeeze(sem(Expripsizehist_dayanim_ep3(currday,:,:),2))';
            % Ep5
            Conmeanhist5 = squeeze(mean(Conripsizehist_dayanim_ep5(currday,:,:),2))'; % Mean along 2nd dimension of animals
            Conerrhist5 = squeeze(sem(Conripsizehist_dayanim_ep5(currday,:,:),2))'; % sem along 2nd dimension of animals
            Expmeanhist5 = squeeze(mean(Expripsizehist_dayanim_ep5(currday,:,:),2))';
            Experrhist5 = squeeze(sem(Expripsizehist_dayanim_ep5(currday,:,:),2))';
            % Combine across epochs3and5
            Conmeanhist = mean([Conmeanhist3;Conmeanhist5],1);
            Conerrhist = sem([Conerrhist3;Conerrhist5],1);
            Expmeanhist = mean([Expmeanhist3;Expmeanhist5],1);
            Experrhist = sem([Experrhist3;Experrhist5],1);
            % Vectors of histogram for current day - for stats
            Conhistvec_ep1_currday = squeeze(Conripsizehist_dayanim_ep1(currday,:,:)); Conhistvec_ep1(currday,:) = Conhistvec_ep1_currday(:);
            Exphistvec_ep1_currday = squeeze(Expripsizehist_dayanim_ep1(currday,:,:)); Exphistvec_ep1(currday,:) = Exphistvec_ep1_currday(:);
            Conhistvec_ep3_currday = squeeze(Conripsizehist_dayanim_ep3(currday,:,:)); Conhistvec_ep3(currday,:) = Conhistvec_ep3_currday(:);
            Exphistvec_ep3_currday = squeeze(Expripsizehist_dayanim_ep3(currday,:,:)); Exphistvec_ep3(currday,:) = Exphistvec_ep3_currday(:);
            Conhistvec_ep5_currday = squeeze(Conripsizehist_dayanim_ep5(currday,:,:)); Conhistvec_ep5(currday,:) = Conhistvec_ep5_currday(:);
            Exphistvec_ep5_currday = squeeze(Expripsizehist_dayanim_ep5(currday,:,:)); Exphistvec_ep5(currday,:) = Exphistvec_ep5_currday(:);
            
            % *********************
            % Con vs Exp
            % *********************
            if epcomb==1
                % Epoch 3 and 5 combined
                %-----------------------
                figure; hold on;
                redimscreen_2horsubplots;
                subplot(1,2,1); hold on;
                plot(ripsizeedgesn,[0,Conmeanhist1],'b.-','Linewidth',2,'MarkerSize',18);
                jbfill(ripsizeedges, [Conmeanhist1+Conerrhist1],...
                    [Conmeanhist1-Conerrhist1],'b','b',1,0.3);
                plot(ripsizeedgesn,[0,Expmeanhist1],'r.-','Linewidth',2,'MarkerSize',18);
                jbfill(ripsizeedges, [Expmeanhist1+Experrhist1],...
                    [Expmeanhist1-Experrhist1],'r','r',1,0.3);
                set(gca,'YLim',[0 1]); set(gca,'XLim',[0 max(ripsizeedges)+1]);
                title(['Pre: Con vs Exp - Ripple Size during Rest. Day' num2str(n)],'FontSize',tfont,'Fontweight','normal');
                ylabel('Cumulative Proportion','FontSize',yfont,'Fontweight','normal');
                xlabel('Ripple Size (stdev)','FontSize',xfont,'Fontweight','normal');
                % ----------- Stats -------------
                % Raw values of ripsize in std combined across animals for all days
                [h_ripsize0,p_ripsize0] = ttest2(Conripsize_ep1_day{currday},Expripsize_ep1_day{currday});
                % ripsizehist in std combined across animals for current day
                [p_ripsizehist0,h_ripsizehist0] = ranksum(Conhistvec_ep1(currday,:),Exphistvec_ep1(currday,:));
                %         if h_ripsize0 == 1
                %             plot(7, 0.5, 'r*','MarkerSize',12);
                %         end
                text(7,0.4,['p = ',num2str(p_ripsize0)],'FontSize',xfont);
                %---------------------------------
                subplot(1,2,2); hold on;
                plot(ripsizeedgesn,[0,Conmeanhist],'b.-','Linewidth',2,'MarkerSize',18);
                jbfill(ripsizeedges, [Conmeanhist+Conerrhist],...
                    [Conmeanhist-Conerrhist],'b','b',1,0.3);
                plot(ripsizeedgesn,[0,Expmeanhist],'r.-','Linewidth',2,'MarkerSize',18);
                jbfill(ripsizeedges, [Expmeanhist+Experrhist],...
                    [Expmeanhist-Experrhist],'r','r',1,0.3);
                set(gca,'YLim',[0 1]); set(gca,'XLim',[0 max(ripsizeedges)+1]);
                title(['Post: Con vs Exp - Ripple Size during Rest. Day' num2str(n)],'FontSize',tfont,'Fontweight','normal');
                ylabel('Cumulative Proportion','FontSize',yfont,'Fontweight','normal');
                xlabel('Ripple Size (stdev)','FontSize',xfont,'Fontweight','normal');
                % ----------- Stats -------------
                % Raw values of ripsize in std combined across animals for all days
                [h_ripsize,p_ripsize] = ttest2([Conripsize_ep3_day{currday}, Conripsize_ep5_day{currday}],[Expripsize_ep3_day{currday}, Expripsize_ep5_day{currday}]);
                % ripsizehist in std combined across animals for current day
                [p_ripsizehist,h_ripsizehist] = ranksum([Conhistvec_ep3(currday,:),Conhistvec_ep5(currday,:)],[Exphistvec_ep3(currday,:),Exphistvec_ep5(currday,:)]);
                if h_ripsize == 1
                    plot(7, 0.5, 'r*','MarkerSize',12);
                end
                text(7,0.4,['p = ',num2str(p_ripsize)],'FontSize',xfont);
                %---------------------------------
            else
                % Epoch 3 and 5 separate
                %-----------------------
                figure; hold on;
                redimscreen_2horsubplots;
                subplot(1,3,1); hold on;
                plot(ripsizeedgesn,[0,Conmeanhist1],'b.-','Linewidth',2,'MarkerSize',18);
                jbfill(ripsizeedges, [Conmeanhist1+Conerrhist1],...
                    [Conmeanhist1-Conerrhist1],'b','b',1,0.3);
                plot(ripsizeedgesn,[0,Expmeanhist1],'r.-','Linewidth',2,'MarkerSize',18);
                jbfill(ripsizeedges, [Expmeanhist1+Experrhist1],...
                    [Expmeanhist1-Experrhist1],'r','r',1,0.3);
                set(gca,'YLim',[0 1]); set(gca,'XLim',[0 max(ripsizeedges)+1]);
                title(['Pre: Con vs Exp - Ripple Size during Rest. Day' num2str(n)],'FontSize',tfont,'Fontweight','normal');
                ylabel('Cumulative Proportion','FontSize',yfont,'Fontweight','normal');
                xlabel('Ripple Size (stdev)','FontSize',xfont,'Fontweight','normal');
                % ----------- Stats -------------
                % Raw values of ripsize in std combined across animals for all days
                [h_ripsize0,p_ripsize0] = ttest2(Conripsize_ep1_day{currday},Expripsize_ep1_day{currday});
                % ripsizehist in std combined across animals for current day
                [p_ripsizehist0,h_ripsizehist0] = ranksum(Conhistvec_ep1(currday,:),Exphistvec_ep1(currday,:));
                %         if h_ripsize0 == 1
                %             plot(7, 0.5, 'r*','MarkerSize',12);
                %         end
                text(7,0.4,['p = ',num2str(p_ripsize0)],'FontSize',xfont);
                %---------------------------------
                subplot(1,3,2); hold on;
                plot(ripsizeedgesn,[0,Conmeanhist3],'b.-','Linewidth',2,'MarkerSize',18);
                jbfill(ripsizeedges, [Conmeanhist3+Conerrhist3],...
                    [Conmeanhist3-Conerrhist3],'b','b',1,0.3);
                plot(ripsizeedgesn,[0,Expmeanhist3],'r.-','Linewidth',2,'MarkerSize',18);
                jbfill(ripsizeedges, [Expmeanhist3+Experrhist3],...
                    [Expmeanhist3-Experrhist3],'r','r',1,0.3);
                set(gca,'YLim',[0 1]); set(gca,'XLim',[0 max(ripsizeedges)+1]);
                title(['Post1: Con vs Exp - Ripple Size during Rest. Day' num2str(n)],'FontSize',tfont,'Fontweight','normal');
                ylabel('Cumulative Proportion','FontSize',yfont,'Fontweight','normal');
                xlabel('Ripple Size (stdev)','FontSize',xfont,'Fontweight','normal');
                % ----------- Stats -------------
                % Raw values of ripsize in std combined across animals for all days
                [h_ripsize1,p_ripsize1] = ttest2(Conripsize_ep3_day{currday},Expripsize_ep3_day{currday});
                % ripsizehist in std combined across animals for current day
                [p_ripsizehist1,h_ripsizehist1] = ranksum(Conhistvec_ep3(currday,:),Exphistvec_ep3(currday,:));
                %         if h_ripsize1 == 1
                %             plot(7, 0.5, 'r*','MarkerSize',12);
                %         end
                text(7,0.4,['p = ',num2str(p_ripsize1)],'FontSize',xfont);
                %---------------------------------
                subplot(1,3,3); hold on;
                plot(ripsizeedgesn,[0,Conmeanhist5],'b.-','Linewidth',2,'MarkerSize',18);
                jbfill(ripsizeedges, [Conmeanhist5+Conerrhist5],...
                    [Conmeanhist5-Conerrhist5],'b','b',1,0.3);
                plot(ripsizeedgesn,[0,Expmeanhist5],'r.-','Linewidth',2,'MarkerSize',18);
                jbfill(ripsizeedges, [Expmeanhist5+Experrhist5],...
                    [Expmeanhist5-Experrhist5],'r','r',1,0.3);
                set(gca,'YLim',[0 1]); set(gca,'XLim',[0 max(ripsizeedges)+1]);
                title(['Post2: Con vs Exp - Ripple Size during Rest. Day' num2str(n)],'FontSize',tfont,'Fontweight','normal');
                ylabel('Cumulative Proportion','FontSize',yfont,'Fontweight','normal');
                xlabel('Ripple Size (stdev)','FontSize',xfont,'Fontweight','normal');
                % ----------- Stats -------------
                % Raw values of ripsize in std combined across animals for all days
                [h_ripsize2,p_ripsize2] = ttest2(Conripsize_ep5_day{currday},Expripsize_ep5_day{currday});
                % ripsizehist in std combined across animals for current day
                [p_ripsizehist2,h_ripsizehist2] = ranksum(Conhistvec_ep5(currday,:),Exphistvec_ep5(currday,:));
                %         if h_ripsize2 == 1
                %             plot(7, 0.5, 'r*','MarkerSize',12);
                %         end
                text(7,0.4,['p = ',num2str(p_ripsize2)],'FontSize',xfont);
                %---------------------------------
            end % if epcomb
            
            % *****************************
            % Con vs Con   and   Exp vs Exp
            % *****************************
            if epcomb==1
                % Epoch 3 and 5 combined
                %-----------------------
                figure; hold on;
                redimscreen_2horsubplots;
                subplot(1,2,1); hold on;
                plot(ripsizeedgesn,[0,Conmeanhist1],'k.-','Linewidth',2,'MarkerSize',18);
                jbfill(ripsizeedges, [Conmeanhist1+Conerrhist1],...
                    [Conmeanhist1-Conerrhist1],'k','k',1,0.3);
                plot(ripsizeedgesn,[0,Conmeanhist],'b.-','Linewidth',2,'MarkerSize',18);
                jbfill(ripsizeedges, [Conmeanhist+Conerrhist],...
                    [Conmeanhist-Conerrhist],'b','b',1,0.3);
                set(gca,'YLim',[0 1]); set(gca,'XLim',[0 max(ripsizeedges)+1]);
                title(['Pre vs Post: Con - Ripple Size during Rest. Day' num2str(n)],'FontSize',tfont,'Fontweight','normal');
                ylabel('Cumulative Proportion','FontSize',yfont,'Fontweight','normal');
                xlabel('Ripple Size (stdev)','FontSize',xfont,'Fontweight','normal');
                % ----------- Stats -------------
                % Raw values of ripsize in std combined across animals for all days
                [h_ripsize0,p_ripsize0] = ttest2(Conripsize_ep1_day{currday},[Conripsize_ep3_day{currday}, Conripsize_ep5_day{currday}]);
                % ripsizehist in std combined across animals for current day
                [p_ripsizehist0,h_ripsizehist0] = ranksum(Conhistvec_ep1(currday,:),[Conhistvec_ep3(currday,:),Conhistvec_ep5(currday,:)]);
                %         if h_ripsize0 == 1
                %             plot(7, 0.5, 'r*','MarkerSize',12);
                %         end
                text(7,0.4,['p = ',num2str(p_ripsize0)],'FontSize',xfont);
                %---------------------------------
                subplot(1,2,2); hold on;
                plot(ripsizeedgesn,[0,Expmeanhist1],'k.-','Linewidth',2,'MarkerSize',18);
                jbfill(ripsizeedges, [Expmeanhist1+Experrhist1],...
                    [Expmeanhist1-Experrhist1],'k','k',1,0.3);
                plot(ripsizeedgesn,[0,Expmeanhist],'r.-','Linewidth',2,'MarkerSize',18);
                jbfill(ripsizeedges, [Expmeanhist+Experrhist],...
                    [Expmeanhist-Experrhist],'r','r',1,0.3);
                set(gca,'YLim',[0 1]); set(gca,'XLim',[0 max(ripsizeedges)+1]);
                title(['Pre vs Post: Exp - Ripple Size during Rest. Day' num2str(n)],'FontSize',tfont,'Fontweight','normal');
                ylabel('Cumulative Proportion','FontSize',yfont,'Fontweight','normal');
                xlabel('Ripple Size (stdev)','FontSize',xfont,'Fontweight','normal');
                % ----------- Stats -------------
                % Raw values of ripsize in std combined across animals for all days
                [h_ripsize,p_ripsize] = ttest2(Expripsize_ep1_day{currday},[Expripsize_ep3_day{currday}, Expripsize_ep5_day{currday}]);
                % ripsizehist in std combined across animals for current day
                [p_ripsizehist,h_ripsizehist] = ranksum(Exphistvec_ep1(currday,:),[Exphistvec_ep3(currday,:),Exphistvec_ep5(currday,:)]);
                %         if h_ripsize == 1
                %             plot(7, 0.5, 'r*','MarkerSize',12);
                %         end
                text(7,0.4,['p = ',num2str(p_ripsize)],'FontSize',xfont);
                %---------------------------------
            else
                % Epoch 3 and 5 separate
                %-----------------------
                figure; hold on;
                redimscreen_2horsubplots;
                subplot(2,2,1); hold on;
                plot(ripsizeedgesn,[0,Conmeanhist1],'k.-','Linewidth',2,'MarkerSize',18);
                jbfill(ripsizeedges, [Conmeanhist1+Conerrhist1],...
                    [Conmeanhist1-Conerrhist1],'k','k',1,0.3);
                plot(ripsizeedgesn,[0,Conmeanhist3],'b.-','Linewidth',2,'MarkerSize',18);
                jbfill(ripsizeedges, [Conmeanhist3+Conerrhist3],...
                    [Conmeanhist3-Conerrhist3],'b','b',1,0.3);
                set(gca,'YLim',[0 1]); set(gca,'XLim',[0 max(ripsizeedges)+1]);
                title(['Pre vs Post1: Con  - Ripple Size during Rest. All Days'],'FontSize',tfont,'Fontweight','normal');
                ylabel('Cumulative Proportion','FontSize',yfont,'Fontweight','normal');
                xlabel('Ripple Size (stdev)','FontSize',xfont,'Fontweight','normal');
                % ----------- Stats -------------
                % Raw values of ripsize in std combined across animals for all days
                [h_ripsize0,p_ripsize0] = ttest2(Conripsize_ep1_day{currday},Conripsize_ep3_day{currday});
                % ripsizehist in std combined across animals for current day
                [p_ripsizehist0,h_ripsizehist0] = ranksum(Conhistvec_ep1(currday,:),Conhistvec_ep3(currday,:));
                %         if h_ripsize0 == 1
                %             plot(7, 0.5, 'r*','MarkerSize',12);
                %         end
                text(7,0.4,['p = ',num2str(p_ripsize0)],'FontSize',xfont);
                %---------------------------------
                subplot(2,2,2); hold on;
                plot(ripsizeedgesn,[0,Conmeanhist1],'k.-','Linewidth',2,'MarkerSize',18);
                jbfill(ripsizeedges, [Conmeanhist1+Conerrhist1],...
                    [Conmeanhist1-Conerrhist1],'k','k',1,0.3);
                plot(ripsizeedgesn,[0,Conmeanhist5],'b.-','Linewidth',2,'MarkerSize',18);
                jbfill(ripsizeedges, [Conmeanhist5+Conerrhist5],...
                    [Conmeanhist5-Conerrhist5],'b','b',1,0.3);
                set(gca,'YLim',[0 1]); set(gca,'XLim',[0 max(ripsizeedges)+1]);
                title(['Pre vs Post2: Con - Ripple Size during Rest. Day' num2str(n)],'FontSize',tfont,'Fontweight','normal');
                ylabel('Cumulative Proportion','FontSize',yfont,'Fontweight','normal');
                xlabel('Ripple Size (stdev)','FontSize',xfont,'Fontweight','normal');
                % ----------- Stats -------------
                % Raw values of ripsize in std combined across animals for all days
                [h_ripsize1,p_ripsize1] = ttest2(Conripsize_ep1_day{currday},Conripsize_ep5_day{currday});
                % ripsizehist in std combined across animals for current day
                [p_ripsizehist1,h_ripsizehist1] = ranksum(Conhistvec_ep1(currday,:),Conhistvec_ep5(currday,:));
                %         if h_ripsize1 == 1
                %             plot(7, 0.5, 'r*','MarkerSize',12);
                %         end
                text(7,0.4,['p = ',num2str(p_ripsize1)],'FontSize',xfont);
                %---------------------------------
                subplot(2,2,3); hold on;
                plot(ripsizeedgesn,[0,Expmeanhist1],'k.-','Linewidth',2,'MarkerSize',18);
                jbfill(ripsizeedges, [Expmeanhist1+Experrhist1],...
                    [Expmeanhist1-Experrhist1],'k','k',1,0.3);
                plot(ripsizeedgesn,[0,Expmeanhist3],'r.-','Linewidth',2,'MarkerSize',18);
                jbfill(ripsizeedges, [Expmeanhist3+Experrhist3],...
                    [Expmeanhist3-Experrhist3],'r','r',1,0.3);
                set(gca,'YLim',[0 1]); set(gca,'XLim',[0 max(ripsizeedges)+1]);
                title(['Pre vs Post1: Exp - Ripple Size during Rest. All Days'],'FontSize',tfont,'Fontweight','normal');
                ylabel('Cumulative Proportion','FontSize',yfont,'Fontweight','normal');
                xlabel('Ripple Size (stdev)','FontSize',xfont,'Fontweight','normal');
                % ----------- Stats -------------
                % Raw values of ripsize in std combined across animals for all days
                [h_ripsize2,p_ripsize2] = ttest2(Expripsize_ep1_day{currday},Expripsize_ep3_day{currday});
                % ripsizehist in std combined across animals for current day
                [p_ripsizehist2,h_ripsizehist2] = ranksum(Exphistvec_ep1(currday,:),Exphistvec_ep3(currday,:));
                %         if h_ripsize2 == 1
                %             plot(7, 0.5, 'r*','MarkerSize',12);
                %         end
                text(7,0.4,['p = ',num2str(p_ripsize2)],'FontSize',xfont);
                %---------------------------------
                subplot(2,2,4); hold on;
                plot(ripsizeedgesn,[0,Expmeanhist1],'k.-','Linewidth',2,'MarkerSize',18);
                jbfill(ripsizeedges, [Expmeanhist1+Experrhist1],...
                    [Expmeanhist1-Experrhist1],'k','k',1,0.3);
                plot(ripsizeedgesn,[0,Expmeanhist5],'r.-','Linewidth',2,'MarkerSize',18);
                jbfill(ripsizeedges, [Expmeanhist5+Experrhist5],...
                    [Expmeanhist5-Experrhist5],'r','r',1,0.3);
                set(gca,'YLim',[0 1]); set(gca,'XLim',[0 max(ripsizeedges)+1]);
                title(['Pre vs Post2: Exp - Ripple Size during Rest. Day' num2str(n)],'FontSize',tfont,'Fontweight','normal');
                ylabel('Cumulative Proportion','FontSize',yfont,'Fontweight','normal');
                xlabel('Ripple Size (stdev)','FontSize',xfont,'Fontweight','normal');
                % ----------- Stats -------------
                % Raw values of ripsize in std combined across animals for all days
                [h_ripsize3,p_ripsize3] = ttest2(Expripsize_ep1_day{currday},Expripsize_ep5_day{currday});
                % ripsizehist in std combined across animals for current day
                [p_ripsizehist3,h_ripsizehist3] = ranksum(Exphistvec_ep1(currday,:),Exphistvec_ep5(currday,:));
                %         if h_ripsize3 == 1
                %             plot(7, 0.5, 'r*','MarkerSize',12);
                %         end
                text(7,0.4,['p = ',num2str(p_ripsize2)],'FontSize',xfont);
                %---------------------------------
            end % end epcomb
            
            
        end % end days
    end % if dodayripsize==1
    
    keyboard;
    
    
end% if figopt2




























