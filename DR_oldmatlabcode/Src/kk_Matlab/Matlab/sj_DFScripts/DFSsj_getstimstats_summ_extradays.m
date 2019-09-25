
% Mainly for extra days - gives speed bef-aft stim ?


% From sj_stimstatssumm1. I want to broaden it to easily plot distributions, 
% eg. speed at stimuln distribution, etc

% sj_stimstatssumm1 and this script gets data saved by the analysis function sj_stim_behstats.m 
% which saves a file called prefix_stimstats in each "prefixdirect" directory
% Instead of loading the file,
% Run the analysis from sj_stim_behstats again here, as DFAsj_stim_behstats
% Analysis function would load "pos","linpos" and "DIO" structure, and would use the iterator "epoch_behaveanal" 

clear; %close all;
runscript = 0; % if is 0, just loads saved processed data and plots stuff 
savedata = 0; % save data option - only works if runscript is also on
figopt1=0; % Figure Options - main summ figs 
dodistr=0; % Figures for speed distributions
dosepdays=0; % Separate figures for speed distributions on each day
dostimisi=0; % Figures for stimISI
doisisepdays=0;% Separate figures for stim isi on each day
doposdistr=0; % Figures for pos distr at stim - all days - If combining posn data across animals
dopossepdays=0; % Figures for position distr on each day - Can combine across animals for each day, or use to plot for 1 animal

dostimlearncorr=0; % Stimrate vs learning rate

doposquant=0; % Position at stimulation quantification - after review 

dovelstim_befaft=1; % Vel before after stimln - after review

% For binning and smoothing position
binsize = 2; % cm square 
stdev = 2; 
timestep = 1/29.97;
threshocc = 0.02; % Threshold occupancy in seconds

savedir = '/data25/sjadhav/RippleInterruption/ProcessedData/';
%savefile = [savedir 'Stimstats_summ'];
%savefile = [savedir 'Stimstats_summ_rev']; % after revision

savefile = [savedir 'Stimstats_summ_extradays']; % Days 9 and 10

% If runscript, run Datafilter and save data
if runscript == 1
    
    %Animal selection
    %-----------------------------------------------------
    %Expanimals = {'REc,''REd','REe','REf'};
    %Conanimals = {'RCa','RCb','RCc','RCd'}
    Expanimals = {'REd','REe','REf'};
    Conanimals = {'RCd'};
    
    %Filter creation
    %--------------------------------------------------------
    % epoch filter
    dayfilter = '9:10'; % Shantanu - I am adding day filter to parse out epoch filter
    epochfilter = 'isequal($type, ''run'')';
    
    % iterator
    iterator = 'epochbehaveanal';
    
    % Create filter    
    Expstimf = createfilter('animal',Expanimals,'days',dayfilter,'epochs',epochfilter,'iterator', iterator);
    Constimf = createfilter('animal',Conanimals,'days',dayfilter,'epochs',epochfilter,'iterator', iterator);
    
    % Set filter function
    % --------------------
    Expstimf = setfilterfunction(Expstimf, 'DFAsj_stim_behstats',{'DIO','pos','linpos'});
    Constimf = setfilterfunction(Constimf, 'DFAsj_stim_behstats',{'DIO','pos','linpos'});
    
    % run analysis
    % -------------
    Expstimf = runfilter(Expstimf); 
    Constimf = runfilter(Constimf);
    
    %--------------------- Finished Filter Function Run -------------------
    
    disp('Finished running filter');
    
    if savedata == 1
        clear figopt1 runscript savedata dosepdays dostimisi doisisepdays doposdistr dopossepdays dostimlearncorr doposquant dovelstim_befaft
        save(savefile);
    end
    
else
    
    load(savefile);
    
end % end runscript

if ~exist('savedata') % Jusr ran and saved data. Dont go below
    return
end

%--------------------- End Run Script -------------------

% ------------------------------
% Figure and Font Sizes

forppr = 1;
% If yes, everything set to redimscreen_figforppr1
% If not, everything set to redimscreen_figforppt1

figdir = '/data25/sjadhav/RippleInterruption/Figures/01AugSep11_RippleDisFigs/StimulationStats/SpeedAtStim/';
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
savefig1=0;
% --------Parameters-------------------------------
lockout=0.25; % in sec
nov=1:3;
fam=5:8;
day_cutoff=8;
lastdays=2;

% ---------Extract Data and Calculate ---------

str=['Exp';'Con'];
%str=['Exp'];
days = unique(Expstimf(1).epochs{1}(:,1));
allepochs = unique(Expstimf(1).epochs{1}(:,2)); % ep 2 and 4

for g = 1:size(str,1)   % Do Exp and Con groups separately

    eval(['totanim = length(',str(g,:),'stimf);']);
    for an=1:totanim % Across animals
         for d=1:length(days)
            currday = days(d);
            for ep = 1:length(allepochs)
                currep = allepochs(ep);
                index = eval(['find( (',str(g,:),'stimf(an).epochs{1}(:,1)==currday) & (',str(g,:),'stimf(an).epochs{1}(:,2)==currep) );']);
                eval([str(g,:),'stimrate(an,d,ep) =',str(g,:),'stimf(an).output{1}(index).stimrate;']);
                eval([str(g,:),'totalstim(an,d,ep) =',str(g,:),'stimf(an).output{1}(index).totalstim;']);
                eval([str(g,:),'totaltime(an,d,ep) =',str(g,:),'stimf(an).output{1}(index).totaltime;']);
                
                eval([str(g,:),'stimtime{an}{d}{ep} =',str(g,:),'stimf(an).output{1}(index).stimtime(:,1);']); % stimstarttimes
                eval([str(g,:),'stimisi{an}{d}{ep} =diff(',str(g,:),'stimtime{an}{d}{ep});']); % get stimisi
                
                eval([str(g,:),'velstim{an}{d}{ep} =',str(g,:),'stimf(an).output{1}(index).velstim;']);
                eval([str(g,:),'posstim{an}{d}{ep} =',str(g,:),'stimf(an).output{1}(index).posstim;']);
                eval([str(g,:),'vel{an}{d}{ep} =',str(g,:),'stimf(an).output{1}(index).vel;']);
                eval([str(g,:),'pos{an}{d}{ep} =',str(g,:),'stimf(an).output{1}(index).pos;']);
                
                               
                eval([str(g,:),'velstim_moving{an}{d}{ep} =',str(g,:),'stimf(an).output{1}(index).velstim_moving;']);
                eval([str(g,:),'velstim_still{an}{d}{ep} =',str(g,:),'stimf(an).output{1}(index).velstim_still;']);
                eval([str(g,:),'vel_moving{an}{d}{ep} =',str(g,:),'stimf(an).output{1}(index).vel_moving;']);
                eval([str(g,:),'vel_still{an}{d}{ep} =',str(g,:),'stimf(an).output{1}(index).vel_still;']);
                eval([str(g,:),'time_moving{an}{d}{ep} =',str(g,:),'stimf(an).output{1}(index).time_moving;']);
                eval([str(g,:),'time_still{an}{d}{ep} =',str(g,:),'stimf(an).output{1}(index).time_still;']);     
                
                % Posn at Stimln
                eval([str(g,:),'wellpos{an}{d}{ep} =',str(g,:),'stimf(an).output{1}(index).wellpos;']);
                eval([str(g,:),'intpos{an}{d}{ep} =',str(g,:),'stimf(an).output{1}(index).intpos;']);
                %eval([str(g,:),'nstimwell{an}{d}{ep} =',str(g,:),'stimf(an).output{1}(index).nstimwell;']);
                %eval([str(g,:),'nstimint{an}{d}{ep} =',str(g,:),'stimf(an).output{1}(index).nstimint;']);
                %eval([str(g,:),'nstimarm{an}{d}{ep} =',str(g,:),'stimf(an).output{1}(index).nstimarm;']);
                %eval([str(g,:),'pstimarm{an}{d}{ep} =',str(g,:),'stimf(an).output{1}(index).pstimarm;']);
                eval([str(g,:),'nstimwell(an,d,ep,:) =',str(g,:),'stimf(an).output{1}(index).nstimwell;']); % 3 values for each epoch
                eval([str(g,:),'nstimint(an,d,ep,:) =',str(g,:),'stimf(an).output{1}(index).nstimint;']); % 3 values for each epoch
                eval([str(g,:),'nstimarm(an,d,ep,:) =',str(g,:),'stimf(an).output{1}(index).nstimarm;']); % 3 values for each epoch
                eval([str(g,:),'pstimarm{an}{d}{ep} =',str(g,:),'stimf(an).output{1}(index).pstimarm;']); % 3 values for each epoch
                
                % Vel befaft
                eval([str(g,:),'velstim_befaft{an}{d}{ep} =',str(g,:),'stimf(an).output{1}(index).velstim_befaft;']);
                
                
            end % end epoch
            % Combine vel vectors for day across epochs
            d;
            eval([str(g,:),'velstim_day{an}{d} =[',str(g,:),'velstim{an}{d}{1},',str(g,:),'velstim{an}{d}{2}];']); % row vector above
            eval([str(g,:),'vel_day{an}{d} =[',str(g,:),'vel{an}{d}{1};',str(g,:),'vel{an}{d}{2}];']); % column vector above
            % Get mean vel for day
            eval([str(g,:),'meanvel_day(an,d)= nanmean(',str(g,:),'vel_day{an}{d});']);
            % Combine posn matr for day across epochs
            eval([str(g,:),'posstim_day{an}{d} =[',str(g,:),'posstim{an}{d}{1};',str(g,:),'posstim{an}{d}{2}];']); % nX2 matrix of x-y positions
            eval([str(g,:),'pos_day{an}{d} =[',str(g,:),'pos{an}{d}{1};',str(g,:),'pos{an}{d}{2}];']); % nX2 matrix of x-y positions
            % Combine stimisi vectors for day across epochs
            eval([str(g,:),'stimisi_day{an}{d} =[',str(g,:),'stimisi{an}{d}{1};',str(g,:),'stimisi{an}{d}{2}];']); % row vector above          
         
             % Take avg "posn at stimln" for the two epochs, AND GET TOTAL STIM
            eval([str(g,:),'avgnstimwell(an,d,:) = mean(squeeze([',str(g,:),'nstimwell(an,d,1,:);',str(g,:),'nstimwell(an,d,2,:)]));']);
            eval([str(g,:),'avgnstimint(an,d,:) = mean(squeeze([',str(g,:),'nstimint(an,d,1,:);',str(g,:),'nstimint(an,d,2,:)]));']);
            eval([str(g,:),'avgnstimarm(an,d,:) = mean(squeeze([',str(g,:),'nstimarm(an,d,1,:);',str(g,:),'nstimarm(an,d,2,:)]));']); % Something wrong here
            eval([str(g,:),'avgtotalstim(an,d) = mean([',str(g,:),'totalstim(an,d,1);',str(g,:),'totalstim(an,d,2)]);']);
            % Get in terms of ratio by dividing by total stim for each day
            eval([str(g,:),'avgnstimwell(an,d,:) =',str(g,:),'avgnstimwell(an,d,:)./',str(g,:),'avgtotalstim(an,d);']);
            eval([str(g,:),'avgnstimint(an,d,:) =',str(g,:),'avgnstimint(an,d,:)./',str(g,:),'avgtotalstim(an,d);']);
            
            % Combine epochs for velstim_befaft
            eval([str(g,:),'velstim_befaft_day{an}{d} =[',str(g,:),'velstim_befaft{an}{d}{1};',str(g,:),'velstim_befaft{an}{d}{2}];']); % [x,300][y,300] -> [x+y,300]
            
            
         end % end day
    end % end anim
    % Average across epochs for parameter vs. day plots
    eval([str(g,:),'stimrate_day = nanmean(',str(g,:),'stimrate,3);']);
    eval([str(g,:),'totalstim_day = nanmean(',str(g,:),'totalstim,3);']);
    eval([str(g,:),'totaltime_day = nanmean(',str(g,:),'totaltime,3);']);
 
    str(g,:)
    
    % Speed Distributions. like ripple size distribution    
    % -------------------------------------------------
    % histogram edges
    xhist =[0:1:50];
    % Initialize to store for alldays for stats
    eval([str(g,:),'vel_alldays=[];']);
    eval([str(g,:),'velstim_alldays=[];']);
    % Loop over each day and get separately
    for d=1:length(days)
        
        %currday=days(d);
        %USE DAY INDEX INSTEAD OF 9 and 10
        currday=d;
        
        % Initialize to store raw values for this day across animals for stats
        eval([str(g,:),'allvel_day{',num2str(currday),'}=[];']);
        eval([str(g,:),'allvelstim_day{',num2str(currday),'}=[];']);
        
        % Get from each animal 
         for ani=1:totanim % Loop over anim
             an=ani
             % Get current data
             eval(['currvel =',str(g,:),'vel_day{an}{currday};']);
             eval(['currvelstim =',str(g,:),'velstim_day{an}{currday};']);
             % Make histogram for current day and epoch 
             h = histc(currvel,xhist); hnorm=h./max(h);
             eval([str(g,:),'velhist_dayanim(currday,an,:) = hnorm;']);
             hh = histc(currvelstim,xhist); hhnorm=hh./max(hh);
             eval([str(g,:),'velstimhist_dayanim(currday,an,:) = hhnorm;']);
             % Store raw values for day stats 
             eval([str(g,:),'allvel_day{',num2str(currday),'}=[',str(g,:),'allvel_day{',num2str(currday),'};currvel];']);
             eval([str(g,:),'allvelstim_day{',num2str(currday),'}=[',str(g,:),'allvelstim_day{',num2str(currday),'},currvelstim];']);
             % Store raw values for alldays for grp stats
             eval([str(g,:),'vel_alldays=[',str(g,:),'vel_alldays;currvel];']);
             eval([str(g,:),'velstim_alldays=[',str(g,:),'velstim_alldays,currvelstim];']);
             % Fractions below thrs for velstim
             frac_below5 = sum(hhnorm(1:7)) ./ sum(hhnorm);
             frac_below10 = sum(hhnorm(1:12)) ./ sum(hhnorm);
             eval([str(g,:),'frac_below5(an,d)=frac_below5;']);
             eval([str(g,:),'frac_below10(an,d)=frac_below10;']);
         end % end animal
    end % end day
    
    % Speed before and after stimulation. velstim_befaft distributions. After review
    % ------------------------------------------------------------------------------
    
    binsize_befaft=100; % 100 ms bins. Speed is in 33 ms resolution. So combine 3 bins together for the velocity vector
    % Initialize to store for alldays 
    eval([str(g,:),'velstim_befaft_alldays=[];']); 
    % Using this for stats doesnt make sense for befaft - There isnt one value to compare. You want to compare PSTHs directly
    % Also, store diff in speed befor and after
    eval([str(g,:),'velstim_diff_alldays=[];']); % Difference in speed 1sec before and after stimln
    eval([str(g,:),'velstim_bef_alldays=[];']); % Before
    eval([str(g,:),'velstim_aft_alldays=[];']); % After
    eval([str(g,:),'velstim_diff2_alldays=[];']); % Difference in speed 2 sec before and after stimln
    
    for d=1:length(days)
        
        %currday=days(d);
        %USE DAY INDEX INSTEAD OF 9 and 10
        currday=d;
        
        % Initialize to store raw values for this day across animals for stats
        eval([str(g,:),'allvelstim_befaft_day{',num2str(currday),'}=[];']); % Again - not necessary for stats
        
         % Get from each animal 
         for ani=1:totanim % Loop over anim
             an=ani;
             
             % Save for each animal separately - initialize for each animal of day=1
             if d==1
                 eval([str(g,:),'velstim_befaft_byanim{',num2str(an),'}=[];']);
                 eval([str(g,:),'velstim_diff_byanim{',num2str(an),'}=[];']);
             end
             
             
             % Get current data
             eval(['currvelstim_befaft =',str(g,:),'velstim_befaft_day{an}{currday};']);
             currvelstim_befaft=reshape(currvelstim_befaft,size(currvelstim_befaft,1),3,100); % 300pts, each pt 33.33ms,
             %currvelstim_befaft=reshape(currvelstim_befaft,size(currvelstim_befaft,1),15,20);
             % For 100ms bins, combine 3 bins. For 500ms, combine 15 bins 
             currvelstim_befaft = squeeze(mean(currvelstim_befaft,2)); 
             % This is the histogram in 100ms bins/0.5sec bins, 5 sec (aka 50 bins) before and after stimulation
             
             % Get speed in 1 sec before and after stim. If 100ms bins, this is mean of 10 bins before and after
             currvelstim_bef = mean(currvelstim_befaft(:,41:50),2); 
             currvelstim_aft = mean(currvelstim_befaft(:,51:60),2); 
             currvelstim_diff = currvelstim_bef-currvelstim_aft; 
             
             currvelstim_bef2 = mean(currvelstim_befaft(:,31:50),2); 
             currvelstim_aft2 = mean(currvelstim_befaft(:,51:70),2); 
             currvelstim_diff2 = currvelstim_bef2-currvelstim_aft2; 
             
             % Save for each day separately
             eval([str(g,:),'allvelstim_befaft_day{',num2str(currday),'}=[',str(g,:),'allvelstim_befaft_day{',num2str(currday),'};currvelstim_befaft];']);
             % And Save for alldays combined
             eval([str(g,:),'velstim_befaft_alldays=[',str(g,:),'velstim_befaft_alldays;currvelstim_befaft];']);
             % And, Save Diff in speed due to stimln for all days
             eval([str(g,:),'velstim_bef_alldays=[',str(g,:),'velstim_bef_alldays;currvelstim_bef];']);
             eval([str(g,:),'velstim_aft_alldays=[',str(g,:),'velstim_aft_alldays;currvelstim_aft];']);
             eval([str(g,:),'velstim_diff_alldays=[',str(g,:),'velstim_diff_alldays;currvelstim_diff];']);
             eval([str(g,:),'velstim_diff2_alldays=[',str(g,:),'velstim_diff2_alldays;currvelstim_diff2];']);
             % Save for each animal separately
             eval([str(g,:),'velstim_befaft_byanim{',num2str(an),'}=[',str(g,:),'velstim_befaft_byanim{',num2str(an),'};currvelstim_befaft];']); % matrix
             eval([str(g,:),'velstim_diff_byanim{',num2str(an),'}=[',str(g,:),'velstim_diff_byanim{',num2str(an),'};currvelstim_diff];']); % vector
             
         end
    end
        
    
    
    % StimISI distributions
    % ---------------------
    % histogram edges
    xhist =[0.25:0.1:5];
    % Initialize to store for alldays for stats
    eval([str(g,:),'stimisi_alldays=[];']);
    % Loop over each day and get separately
    for d=1:length(days)
        
        %currday=days(d);
        %USE DAY INDEX INSTEAD OF 9 and 10
        currday=d;
        
        % Initialize to store raw values for this day across animals for stats
        eval([str(g,:),'allstimisi_day{',num2str(currday),'}=[];']);
        
        % Get from each animal 
         for ani=1:totanim % Loop over anim
             an=ani;
             % Get current data
             eval(['currstimisi =',str(g,:),'stimisi_day{an}{currday};']);
             % Make histogram for current day and epoch 
             hh = histc(currstimisi,xhist); hhnorm=hh./max(hh);
             eval([str(g,:),'stimisihist_dayanim(currday,an,:) = hhnorm;']);
             % Store raw values for day stats 
             eval([str(g,:),'allstimisi_day{',num2str(currday),'}=[',str(g,:),'allstimisi_day{',num2str(currday),'};currstimisi];']);
             % Store raw values for alldays for grp stats
             eval([str(g,:),'stimisi_alldays=[',str(g,:),'stimisi_alldays;currstimisi];']);      
         end % end animal
    end % end day
    
    
    
    % Position matrices. Can I Combine across animals for REd,REe,REf and RCb,RCc,RCd - same tracks 
    % --------------------
    % Initialize to store for alldays 
    eval([str(g,:),'pos_alldays=[];']);
    eval([str(g,:),'posstim_alldays=[];']);
    % Loop over each day and get separately
    for d=1:length(days)
                
        %currday=days(d);
        %USE DAY INDEX INSTEAD OF 9 and 10
        currday=d;
        
        % Initialize to store raw values for this day across animals for stats
        eval([str(g,:),'allpos_day{',num2str(currday),'}=[];']);
        eval([str(g,:),'allposstim_day{',num2str(currday),'}=[];']);
        
        % Get from each animal 
        % for ani=2:totanim % Loop over anim
        % an=ani-1;
         for ani=1:1 % Loop over anim
             an=ani;
             % Get current data
             eval(['currpos =',str(g,:),'pos_day{an}{currday};']);
             eval(['currposstim =',str(g,:),'posstim_day{an}{currday};']);
             
             % Store raw values for day across animals
             eval([str(g,:),'allpos_day{',num2str(currday),'}=[',str(g,:),'allpos_day{',num2str(currday),'};currpos];']); % nX2 matrix of x-y positions
             eval([str(g,:),'allposstim_day{',num2str(currday),'}=[',str(g,:),'allposstim_day{',num2str(currday),'};currposstim];']); % nX2 matrix of x-y positions
             % Store raw values for alldays for grp stats
             eval([str(g,:),'pos_alldays=[',str(g,:),'pos_alldays;currpos];']);
             eval([str(g,:),'posstim_alldays=[',str(g,:),'posstim_alldays;currposstim];']);
          end % end animal
    end % end day
       
end % end str=Exp or Con
    

%****************************************
% Figures -
%****************************************

if figopt1 == 1
    
    % First make plots similar to sj_stimstats_summ
    % ---------------------------------------------
    
    % ************* Stim Rate ******************
    % 1) Avg Stim Rate for each Animal
    figure; hold on;
    if forppr==1
        redimscreen_figforppr1;
    else
        redimscreen_figforppt1;
    end
    
    Exp_animstimrate = mean(Expstimrate_day,2);
    Con_animstimrate = mean(Constimrate_day,2);
    plot(1*ones(size(Exp_animstimrate)), Exp_animstimrate,'ro','MarkerSize',12,'Linewidth',3);
    plot(2*ones(size(Con_animstimrate)), Con_animstimrate,'bd','MarkerSize',12,'Linewidth',3);
    [p_expcon_animstimrate,h_expcon_animstimrate] = ranksum(Con_animstimrate,Exp_animstimrate); %p=0.6857
    
    axis([0.5 2.5 0 2]);
    set(gca,'XTick',[1 2],'XTickLabel',{'Exp';'Con'}','FontSize',xfont,'Fontweight','normal');
    title('Mean Stimlation Rate','FontSize',tfont,'Fontweight','normal');
    ylabel('Mean Stimlation Rate','FontSize',yfont,'Fontweight','normal');
    set(gca,'YLim',[0 2.2]);
    if savefig1==1,
        figfile = [figdir,'ConvsExp_AnimStimRate'];
        print('-dpdf', figfile);
        print('-djpeg', figfile);
        saveas(gcf,figfile,'fig');
    end
    
    % 2) Avg Stim Rate for Group
    figure; hold on;
    if forppr==1
        redimscreen_figforppr1;
    else
        redimscreen_figforppt1;
    end
    
    Exp_meanstimrate = mean(mean(Expstimrate_day));
    Exp_errstimrate = sem(sem(Expstimrate_day));
    Con_meanstimrate = mean(mean(Constimrate_day));
    Con_errstimrate = sem(sem(Constimrate_day));
  
    errorbar([1],[Exp_meanstimrate],[Exp_errstimrate],'ro','MarkerSize',12,'Linewidth',3);
    errorbar([2],[Con_meanstimrate],[Con_errstimrate],'bd','MarkerSize',12,'Linewidth',3);
    [h_expcon_meanstimrate,p_expcon_meanstimrate] = ttest2(Expstimrate_day(:),Constimrate_day(:)); %p=0.3340
    
    axis([0.5 2.5 0 1.25]);
    set(gca,'XTick',[1 2],'XTickLabel',{'Exp';'Con'}','FontSize',xfont,'Fontweight','normal');
    title('Mean Stimlation Rate','FontSize',tfont,'Fontweight','normal');
    ylabel('Mean Stimlation Rate','FontSize',yfont,'Fontweight','normal');
    if savefig1==1,
        figfile = [figdir,'ConvsExp_MeanStimRate'];
        print('-dpdf', figfile);
        print('-djpeg', figfile);
        saveas(gcf,figfile,'fig');
    end
    
    % 3) Avg Stim Rate vs day
    figure; hold on;
    if forppr==1
        redimscreen_figforppr1;
    else
        redimscreen_figforppt1;
    end
    
    Exp_daystimrate = mean(Expstimrate_day,1); Exp_daystimrateerr = sem(Expstimrate_day,1);
    Con_daystimrate = mean(Constimrate_day,1); Con_daystimrateerr = sem(Constimrate_day,1);
  
    errorbar(days,Exp_daystimrate,Exp_daystimrateerr,'ro-','MarkerSize',12,'Linewidth',3);
    errorbar(days+0.1,Con_daystimrate,Con_daystimrateerr,'bd-','MarkerSize',12,'Linewidth',3);
    [panovarm_stimrate] = anova_rm({Expstimrate_day Constimrate_day},'off');
    %pgrp=0.7247, ptime=0.0064, pint=0.1354
    
    % Make Plot presentable
    title('Stimulation Rate vs Days','FontSize',tfont,'Fontweight','normal');
    xlabel('Day','FontSize',xfont,'Fontweight','normal');
    ylabel('Stimulation Rate vs Days','FontSize',yfont,'Fontweight','normal');
    %legend('Location','NorthEast');
    axis([0.5 max(days)+0.5 0 2]); 
    if savefig1==1,
        figfile = [figdir,'ConvsExp_MeanstimrateVsDays'];
        print('-dpdf', figfile);
        print('-djpeg', figfile);
        saveas(gcf,figfile,'fig');
    end
    
    
    % ************* Mean speed ******************
    % 1) Mean Speed for each animal
    figure; hold on;
    if forppr==1
        redimscreen_figforppr1;
    else
        redimscreen_figforppt1;
    end
    Exp_animmeanvel = mean(Expmeanvel_day,2);
    Con_animmeanvel = mean(Conmeanvel_day,2);
    plot(1*ones(size(Exp_animmeanvel)), Exp_animmeanvel,'ro','MarkerSize',12,'Linewidth',3);
    plot(2*ones(size(Con_animmeanvel)), Con_animmeanvel,'bd','MarkerSize',12,'Linewidth',3);
    [p_expcon_animmeanvel,h_expcon_animmeanvel] = ranksum(Con_animmeanvel,Exp_animmeanvel); %p=0.4857
     
    axis([0.5 2.5 0 15]);
    set(gca,'XTick',[1 2],'XTickLabel',{'Exp';'Con'}','FontSize',xfont,'Fontweight','normal');
    title('Mean Speed','FontSize',tfont,'Fontweight','normal');
    ylabel('Mean Speed','FontSize',yfont,'Fontweight','normal');
    if savefig1==1,
        figfile = [figdir,'ConvsExp_MeanSpeed'];
        print('-dpdf', figfile);
        print('-djpeg', figfile);
        saveas(gcf,figfile,'fig');
    end
    
    % 2) Avg Speed vs day
    figure; hold on;
    if forppr==1
        redimscreen_figforppr1;
    else
        redimscreen_figforppt1;
    end
    
    Exp_daymeanvel = mean(Expmeanvel_day,1); Exp_daymeanvelerr = sem(Expmeanvel_day,1);
    Con_daymeanvel = mean(Conmeanvel_day,1); Con_daymeanvelerr = sem(Conmeanvel_day,1);
    errorbar(days,Exp_daymeanvel,Exp_daymeanvelerr,'ro-','MarkerSize',12,'Linewidth',3);
    errorbar(days+0.1,Con_daymeanvel,Con_daymeanvelerr,'bd-','MarkerSize',12,'Linewidth',3);
    [panovarm_meanvel] = anova_rm({Expmeanvel_day Conmeanvel_day},'off');
    %pgrp=0.3289, ptime=0, pint=0.3288
    
    % Make Plot presentable
    title('Mean Speed vs Days','FontSize',tfont,'Fontweight','normal');
    xlabel('Day','FontSize',xfont,'Fontweight','normal');
    ylabel('Mean Speed','FontSize',yfont,'Fontweight','normal');
    %legend('Location','NorthEast');
    axis([0.5 max(days)+0.5 0 15]); 
    if savefig1==1,
        figfile = [figdir,'ConvsExp_MeanSpeedVsDays'];
        print('-dpdf', figfile);
        print('-djpeg', figfile);
        saveas(gcf,figfile,'fig');
    end
    
    % ************* Speed at stim - Fraction below thrs ******************
    % 1) Mean for each animal
    figure; hold on;
    redimscreen_2horsubplots; jit=[0,0.05,0.1,0.15];
    subplot(1,2,1); hold on;
    Exp_animfracbelow5 = mean(Expfrac_below5,2);
    Con_animfracbelow5 = mean(Confrac_below5,2);
    plot(1+jit, Exp_animfracbelow5,'ro','MarkerSize',12,'Linewidth',3);
    plot(2+jit, Con_animfracbelow5,'bd','MarkerSize',12,'Linewidth',3);
    [p_expcon_animfracbelow5,h_expcon_animfracbelow5] = ranksum(Con_animfracbelow5,Exp_animfracbelow5); %p=1
    axis([0.5 2.5 0 1]);
    set(gca,'XTick',[1 2],'XTickLabel',{'Exp';'Con'}','FontSize',xfont,'Fontweight','normal');
    title('Prop of stim below 5 cm/s','FontSize',tfont,'Fontweight','normal');
    ylabel('Proportion of stimulations','FontSize',yfont,'Fontweight','normal');
    subplot(1,2,2); hold on;
    Exp_animfracbelow10 = mean(Expfrac_below10,2);
    Con_animfracbelow10 = mean(Confrac_below10,2);
    plot(1+jit, Exp_animfracbelow10,'ro','MarkerSize',12,'Linewidth',3);
    plot(2+jit, Con_animfracbelow10,'bd','MarkerSize',12,'Linewidth',3);
    [p_expcon_animfracbelow10,h_expcon_animfracbelow10] = ranksum(Con_animfracbelow10,Exp_animfracbelow10); %p=0.3429
    axis([0.5 2.5 0 1]);
    set(gca,'XTick',[1 2],'XTickLabel',{'Exp';'Con'}','FontSize',xfont,'Fontweight','normal');
    title('Prop of stim below 10 cm/s','FontSize',tfont,'Fontweight','normal');
    ylabel('Proportion of stimulations','FontSize',yfont,'Fontweight','normal');
    
    if savefig1==1,
        figfile = [figdir,'ConvsExp_Fracbelowthrs'];
        print('-dpdf', figfile);
        print('-djpeg', figfile);
        saveas(gcf,figfile,'fig');
    end
    
    % 2) Avg frac vs day
    figure; hold on;
    redimscreen_2horsubplots;
    subplot(1,2,1); hold on;
    Exp_dayfracbelow5 = mean(Expfrac_below5,1); Exp_dayfracbelow5err = sem(Expfrac_below5,1);
    Con_dayfracbelow5 = mean(Confrac_below5,1); Con_dayfracbelow5err = sem(Confrac_below5,1);
    errorbar(days,Exp_dayfracbelow5,Exp_dayfracbelow5err,'ro-','MarkerSize',12,'Linewidth',3);
    errorbar(days+0.1,Con_dayfracbelow5,Con_dayfracbelow5err,'bd-','MarkerSize',12,'Linewidth',3);
    [panovarm_fracbelow5] = anova_rm({Expfrac_below5 Confrac_below5},'off');
    %pgrp=0.8253, ptime=0.9434, pint=0.1211
    title('Prop of stim below 5 cm/s vs Days','FontSize',tfont,'Fontweight','normal');
    xlabel('Day','FontSize',xfont,'Fontweight','normal');
    ylabel('Prop of stim below 5 cm/s','FontSize',yfont,'Fontweight','normal');
    %legend('Location','NorthEast');
    axis([0.5 max(days)+0.5 0 1]); 
    
    subplot(1,2,2); hold on;
    Exp_dayfracbelow10 = mean(Expfrac_below10,1); Exp_dayfracbelow10err = sem(Expfrac_below10,1);
    Con_dayfracbelow10 = mean(Confrac_below10,1); Con_dayfracbelow10err = sem(Confrac_below10,1);
    errorbar(days,Exp_dayfracbelow10,Exp_dayfracbelow10err,'ro-','MarkerSize',12,'Linewidth',3);
    errorbar(days+0.1,Con_dayfracbelow10,Con_dayfracbelow10err,'bd-','MarkerSize',12,'Linewidth',3);
    [panovarm_fracbelow10] = anova_rm({Expfrac_below10 Confrac_below10},'off');
    %pgrp=0.5688, ptime=0.1143, pint=0.1606
    title('Prop of stim below 10 cm/s vs Days','FontSize',tfont,'Fontweight','normal');
    xlabel('Day','FontSize',xfont,'Fontweight','normal');
    ylabel('Prop of stim below 10 cm/s','FontSize',yfont,'Fontweight','normal');
    %legend('Location','NorthEast');
    axis([0.5 max(days)+0.5 0 1]); 
    
    if savefig1==1,
        figfile = [figdir,'ConvsExp_FracbelowthrsVsDays'];
        print('-dpdf', figfile);
        print('-djpeg', figfile);
        saveas(gcf,figfile,'fig');
    end
    
    
end % end figopt1


% *************************************
% Distributions - Combined for all days
% *************************************
  
xhist =[0:1:50];
if dodistr==1
   
     Conmeanvelhist = mean(squeeze(mean(Convelhist_dayanim,2))); % Mean along animals, then days
     Conmeanvelstimhist = mean(squeeze(mean(Convelstimhist_dayanim,2))); % Mean along animals, then days
     Expmeanvelhist = mean(squeeze(mean(Expvelhist_dayanim,2))); % Mean along animals, then days
     Expmeanvelstimhist = mean(squeeze(mean(Expvelstimhist_dayanim,2))); % Mean along animals, then days
     Conerrvelhist = sem(squeeze(sem(Convelhist_dayanim,2))); 
     Conerrvelstimhist = sem(squeeze(sem(Convelstimhist_dayanim,2)));
     Experrvelhist = sem(squeeze(sem(Expvelhist_dayanim,2))); 
     Experrvelstimhist = sem(squeeze(sem(Expvelstimhist_dayanim,2))); 
     
     % Vectors of histogram - for stats
     Convelhistvec_alldays = Convelhist_dayanim(:);
     Expvelhistvec_alldays = Expvelhist_dayanim(:);
     Convelstimhistvec_alldays = Convelstimhist_dayanim(:);
     Expvelstimhistvec_alldays = Expvelstimhist_dayanim(:);
     
     % All speeds during behavior
     figure; hold on;
     if forppr==1
         redimscreen_figforppr1;
     else
         redimscreen_figforppt1;
     end
     plot(xhist,Conmeanvelhist,'b.-','Linewidth',3,'MarkerSize',18);
     plot(xhist,Expmeanvelhist,'r.-','Linewidth',3,'MarkerSize',18);
     jbfill(xhist, [Conmeanvelhist+Conerrvelhist],[Conmeanvelhist-Conerrvelhist],'b','b',1,0.3);
     jbfill(xhist, [Expmeanvelhist+Experrvelhist],[Expmeanvelhist-Experrvelhist],'r','r',1,0.3);
     set(gca,'YLim',[0 1.02]); set(gca,'XLim',[-1 max(xhist)+1]);
     title(['Con and Exp - Speed distr during beh. All Days'],'FontSize',16,'Fontweight','normal');
     ylabel('Proportion','FontSize',yfont,'Fontweight','normal');
     xlabel('Speed(cm/s)','FontSize',xfont,'Fontweight','normal');
     % ----------- Stats -------------
     [h_vel,p_vel] = ttest2(Convel_alldays,Expvel_alldays);
     [h_velhist,p_velhist] = kstest2(Convelhistvec_alldays,Expvelhistvec_alldays);
     
     % Speeds at stimulation
     figure; hold on;
     if forppr==1
         redimscreen_figforppr1;
     else
         redimscreen_figforppt1;
     end
     plot(xhist,Conmeanvelstimhist,'b.-','Linewidth',3,'MarkerSize',18);
     plot(xhist,Expmeanvelstimhist,'r.-','Linewidth',3,'MarkerSize',18);
     jbfill(xhist, [Conmeanvelstimhist+Conerrvelstimhist],[Conmeanvelstimhist-Conerrvelstimhist],'b','b',1,0.3);
     jbfill(xhist, [Expmeanvelstimhist+Experrvelstimhist],[Expmeanvelstimhist-Experrvelstimhist],'r','r',1,0.3);
     set(gca,'YLim',[0 1.02]); set(gca,'XLim',[0 max(xhist)+1]);
     title(['Con and Exp - Speed at stim. All Days'],'FontSize',16,'Fontweight','normal');
     ylabel('Proportion','FontSize',yfont,'Fontweight','normal');
     xlabel('Speed(cm/s)','FontSize',xfont,'Fontweight','normal');
     % ----------- Stats -------------
     [h_velstim,p_velstim] = ttest2(Convelstim_alldays,Expvelstim_alldays);
     [h_velstimhist,p_velstimhist] = kstest2(Convelstimhistvec_alldays,Expvelstimhistvec_alldays);
     
end % end dodistr


% ***************************************
% Distributions - Separately for all days
% ***************************************

if dosepdays==1
    
    for n = 1:length(days)
            currday = days(n);
            Conmeanvelhist = squeeze(mean(Convelhist_dayanim(currday,:,:),2))'; % Mean along 2nd dimension of animals
            Conerrvelhist = squeeze(sem(Convelhist_dayanim(currday,:,:),2))'; % sem along 2nd dimension of animals
            Conmeanvelstimhist = squeeze(mean(Convelstimhist_dayanim(currday,:,:),2))'; 
            Conerrvelstimhist = squeeze(sem(Convelstimhist_dayanim(currday,:,:),2))';
            Expmeanvelhist = squeeze(mean(Expvelhist_dayanim(currday,:,:),2))'; % Mean along 2nd dimension of animals
            Experrvelhist = squeeze(sem(Expvelhist_dayanim(currday,:,:),2))'; % sem along 2nd dimension of animals
            Expmeanvelstimhist = squeeze(mean(Expvelstimhist_dayanim(currday,:,:),2))'; 
            Experrvelstimhist = squeeze(sem(Expvelstimhist_dayanim(currday,:,:),2))';
            
            % Vectors of histogram for current day - for stats
            Convelhistvec_currday = squeeze(Convelhist_dayanim(currday,:,:)); Convelhistvec(currday,:) = Convelhistvec_currday(:);
            Convelstimhistvec_currday = squeeze(Convelstimhist_dayanim(currday,:,:)); Convelstimhistvec(currday,:) = Convelstimhistvec_currday(:);
            Expvelhistvec_currday = squeeze(Expvelhist_dayanim(currday,:,:)); Expvelhistvec(currday,:) = Expvelhistvec_currday(:);
            Expvelstimhistvec_currday = squeeze(Expvelstimhist_dayanim(currday,:,:)); Expvelstimhistvec(currday,:) = Expvelstimhistvec_currday(:);
            
            % Plot for current day
            figure; hold on;
            redimscreen_2horsubplots;
            subplot(1,2,1); hold on;
            plot(xhist,Conmeanvelhist,'b.-','Linewidth',3,'MarkerSize',18);
            plot(xhist,Expmeanvelhist,'r.-','Linewidth',3,'MarkerSize',18);
            jbfill(xhist, [Conmeanvelhist+Conerrvelhist],[Conmeanvelhist-Conerrvelhist],'b','b',1,0.3);
            jbfill(xhist, [Expmeanvelhist+Experrvelhist],[Expmeanvelhist-Experrvelhist],'r','r',1,0.3);
            set(gca,'YLim',[0 1.02]); set(gca,'XLim',[-1 max(xhist)+1]);
            title(['Con and Exp - Speed distr during beh. Day ',num2str(currday)],'FontSize',16,'Fontweight','normal');
            ylabel('Proportion','FontSize',yfont,'Fontweight','normal');
            xlabel('Speed(cm/s)','FontSize',xfont,'Fontweight','normal');
            % ----------- Stats -------------
            [h_velday(n),p_velday(n)] = ttest2(Conallvel_day{currday},Expallvel_day{currday});
            [h_velhistday(n),p_velhistday(n)] = kstest2(Convelhistvec(currday,:),Expvelhistvec(currday,:));
            
            subplot(1,2,2); hold on;
            plot(xhist,Conmeanvelstimhist,'b.-','Linewidth',3,'MarkerSize',18);
            plot(xhist,Expmeanvelstimhist,'r.-','Linewidth',3,'MarkerSize',18);
            jbfill(xhist, [Conmeanvelstimhist+Conerrvelstimhist],[Conmeanvelstimhist-Conerrvelstimhist],'b','b',1,0.3);
            jbfill(xhist, [Expmeanvelstimhist+Experrvelstimhist],[Expmeanvelstimhist-Experrvelstimhist],'r','r',1,0.3);
            set(gca,'YLim',[0 1.02]); set(gca,'XLim',[-1 max(xhist)+1]);
            title(['Con and Exp - Speed at stim. Day ',num2str(currday)],'FontSize',16,'Fontweight','normal');
            ylabel('Proportion','FontSize',yfont,'Fontweight','normal');
            xlabel('Speed(cm/s)','FontSize',xfont,'Fontweight','normal');
            % ----------- Stats -------------
            [h_velstimday(n),p_velstimday(n)] = ttest2(Conallvelstim_day{currday},Expallvelstim_day{currday});
            [h_velstimhistday(n),p_velstimhistday(n)] = kstest2(Convelstimhistvec(currday,:),Expvelstimhistvec(currday,:));        
    end % end days
end % end dosepdays


% *************************************
% Stim ISI Distributions - All Days
% *************************************

xhist=[0.25:0.1:5];
if dostimisi==1
   
     Conmeanstimisihist = mean(squeeze(mean(Constimisihist_dayanim,2))); % Mean along animals, then days
     Expmeanstimisihist = mean(squeeze(mean(Expstimisihist_dayanim,2))); % Mean along animals, then days
     Conerrstimisihist = sem(squeeze(sem(Constimisihist_dayanim,2)));
     Experrstimisihist = sem(squeeze(sem(Expstimisihist_dayanim,2))); 
     
     % Vectors of histogram - for stats
     Constimisihistvec_alldays = Constimisihist_dayanim(:);
     Expstimisihistvec_alldays = Expstimisihist_dayanim(:);

     % StimISI
     figure; hold on;
     if forppr==1
         redimscreen_figforppr1;
     else
         redimscreen_figforppt1;
     end
     plot(xhist,Conmeanstimisihist,'b.-','Linewidth',1,'MarkerSize',18);
     plot(xhist,Expmeanstimisihist,'r.-','Linewidth',1,'MarkerSize',18);
     jbfill(xhist, [Conmeanstimisihist+Conerrstimisihist],[Conmeanstimisihist-Conerrstimisihist],'b','b',1,0.3);
     jbfill(xhist, [Expmeanstimisihist+Experrstimisihist],[Expmeanstimisihist-Experrstimisihist],'r','r',1,0.3);
     set(gca,'YLim',[0 1.02]); set(gca,'XLim',[0 max(xhist)+1]);
     set(gca,'XLim',[0 3.05]);
     title(['Con and Exp - Stim ISI. All Days'],'FontSize',16,'Fontweight','normal');
     ylabel('Proportion','FontSize',yfont,'Fontweight','normal');
     xlabel('Stim ISI (sec)','FontSize',xfont,'Fontweight','normal');
     % ----------- Stats -------------
     [h_stimisi,p_stimisi] = ttest2(Constimisi_alldays,Expstimisi_alldays);
     [h_stimisihist,p_stimisihist] = kstest2(Constimisihistvec_alldays,Expstimisihistvec_alldays);
     
end % end dostimisi

% ***************************************
% Stim ISI Distributions - Separately for days
% ***************************************

if doisisepdays==1
    
    for n = 1:length(days)
            currday = days(n);
            Conmeanstimisihist = squeeze(mean(Constimisihist_dayanim(currday,:,:),2))'; 
            Conerrstimisihist = squeeze(sem(Constimisihist_dayanim(currday,:,:),2))';
            Expmeanstimisihist = squeeze(mean(Expstimisihist_dayanim(currday,:,:),2))'; 
            Experrstimisihist = squeeze(sem(Expstimisihist_dayanim(currday,:,:),2))';
            
            % Vectors of histogram for current day - for stats
            Constimisihistvec_currday = squeeze(Constimisihist_dayanim(currday,:,:)); Constimisihistvec(currday,:) = Constimisihistvec_currday(:);
            Expstimisihistvec_currday = squeeze(Expstimisihist_dayanim(currday,:,:)); Expstimisihistvec(currday,:) = Expstimisihistvec_currday(:);
            
            % Plot for current day
            figure; hold on;
            if forppr==1
                redimscreen_figforppr1;
            else
                redimscreen_figforppt1;
            end
            plot(xhist,Conmeanstimisihist,'b.-','Linewidth',1,'MarkerSize',18);
            plot(xhist,Expmeanstimisihist,'r.-','Linewidth',1,'MarkerSize',18);
            jbfill(xhist, [Conmeanstimisihist+Conerrstimisihist],[Conmeanstimisihist-Conerrstimisihist],'b','b',1,0.3);
            jbfill(xhist, [Expmeanstimisihist+Experrstimisihist],[Expmeanstimisihist-Experrstimisihist],'r','r',1,0.3);
            set(gca,'YLim',[0 1.02]); set(gca,'XLim',[-1 max(xhist)+1]);
            title(['Con and Exp - Speed at stim. Day ',num2str(currday)],'FontSize',16,'Fontweight','normal');
            ylabel('Proportion','FontSize',yfont,'Fontweight','normal');
            xlabel('Speed(cm/s)','FontSize',xfont,'Fontweight','normal');
            % ----------- Stats -------------
            [h_stimisiday(n),p_stimisiday(n)] = ttest2(Conallstimisi_day{currday},Expallstimisi_day{currday});
            [h_stimisihistday(n),p_stimisihistday(n)] = kstest2(Constimisihistvec(currday,:),Expstimisihistvec(currday,:));        
    end % end days
end % end dosepdays












% *************************************
% Position Distributions - Combined for all days
% *************************************
   
if binsize==1, maxlim = 115; end
if binsize==2, maxlim = 57; end
if doposdistr==1
    
    str=['Exp';'Con'];
   for g = 1:size(str,1)       
        % Do density of stimulations if multiple animals
        % ------------------------------------------------
        % Gaussian smooth stim density
        % Call posn at stimln spikes - use same code as openfieldrate
        eval(['tmpposition = ',str(g,:),'pos_alldays;']);
        eval(['tmpspikes = ',str(g,:),'posstim_alldays;']);
        minx = floor(min(tmpposition(:,1))) - 1;
        maxx = ceil(max(tmpposition(:,1))) + 1;
        binx = (minx:binsize:maxx);
        miny = floor(min(tmpposition(:,2))) - 1;
        maxy = ceil(max(tmpposition(:,2))) + 1;
        biny = (miny:binsize:maxy);
        
        [spikes BX BY] = hist2(tmpspikes(:,1), tmpspikes(:,2), binx, biny);
        eval([str(g,:),'posnstim = spikes;']); 
        % Get occupancy also
        [occupancy xticks yticks] = hist2(tmpposition(:,1), tmpposition(:,2), binx, biny);
        eval([str(g,:),'posocc = occupancy;']); 
        
        % Forget rate - just get nstim
        %         % Get occupancy in case you want stim rate
        %         [occupancy xticks yticks] = hist2(tmpposition(:,1), tmpposition(:,2), binx, biny);
        %         tmpspikes = Expposstim_alldays;
        %         %1) Get stimrate from stim and occupancy
        %         nonzero = find(occupancy ~= 0);
        %         spikerate = zeros(size(spikes));
        %         spikerate(nonzero) = spikes(nonzero) ./(timestep* occupancy(nonzero) );
        %         %2) Smooth stimrate and occupancy
        %         g = gaussian2(std,(6*std));
        %         smoothedstimrate = filter2(g,(spikerate)); % is this the right filter?
        %         %smoothedoccupancy = [];
        %         %smoothedoccupancy = zeros(size(output.spikes));
        %         smoothedoccupancy = filter2(g, occupancy);
        %         % 3) Turn occupancy to seconds and set spikerate wherever occupancy
        %         %is < threshold occupancy in seconds to 0
        %         occupancy = timestep*occupancy;
        %         smoothedoccupancy = timestep*smoothedoccupancy;
        %         %zero = find(smoothedoccupancy == 0);
        %         zero = find(smoothedoccupancy <= threshocc);
        %         %smoothedspikerate(zero) = -2; %no occupancy is negative and diff/darker from occupancy but no spikes
        %         smoothedspikerate(zero) = -1;
        %         eval([str(g,:),'posstimrate = smoothedspikerate;']);
    end
    
    % Plot map of nstim
    figure; hold on;
    redimscreen_2horsubplots;
    
    subplot(1,2,1); hold on;
    
    Expposnstim_norm = Expposnstim./max(max(Expposnstim));
    % Set all occupancy<thrsocc to just below 0. bounds will be from this no to 1. timestep=0.0334. 30pts~1 sec. 6pts=0.2sec
    Expposnstim_norm(find( Expposocc < 6))=-0.01; 
    cmap = jet(1024); cmap(1,:) = 1; colormap(cmap);
    bounds = [-0.01,1];
    h = imagesc(Expposnstim_norm,bounds); 
    ch = colorbar;
    title(['Exp NStim-Pos. All Days'],'FontSize',16,'Fontweight','normal');
    ylabel('y-pos','FontSize',yfont,'Fontweight','normal');
    xlabel('x-pos','FontSize',xfont,'Fontweight','normal');
    
    subplot(1,2,2); hold on;
    Conposnstim_norm = Conposnstim./max(max(Conposnstim));
    % Set all occupancy<thrsocc to just below 0. bounds will be from this no to 1
    Conposnstim_norm(find( Conposocc < 6))=-0.01; 
    cmap = jet(1024); cmap(1,:) = 1; colormap(cmap);
    bounds = [-0.01,1];
    h = imagesc(Conposnstim_norm,bounds);
    ch = colorbar;
    title(['Con NStim-Pos. All Days'],'FontSize',16,'Fontweight','normal');
    ylabel('y-pos','FontSize',yfont,'Fontweight','normal');
    xlabel('x-pos','FontSize',xfont,'Fontweight','normal');
           
    if savefig1==1,
        figfile = [figdir,'REdAndRCb_PosAtStimMap_AllDays'];
        print('-dpdf', figfile); print('-djpeg', figfile); saveas(gcf,figfile,'fig');
    end
    
end % doposdistr



if dopossepdays==1
    str=['Exp';'Con'];
    % Make 1 figure for control and 1 for Exp
    % These are for dots
%     figure(50); hold on; redimscreen;
%     figure(100); hold on; redimscreen;
    
    % These are for maps
    figure(51); hold on; redimscreen;
    figure(101); hold on; redimscreen;
    
    for n = 1:length(days)
        currday = days(n);       
        for g = 1:size(str,1)
            % Call posn at stimln spikes - use same code as openfieldrate
            eval(['tmpposition = ',str(g,:),'allpos_day{currday};']);
            eval(['tmpspikes = ',str(g,:),'allposstim_day{currday};']);
            
            % Can do density of stimulations or stim rate or norm stim
            % Do density of stimulations if multiple animals
            % Gaussian smooth stim density
            minx = floor(min(tmpposition(:,1))) - 1;
            maxx = ceil(max(tmpposition(:,1))) + 1;
            binx = (minx:binsize:maxx);
            miny = floor(min(tmpposition(:,2))) - 1;
            maxy = ceil(max(tmpposition(:,2))) + 1;
            biny = (miny:binsize:maxy);
            [spikes BX BY] = hist2(tmpspikes(:,1), tmpspikes(:,2), binx, biny);
            eval([str(g,:),'posnstim = spikes;']);
            
            eval([str(g,:),'allposcurr = tmpposition;']);
            eval([str(g,:),'posstimcurr = tmpspikes;']);
            
            % Get occupancy also
            [occupancy xticks yticks] = hist2(tmpposition(:,1), tmpposition(:,2), binx, biny);
            eval([str(g,:),'posocc = occupancy;']); 

        end
        
        % Plot map of nstim
        
        figure(51); hold on;
        subplot(2,4,n); hold on;
        
        Conposnstim_norm = Conposnstim./max(max(Conposnstim));
        % Set all occupancy<thrsocc to just below 0. bounds will be from this no to 1. timestep=0.0334. 30pts~1 sec. 6pts=0.2sec
        Conposnstim_norm(find( Conposocc < 3))=-0.01;
        cmap = jet(1024); cmap(1,:) = 1; colormap(cmap);
        bounds = [-0.01,1];
        h = imagesc(Conposnstim_norm,bounds);
        %ch = colorbar;
        set(gca,'YLim',[0 maxlim]);
        set(gca,'XLim',[0 maxlim]);
        axis off
        %title(['Con NStim-Pos. All Days'],'FontSize',16,'Fontweight','normal');
        
    
        figure(101); hold on;
        subplot(2,4,n); hold on;
        
        Expposnstim_norm = Expposnstim./max(max(Expposnstim));
        % Set all occupancy<thrsocc to just below 0. bounds will be from this no to 1. timestep=0.0334. 30pts~1 sec. 6pts=0.2sec
        Expposnstim_norm(find( Expposocc < 3))=-0.01;
        cmap = jet(1024); cmap(1,:) = 1; colormap(cmap);
        bounds = [-0.01,1];
        h = imagesc(Expposnstim_norm,bounds);
        %ch = colorbar;
        set(gca,'YLim',[0 maxlim]);
        set(gca,'XLim',[0 maxlim]);
        axis off
        %title(['Exp NStim-Pos. All Days'],'FontSize',16,'Fontweight','normal');
        mno=1;
        
        
        % These are for Dots
%         figure(50); hold on;
%         subplot(2,4,n); hold on;
%         plot(Conallposcurr(:,1),Conallposcurr(:,2),'Color',Clgy,'LineWidth',2);
%         plot(Conposstimcurr(:,1),Conposstimcurr(:,2),'bd','MarkerSize',4);
%         set(gca,'YLim',[0 110]);
%         set(gca,'XLim',[35 145]);
%         axis off
%         if n==1
%             title(['Con Pos at Stim-Pos. Days 1 to 8'],'FontSize',14,'Fontweight','normal');
%         end
%         ylabel('y-pos','FontSize',yfont,'Fontweight','normal');
%         xlabel('x-pos','FontSize',xfont,'Fontweight','normal');
%         figure(100); hold on;
%         subplot(2,4,n); hold on;
%         plot(Expallposcurr(:,1),Expallposcurr(:,2),'Color',Clgy,'LineWidth',2);
%         plot(Expposstimcurr(:,1),Expposstimcurr(:,2),'rd','MarkerSize',4);
%         set(gca,'YLim',[0 110]);
%         set(gca,'XLim',[35 145]);
%         axis off
%         if n==1
%             title(['Exp Pos at Stim-Pos. Days 1 to 8'],'FontSize',14,'Fontweight','normal');
%         end
  
    end % end day
    
    if savefig1==1,
        figure(51); hold on;
        figfile = [figdir,'RCb_PosAtStim_Map2'];
        print('-dpdf', figfile); print('-djpeg', figfile); saveas(gcf,figfile,'fig');
        
        figure(101); hold on;
        figfile = [figdir,'REd_PosAtStim_Map2'];
        print('-dpdf', figfile); print('-djpeg', figfile); saveas(gcf,figfile,'fig');
    end
    
    
end % dopossepdays


if dostimlearncorr==1
    
    % stimrate
    Exp_animstimrate = mean(Expstimrate_day,2);
    Con_animstimrate = mean(Constimrate_day,2);
    
    % learning trial
    Exp_lt = [55;197;92;46];
    Con_lt = [18;11;18;15];
    % learning day
    Exp_ld = [3;8;4;3];
    Con_ld = [1;1;2;1];
    
    [rt, pt] = corrcoef(Exp_animstimrate, Exp_lt); % r=-0.30, p=0.6992
    [rd, pd] = corrcoef(Exp_animstimrate, Exp_ld); % r=-0.33, p=0.6727
    [rtc, ptc] = corrcoef(Con_animstimrate, Con_lt); % r=0.47, p=0.5334
    [rdc, pdc] = corrcoef(Con_animstimrate, Con_ld); % r=-0.51, p=0.4856
    
    figure; hold on;
    redimscreen_2horsubplots;
    subplot(1,2,1); hold on;
    plot(Exp_animstimrate,Exp_lt,'ro','Markersize',12,'LineWidth',3); % r=-0.30, p=0.70, corr
    plot(Con_animstimrate,Con_lt,'bd','Markersize',12,'LineWidth',3); % r=0.47, p=0.53
    set(gca,'YLim',[0 201]);
    set(gca,'XLim',[0 2.1]);
    ylabel('Learning Trial','FontSize',yfont,'Fontweight','normal');
    xlabel('Mean Stimulation Rate (Hz)','FontSize',xfont,'Fontweight','normal');
    title('Learning Trial vs. Stim Rate');
    
    subplot(1,2,2); hold on;
    plot(Exp_animstimrate,Exp_ld,'ro','Markersize',12,'LineWidth',3); % r=-0.33, p=0.67
    plot(Con_animstimrate,Con_ld,'bd','Markersize',12,'LineWidth',3); % r=-0.51, p=0.49
    set(gca,'YLim',[0 8.1]);
    set(gca,'XLim',[0 2.1]);
    ylabel('Learning Day','FontSize',yfont,'Fontweight','normal');
    xlabel('Mean Stimulation Rate (Hz)','FontSize',xfont,'Fontweight','normal');
    title('Learning Day vs. Stim Rate');
    if savefig1==1,
        figfile = [figdir,'LearningRatevsStimRate'];
        print('-dpdf', figfile); print('-djpeg', figfile); saveas(gcf,figfile,'fig');
    end
    
    
end
    




% Position at Stimulation Quantification - After Review
% -----------------------------------------------------

if doposquant==1,
   
    
    %   ****************** Posn at Stim - Fraction stimln at wells and Ctr Int ******************
    % Take mean across all 8 days
    Conavgnstimwell_alldays = squeeze(mean(Conavgnstimwell,2));
    Expavgnstimwell_alldays = squeeze(mean(Expavgnstimwell,2));
    
    % Combine all 3 wells and maybe also plot all 3 wells separately
    Conavgnstimwell_alldays_allwells = sum(Conavgnstimwell_alldays,2); % For each animal
    Expavgnstimwell_alldays_allwells = sum(Expavgnstimwell_alldays,2); % For each animal
    

    % Plot All Wells Combined
    figure; hold on;
    if forppr==1
        redimscreen_figforppr1;
    else
        redimscreen_figforppt1;
    end
    jit=[0,0.05,0.1,0.15];
    plot(1+jit, Expavgnstimwell_alldays_allwells,'ro','MarkerSize',12,'Linewidth',3);
    plot(2+jit, Conavgnstimwell_alldays_allwells,'bd','MarkerSize',12,'Linewidth',3);
    [p_expcon_allwells,h_expcon_allwells] = ranksum(Conavgnstimwell_alldays_allwells,Expavgnstimwell_alldays_allwells); %p=1
    % p=0.8857
    axis([0.5 2.5 0 1]);
    set(gca,'XTick',[1 2],'XTickLabel',{'Exp';'Con'}','FontSize',xfont,'Fontweight','normal');
    title('Fraction of stimulations at well','FontSize',tfont,'Fontweight','normal');
    ylabel('Fraction of stimulations','FontSize',yfont,'Fontweight','normal');
    
    if savefig1==1,
        figfile = [figdir,'FractionStimAtWells'];
        print('-dpdf', figfile);
        print('-djpeg', figfile);
        saveas(gcf,figfile,'fig');
    end
    
    % Get Center Intersection - Center index is 1
    Conavgnstimint_alldays = squeeze(mean(Conavgnstimint,2));
    Expavgnstimint_alldays = squeeze(mean(Expavgnstimint,2));
    Conavgnstimint_alldays_ctr = Conavgnstimwell_alldays(:,1); % For each animal
    Expavgnstimint_alldays_ctr = Expavgnstimwell_alldays(:,1); % For each animal
    % Plot AllCtrIntersection 
    figure; hold on;
    if forppr==1
        redimscreen_figforppr1;
    else
        redimscreen_figforppt1;
    end
    jit=[0,0.05,0.1,0.15];
    plot(1+jit, Expavgnstimint_alldays_ctr,'ro','MarkerSize',12,'Linewidth',3);
    plot(2+jit, Conavgnstimint_alldays_ctr,'bd','MarkerSize',12,'Linewidth',3);
    [p_expcon_ctrint,h_expcon_ctrint] = ranksum(Conavgnstimint_alldays_ctr,Expavgnstimint_alldays_ctr); %p=1
    % p=1
    axis([0.5 2.5 0 1]);
    set(gca,'XTick',[1 2],'XTickLabel',{'Exp';'Con'}','FontSize',xfont,'Fontweight','normal');
    title('Fraction of stimulations at Ctr Int','FontSize',tfont,'Fontweight','normal');
    ylabel('Fraction of stimulations','FontSize',yfont,'Fontweight','normal');
    
    if savefig1==1,
       figfile = [figdir,'FractionStimAtCtrInt'];
       print('-dpdf', figfile);
       print('-djpeg', figfile);
       saveas(gcf,figfile,'fig');
    end
    
   
    %   ****************** Posn at Stim - vs Days ******************
    
    % Sum across wells
    Conavgnstimwell_days = squeeze(sum(Conavgnstimwell,3)); % nanim X ndays
    Expavgnstimwell_days = squeeze(sum(Expavgnstimwell,3));
    Conavgnstimwell_daymean = mean(Conavgnstimwell_days); Conavgnstimwell_dayerr = sem(Conavgnstimwell_days);
    Expavgnstimwell_daymean = mean(Expavgnstimwell_days); Expavgnstimwell_dayerr = sem(Expavgnstimwell_days);
    
    figure; hold on;
    if forppr==1
        redimscreen_figforppr1;
    else
        redimscreen_figforppt1;
    end
    errorbar(days,Expavgnstimwell_daymean,Expavgnstimwell_dayerr,'ro-','MarkerSize',12,'Linewidth',3);
    errorbar(days+0.1,Conavgnstimwell_daymean,Conavgnstimwell_dayerr,'bd-','MarkerSize',12,'Linewidth',3);
    [panovarm_stimrate] = anova_rm({Expavgnstimwell_days Conavgnstimwell_days},'off');
    %pgrp=0.6783, ptime=0.00, pint=0.3832
    
    % Make Plot presentable
    title('Fraction Stim At Wells vs Days','FontSize',tfont,'Fontweight','normal');
    xlabel('Day','FontSize',xfont,'Fontweight','normal');
    ylabel('Fraction Stim at Wells','FontSize',yfont,'Fontweight','normal');
    %legend('Location','NorthEast');
    axis([0.5 max(days)+0.5 0 1]); 
    if savefig1==1,
        figfile = [figdir,'FractionStimAtWellsVsDays'];
        print('-dpdf', figfile);
        print('-djpeg', figfile);
        saveas(gcf,figfile,'fig');
    end
    
    
    % For Ctr Arm - Center index is 1
    Conavgnstimint_ctrdays = squeeze(Conavgnstimint(:,:,1));
    Expavgnstimint_ctrdays = squeeze(Expavgnstimint(:,:,1));
    Conavgnstimint_ctrdaymean = mean(Conavgnstimint_ctrdays); Conavgnstimint_ctrdayerr = sem(Conavgnstimint_ctrdays);
    Expavgnstimint_ctrdaymean = mean(Expavgnstimint_ctrdays); Expavgnstimint_ctrdayerr = sem(Expavgnstimint_ctrdays);

    figure; hold on;
    if forppr==1
        redimscreen_figforppr1;
    else
        redimscreen_figforppt1;
    end
    errorbar(days,Expavgnstimint_ctrdaymean,Expavgnstimint_ctrdayerr,'ro-','MarkerSize',12,'Linewidth',3);
    errorbar(days+0.1,Conavgnstimint_ctrdaymean,Conavgnstimint_ctrdayerr,'bd-','MarkerSize',12,'Linewidth',3);
    [panovarm_ctr] = anova_rm({Expavgnstimint_ctrdays Conavgnstimint_ctrdays},'off');
    %pgrp=0.8182, ptime=0.0173, pint=0.5934
    
    % Make Plot presentable
    title('Fraction Stim At Ctr Int vs Days','FontSize',tfont,'Fontweight','normal');
    xlabel('Day','FontSize',xfont,'Fontweight','normal');
    ylabel('Fraction Stim at Ctr Int','FontSize',yfont,'Fontweight','normal');
    %legend('Location','NorthEast');
    axis([0.5 max(days)+0.5 0 0.3]); 
    if savefig1==1,
        figfile = [figdir,'FractionStimAtCtrIntVsDays'];
        print('-dpdf', figfile);
        print('-djpeg', figfile);
        saveas(gcf,figfile,'fig');
    end
    
    
end


% Speed Before After Stim - After Review
% -----------------------------------------------------

if dovelstim_befaft==1,
       
    % Plot mean and error bars Before and After stimulation 
    Conmean_befaft = mean(Convelstim_befaft_alldays);  % n = 2630
    Conerr_befaft = sem(Convelstim_befaft_alldays);
    Expmean_befaft = mean(Expvelstim_befaft_alldays); % n = 4978 
    Experr_befaft = sem(Expvelstim_befaft_alldays);
    
    % Vectors of histogram - for stats
    Convec_alldays = Convelstim_befaft_alldays(:);
    Expvec_alldays = Convelstim_befaft_alldays(:);
    
    xhist = [-5:0.1:4.9];
    figure; hold on;
    if forppr==1
        redimscreen_figforppr1;
    else
        redimscreen_figforppt1;
    end
    
    plot(xhist,Conmean_befaft,'b-','Linewidth',3,'MarkerSize',2);
    plot(xhist,Expmean_befaft,'r-','Linewidth',3,'MarkerSize',2);
    jbfill(xhist, [Conmean_befaft+Conerr_befaft],[Conmean_befaft-Conerr_befaft],'b','b',1,0.3);
    jbfill(xhist, [Expmean_befaft+Experr_befaft],[Expmean_befaft-Experr_befaft],'r','r',1,0.3);
    plot(0*ones(size([0:1:15])),[0:1:15],'k--','LineWidth',2);
    plot(-1*ones(size([0:1:15])),[0:1:15],'k--','LineWidth',1);
    plot(1*ones(size([0:1:15])),[0:1:15],'k--','LineWidth',1);
    plot(-2*ones(size([0:1:15])),[0:1:15],'k--','LineWidth',1);
    plot(2*ones(size([0:1:15])),[0:1:15],'k--','LineWidth',1);
    set(gca,'YLim',[0 15]); 
    set(gca,'XLim',[-5 +5]);
    title(['Con and Exp - Speed bef aft Stim. All Days'],'FontSize',16,'Fontweight','normal');
    ylabel('Speed (cm/sec)','FontSize',yfont,'Fontweight','normal');
    xlabel('Time from stimulation(sec)','FontSize',xfont,'Fontweight','normal');
    % ----------- Stats -------------
    %[h_vel,p_vel] = ttest2(Convel_alldays,Expvel_alldays);
    %[h_velhist,p_velhist] = kstest2(Convelhistvec_alldays,Expvelhistvec_alldays);

    % Now plot the curves with mean removed
    figure; hold on;
    if forppr==1
        redimscreen_figforppr1;
    else
        redimscreen_figforppt1;
    end
    Conmean_befaftnorm = Conmean_befaft - mean(Conmean_befaft);
    Expmean_befaftnorm = Expmean_befaft - mean(Expmean_befaft);
    plot(xhist,Conmean_befaftnorm,'b-','Linewidth',3,'MarkerSize',2);
    plot(xhist,Expmean_befaftnorm,'g-','Linewidth',3,'MarkerSize',2);
    jbfill(xhist, [Conmean_befaftnorm+Conerr_befaft],[Conmean_befaftnorm-Conerr_befaft],'b','b',1,1);
    jbfill(xhist, [Expmean_befaftnorm+Experr_befaft],[Expmean_befaftnorm-Experr_befaft],'g','g',1,1);
    
    CEmean = Expmean_befaftnorm;
    CEerr = Experr_befaft;
   
    plot(xhist,CEmean,'k-','Linewidth',3,'MarkerSize',2);
    jbfill(xhist, [CEmean+CEerr],[CEmean-CEerr],'k','k',1,1);
    
    plot(0*ones(size([-5:1:5])),[-5:1:5],'k--','LineWidth',2);
    plot(-1*ones(size([-5:1:5])),[-5:1:5],'k--','LineWidth',1);
    plot(1*ones(size([-5:1:5])),[-5:1:5],'k--','LineWidth',1);
    plot(-2*ones(size([-5:1:5])),[-5:1:5],'k--','LineWidth',1);
    plot(2*ones(size([-5:1:5])),[-5:1:5],'k--','LineWidth',1);
    set(gca,'YLim',[-3 3]); 
    set(gca,'XLim',[-5.1 5.2]);
    title(['Con and Exp - Norm Speed bef aft Stim. All Days'],'FontSize',16,'Fontweight','normal');
    ylabel('Normalized Speed (cm/sec)','FontSize',yfont,'Fontweight','normal');
    xlabel('Time from stimulation (sec)','FontSize',xfont,'Fontweight','normal');
    
    if savefig1==1,
        figfile = [figdir,'velstim_befaft_withnostim_axis2'];
        print('-dpdf', figfile);
        print('-djpeg', figfile);
        saveas(gcf,figfile,'fig');
    end
    
    
    % Difference in speeds Before and after Stimln - 1 sec window
    figure; hold on;
    if forppr==1
        redimscreen_figforppr1;
    else
        redimscreen_figforppt1;
    end
    bar(1,mean(Expvelstim_diff_alldays),'r');
    bar(2,mean(Convelstim_diff_alldays),'b');
    errorbar(1,mean(Expvelstim_diff_alldays),sem(Expvelstim_diff_alldays),'k.-','LineWidth',2);
    errorbar(2,mean(Convelstim_diff_alldays),sem(Convelstim_diff_alldays),'k.-','LineWidth',2);
    
    set(gca,'YLim',[-3 1]); 
    set(gca,'XTick',[1 2]);
    ylabel('Speed Difference (cm/sec)','FontSize',yfont,'Fontweight','normal');
     
    if savefig1==1,
        figfile = [figdir,'velstim_diff_2swin'];
        print('-dpdf', figfile);
        print('-djpeg', figfile);
        saveas(gcf,figfile,'fig');
    end
    
    % Plot for 2 sec window also
    figure; hold on;
    if forppr==1
        redimscreen_figforppr1;
    else
        redimscreen_figforppt1;
    end
    bar(1,mean(Expvelstim_diff2_alldays),'r');
    bar(2,mean(Convelstim_diff2_alldays),'b');
    errorbar(1,mean(Expvelstim_diff2_alldays),sem(Expvelstim_diff2_alldays),'k','LineWidth',2);
    errorbar(2,mean(Convelstim_diff2_alldays),sem(Convelstim_diff2_alldays),'k','LineWidth',2);
    
    set(gca,'YLim',[-3 1]); 
    set(gca,'XTick',[1 2]);
    
    % Plot the Diff By Anim, Con vs Exp grps
    % ---------------------------------------
    figure; hold on;
    if forppr==1
        redimscreen_figforppr1;
    else
        redimscreen_figforppt1;
    end
    
    for a=1:length(Expvelstim_diff_byanim)
        Expmean_velstim_diff_byanim(a) = mean(Expvelstim_diff_byanim{a});
        Experr_velstim_diff_byanim(a) = sem(Expvelstim_diff_byanim{a});
    end
    for a=1:length(Convelstim_diff_byanim)
        Conmean_velstim_diff_byanim(a) = mean(Convelstim_diff_byanim{a});
        Conerr_velstim_diff_byanim(a) = sem(Convelstim_diff_byanim{a});
    end
    jit=[0,0.05,0.1,0.15];
    plot(3+jit(1:length(Expvelstim_diff_byanim)), Expmean_velstim_diff_byanim+0.1,'ks','MarkerSize',12,'Linewidth',3);
    plot(2+jit(1:length(Convelstim_diff_byanim)), Conmean_velstim_diff_byanim,'bd','MarkerSize',12,'Linewidth',3);
    [p_velstim_diff,h_velstim_diff] = ranksum(Expmean_velstim_diff_byanim,Conmean_velstim_diff_byanim); %p=
    axis([0.5 3.5 -2 0.5]);
    set(gca,'XTick',[1 2],'XTickLabel',{'Exp';'Con'}','FontSize',xfont,'Fontweight','normal');
    title('Speed Difference at Stim - 1 sec window','FontSize',tfont,'Fontweight','normal');
    ylabel('Speed difference','FontSize',yfont,'Fontweight','normal');
    
    if savefig1==1,
       figfile = [figdir,'velstim_diff_byanim_1swin_withnostim'];
       print('-dpdf', figfile);
       print('-djpeg', figfile);
       saveas(gcf,figfile,'fig');
    end
    
end



    
keyboard;

   
   
















