
% Shantanu - 17Apr2012. From DFSsj_stimstats_summ.
% Reducing the plotting code to just the Position at Stimulation, for David, etc.


% From sj_stimstatssumm1. I want to broaden it to easily plot distributions,
% eg. speed at stimuln distribution, etc

% sj_stimstatssumm1 and this script gets data saved by the analysis function sj_stim_behstats.m
% which saves a file called prefix_stimstats in each "prefixdirect" directory
% Instead of loading the file,
% Run the analysis from sj_stim_behstats again here, as DFAsj_stim_behstats
% Analysis function would load "pos","linpos" and "DIO" structure, and would use the iterator "epoch_behaveanal"

clear; close all;
runscript = 0; % if is 0, just loads saved processed data and plots stuff
savedata = 0; % save data option - only works if runscript is also on
figopt1=0; % Figure Options - main summ figs
dodistr=0; % Figures for speed distributions
dosepdays=0; % Separate figures for speed distributions on each day
dostimisi=0; % Figures for stimISI
doisisepdays=0;% Separate figures for stim isi on each day
doposdistr=0; % Figures for pos distr at stim - all days - If combining posn data across animals
dostimlearncorr=0; % Stimrate vs learning rate
dovelstim_befaft=1; % Vel before after stimln - after review

dopossepdays=1; % Figures for position distr on each day - Can combine across animals for each day, or use to plot for 1 animal
doposquant=1; % Position at stimulation quantification - after review

% For binning and smoothing position
binsize = 2; % cm square
stdev = 2;
timestep = 1/29.97;
threshocc = 0.02; % Threshold occupancy in seconds

savedir = '/data25/sjadhav/RippleInterruption/ProcessedData/';
%savefile = [savedir 'Stimstats_summ'];
%savefile = [savedir 'Stimstats_summ_rev']; % after revision
savefile = [savedir 'Stimstats_summ_rev2']; %revision2 - more animals in Exp group
% If runscript, run Datafilter and save data
if runscript == 1
    
    %Animal selection
    %-----------------------------------------------------
    Expanimals = {'REc','REd','REe','REf','REg','REh'};
    Conanimals = {'RCa','RCb','RCc','RCd'};
    
    %Filter creation
    %--------------------------------------------------------
    % epoch filter
    dayfilter = '1:8'; % Shantanu - I am adding day filter to parse out epoch filter
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

figdir = '/data25/sjadhav/RippleInterruption/Figures/01AugSep11_RippleDisFigs/StimulationStats/SpeedAtStim/'; %/SpeedAtStim/
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
                
                
                % Posn at Stimln
                %eval([str(g,:),'wellpos{an}{d}{ep} =',str(g,:),'stimf(an).output{1}(index).wellpos;']);
                %eval([str(g,:),'intpos{an}{d}{ep} =',str(g,:),'stimf(an).output{1}(index).intpos;']);
                %eval([str(g,:),'nstimwell(an,d,ep,:) =',str(g,:),'stimf(an).output{1}(index).nstimwell;']); % 3 values for each epoch
                %eval([str(g,:),'nstimint(an,d,ep,:) =',str(g,:),'stimf(an).output{1}(index).nstimint;']); % 3 values for each epoch
                %eval([str(g,:),'nstimarm(an,d,ep,:) =',str(g,:),'stimf(an).output{1}(index).nstimarm;']); % 3 values for each epoch
                %eval([str(g,:),'pstimarm{an}{d}{ep} =',str(g,:),'stimf(an).output{1}(index).pstimarm;']); % 3 values for each epoch
                
            end % end epoch
            
            % Combine posn matr for day across epochs
            eval([str(g,:),'posstim_day{an}{d} =[',str(g,:),'posstim{an}{d}{1};',str(g,:),'posstim{an}{d}{2}];']); % nX2 matrix of x-y positions
            eval([str(g,:),'pos_day{an}{d} =[',str(g,:),'pos{an}{d}{1};',str(g,:),'pos{an}{d}{2}];']); % nX2 matrix of x-y positions
            
            %              % Take avg "posn at stimln" for the two epochs, AND GET TOTAL STIM
            %             eval([str(g,:),'avgnstimwell(an,d,:) = mean(squeeze([',str(g,:),'nstimwell(an,d,1,:);',str(g,:),'nstimwell(an,d,2,:)]));']);
            %             eval([str(g,:),'avgnstimint(an,d,:) = mean(squeeze([',str(g,:),'nstimint(an,d,1,:);',str(g,:),'nstimint(an,d,2,:)]));']);
            %             eval([str(g,:),'avgnstimarm(an,d,:) = mean(squeeze([',str(g,:),'nstimarm(an,d,1,:);',str(g,:),'nstimarm(an,d,2,:)]));']); % Something wrong here for arm
            %             eval([str(g,:),'avgtotalstim(an,d) = mean([',str(g,:),'totalstim(an,d,1);',str(g,:),'totalstim(an,d,2)]);']);
            %             % Get in terms of ratio by dividing by total stim for each day
            %             eval([str(g,:),'avgnstimwell(an,d,:) =',str(g,:),'avgnstimwell(an,d,:)./',str(g,:),'avgtotalstim(an,d);']);
            %             eval([str(g,:),'avgnstimint(an,d,:) =',str(g,:),'avgnstimint(an,d,:)./',str(g,:),'avgtotalstim(an,d);']);
            
        end % end day
    end % end anim
    % Average across epochs for parameter vs. day plots
    eval([str(g,:),'stimrate_day = nanmean(',str(g,:),'stimrate,3);']);
    eval([str(g,:),'totalstim_day = nanmean(',str(g,:),'totalstim,3);']);
    eval([str(g,:),'totaltime_day = nanmean(',str(g,:),'totaltime,3);']);
    
    
    % Position matrices. Can I Combine across animals for REd,REe,REf and RCb,RCc,RCd - same tracks
    % --------------------
    % Initialize to store for alldays
    eval([str(g,:),'pos_alldays=[];']);
    eval([str(g,:),'posstim_alldays=[];']);
    % Loop over each day and get separately
    for d=1:length(days)
        currday=days(d);
        % Initialize to store raw values for this day across animals for stats
        eval([str(g,:),'allpos_day{',num2str(currday),'}=[];']);
        eval([str(g,:),'allposstim_day{',num2str(currday),'}=[];']);
        
        % Get from each animal
        % for ani=2:totanim % Loop over anim
        % an=ani-1;
        for ani=2:2 % Loop over anim
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
    
    
    % *************************************
    % Position Distributions - Combined for all days
    % *************************************
    
    if binsize==1, maxlim = 115; end
    if binsize==2, maxlim = 57; end
    
    
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
        
    end % doposdistr
    
    
    
end











