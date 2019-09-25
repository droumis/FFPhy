

%Filter for rip modulated or unmodulated PFC cells then call getspatialinfo_DR
warning('off','all');

clear;
runscript = 0;
savedata = 0; % save data option - only works if runscript is also on
% figopt1 = 1; % Figure Options - Individual cells
plotstuff = 1; % needs to be 1 if you want to plot anything.. the other plot params are not mutually exlusive
ploteachan = 1; % plot figs for each animal
plotcomban = 1; %plot figs for comb animals
plotoverdays = 1; %plot spat info over days.. working on this now
modUnmod = [1 2]; % [1 2] if you want both mod and unmod. 1 is mod only. 2 if unmod only. i suggest using both
whichep = 2; % 2 is last w ep, 1 is first w ep, 0 is all w eps % Shantanu!! need to change this when you combine
savedir = '/mnt/data25/sjadhav/HPExpt/ProcessedDataDR/';
% savefile = [savedir 'HP_ripmod_PFC']; area = 'PFC'; clr = 'b'; % PFC
% savefile = [savedir 'HP_ripmodUnmod_PFC']; area = 'PFC'; %clrunmod = 'r'; clrmod = 'b'; % PFC  %%%% make this into a for loop to collect mod and unmod
savefile = [savedir 'testcombine']; area = 'PFC'; %clrunmod = 'r'; clrmod = 'b'; % PFC  %%%% make this into a for loop to collect mod and unmod

savefig1=0;

% Plot options
plotanimidx =  []; % To pick animals for plotting
plotdays = []; % If you only load data when runscript=0 and savedata=0, then this field will supplant days

%params from AS.. which do I need?
Veqn = '>3';
minV=str2num(Veqn(end));
% minVPF = 3 %cm/sec
% minPeakPF = 3
mintime = 3;
% n = 1
% inclStates = 6
% traj = [3 4]
%wtrack = 0; %
traj = [1:4] ;%[1:4]; %[1 4; 2 3];   %%% is this the right way to do it? hmm. Yes, it is.
%armlength = 66 %70 Wtrack %66 6 arm maze;
% maxstage = [12 3]
% comp = 1
% correct = 0
% excluderepeat = 1
% fnc = 3


% If runscript, run Datafilter and save data
if runscript == 1
    for i =  modUnmod;
        %Animal selection
        %-----------------------------------------------------
                animals = {'HPa' 'HPb' 'HPc'};
        %         animals = {'HPc'};
%         animals = {'HPb'};
        
        
        %Filter creation
        %-----------------------------------------------------
        
        % Epoch filter
        % -------------
        dayfilter = '1:8';
        % Either Only do 1st w-track. 2 or 1 epochs per day
        % Or do Wtr1 and Wtr2, 2 epochs per day
        runepochfilter = 'isequal($environment, ''wtr1'') || isequal($environment, ''wtr2'')';
        % sleepepochfilter = 'isequal($type, ''sleep'')'; % Only pre and post sleep marked as sleep
        % sleepepochfilter = 'isequal($environment, ''postsleep'')'; % Only pre and post sleep marked as sleep
        
        % %Cell filter
        % %-----------
        if i < 2;
            placecellfilter = '(strcmp($area, ''PFC'') && strcmp($ripmodtag, ''y''))';   % Ripple mod
        else
            placecellfilter = '(strcmp($area, ''PFC'') && strcmp($ripmodtag, ''n''))'; % Ripple unmod
        end
        
        % Time filter -
        %%-----------
        % I i think this its a good idea to exclude close-to-well firing as
        % well as very low velocity times..
        %     timefilter =  {{'getdistclosestwell', '($distwell >= 10)'} {'getlinvelocity', ['((abs($velocity) ',Veqn,'))'] } };
        
        % exclude ripple firing !! why didn't annabelle have an exlcude rip
        % filter? i think bc she had a velocity filter above which may have
        % been sufficient.
        riptetfilter = '(isequal($descrip, ''riptet''))';
        %     timefilter = {{'DFTFsj_getriptimes','($nripples == 0)','tetfilter',riptetfilter, 'minthresh',3 }};
        timefilter = { {'DFTFsj_getriptimes','($nripples == 0)','tetfilter',riptetfilter,'minthresh',3} };
        
        %     timefilter{length(timefilter)+1} = {'DFTFsj_getriptimes','($nripples == 0)','tetfilter',riptetfilter,'minthresh',3};
        
        
        % extra timefilters from AS
        %%--------------------------
        % i think she's looking at correct vs incorrect... i don't need this
        % right now
        %     if correct == 1 && wtrack == 0
        %         timefilter{length(timefilter)+1} = {'filtergetcorrecttraj', '($state > 0)', 6, minV, 0, 10, 150} ;
        %     elseif correct == 1 && wtrack == 1
        %         timefilter{length(timefilter)+1} = {'filtergetcorrecttraj', '($state > 0)', 6, minV, 0, 10, 150, 'correctorder', [2 1 3]} ;
        %     end
        
        
        %  I don't think i need this... but not positive what it does..
        %     if i < 10
        %         timefilter{length(timefilter)+1} = {'getcalctaskstage', ['($includebehave == ', num2str(i), ')'], comp} ;
        %     elseif i == 12
        %         timefilter{length(timefilter)+1} = {'getcalctaskstage', ['($includebehave == 1) | ($includebehave == 2)'], comp} ;
        %     end
        
        
        % Iterator
        % --------
        iterator = 'singlecellanal';
        
        % Filter creation
        % ----------------
        spatf = createfilter('animal',animals,'days',dayfilter,'epochs',runepochfilter,'cells',placecellfilter,'excludetime', timefilter,'iterator',iterator);
        
        %do i need this?
        spatf = testexcludetimes(spatf, mintime); %removes epochs from analysis if all epoch excluded by excludetimes, mintime = 30
        
        disp('Done Filter Creation');
        
        % Set analysis function
        % ----------------------
        
        %spatf = setfilterfunction(spatf, 'getrawdata_DR', {'linpos', 'spikes'},6,  minV); % without any varargins. use for rawdata collect
        spatf = setfilterfunction(spatf, 'getspatialinfo_DR', {'linpos', 'spikes'},6,  minV, 'peakrate', 3,'appendindex', 1, 'incltraj', traj); %use this one for non combined
        
        
        
        % Run analysis
        % ------------ %% create data struct for both mod, unmod if specified
        % above.
        spatinfo{i} = runfilter(spatf);
    end
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


% ---------------------------------
% Combining output across epochs 
% ---------------------------------

% First, gather data
cnt=0; % Count how many cells will be kept 
allanimindex=[]; allspatinfo=[];
 for an = 1:length(spatinfo)
        for i=1:length(spatinfo(an).output{1})
            % Check for empty output - If Cell defined in epoch and Nspks in ripple response wndow > 0
            if (~isempty(spatinfo(an).output{1}(i).spatinfo_traj))
                cnt = cnt+1;
                anim_index{an}(cnt,:) = spatinfo(an).output{1}(i).index; % Animals kept separate
                % Only indexes
                animindex=[an spatinfo(an).output{1}(i).index]; % Put animal index in front
                allanimindex = [allanimindex; animindex]; % Collect all Anim Day Epoch Tet Cell Index
                % Data
                allspatinfo(:,cnt) = spatinfo(an).output{1}(i).spatinfo_traj; % Only get raster and histogram response
            end
        end
 end
 
% Consolidate single cells across epochs.
% -------------------------------- 

% A) Here, all animals are collapsed together. Can also keep animals seperate
% ----------------------------------------------------------------------

 dummyindex=allanimindex;  % all anim-day-epoch-tet-cell indices
 dayspatinfo_idx=[]; dayspatinfo=[];
 for i=1:size(allspatinfo,1)
     animdaytetcell=allanimindex(i,[1 2 4 5]);
     ind=[];
     while rowfind(animdaytetcell,dummyindex(:,[1 2 4 5]))~=0          % collect all rows (epochs)
         ind = [ind rowfind(animdaytetcell,dummyindex(:,[1 2 4 5]))];        % finds the first matching row
         dummyindex(rowfind(animdaytetcell,dummyindex(:,[1 2 4 5])),:)=[0 0 0 0 0]; % after adding index, remove the corresponding row
         % so you could find the next one
     end
     
     % Gather everything for the current cell across epochs
     currspatinfo=[];
     for r=ind
         currspatinfo = [currspatinfo, currspatinfo(:,r)];
     end    
     
     % Save as variables or in a structure
     dayspatinfo_idx(i,:) = animdaytetcell;
     dayspatinfo(:,i) = mean(currspatinfo,2); % Mean across epochs for all trajs
     
     dayspatinfo_struct(i).index = animdaytetcell;
     dayspatinfo_struct(i).spatinfo_traj = mean(currspatinfo,2); % Mean across epochs for all trajs     
 end
 
% B) Here, animals are kept seperate
% -----------------------------------             
for an = 1:length(animals)
    dummyindex = anim_index{an};     % collect all day-epoch-tet-cell indices
    for i=1:length(spatinfo(an).output{1})
        daytetcell=spatinfo.output{1}(i).index([1 3 4]);
        ind=[];
        while rowfind(daytetcell,dummyindex(:,[1 3 4]))~=0          % collect all rows (epochs)
            ind = [ind rowfind(daytetcell,dummyindex(:,[1 3 4]))];
            dummyindex(rowfind(daytetcell,dummyindex(:,[1 3 4])),:)=[0 0 0 0];
        end
        % Gather everything for the current cell across epochs
        currspatinfo=[];
        for r=ind
            currspatinfo = [currspatinfo, spatinfo.output{1}(r).spatinfo_traj];
        end
        
        % Save as variables or in a structure
        animdayspatinfo_idx(an,i,:) = daytetcell;
        animdayspatinfo(an,i,:) = daytetcell; % Mean across epochs for all trajs
        
        animdayspatinfo_struct(an).index(i,:) = daytetcell;
        animdayspatinfo_struct(an).spatinfo_traj(:,i) = mean(currspatinfo,2); % Mean across epochs for all trajs
    end
end
                




% ----------------------------------
% PLOT!
% --------------------------------------------------------------------



% plot( each animal) do this first//
if plotstuff ==1
    close all
    if ploteachan == 1
        
        for i = 1: length(animals);
            spatmoddata = (spatinfo{1,1}(1,i).output{1}); %houskeeping
            spatunmoddata = (spatinfo{1,2}(1,i).output{1});
            ndays = length(unique(spatinfo{1}(1,i).epochs{1}(:,1)));
            %use last, first, or all epochs (all epochs means treating the
            %same cells across epochs within a day as different cells,
            %which will double the n of cells = bad
            if whichep == 2 %use last ep
                tit = 'Last W Epoch';
                mod = spatmoddata(find((spatmoddata(:,1)==1 & spatmoddata(:,2)>5) | (spatmoddata(:,1)>1 & spatmoddata(:,2)>3)),6); %the last run ep
                unmod = spatunmoddata(find((spatunmoddata(:,1)==1 & spatunmoddata(:,2)>5) | (spatunmoddata(:,1)>1 & spatunmoddata(:,2)>3)),6); %the last run ep
                
                %day chunking
                daymodall = [spatmoddata(find((spatmoddata(:,1)==1 & spatmoddata(:,2)>5) | (spatmoddata(:,1)>1 & spatmoddata(:,2)>3)),1)...
                    spatmoddata(find((spatmoddata(:,1)==1 & spatmoddata(:,2)>5) | (spatmoddata(:,1)>1 & spatmoddata(:,2)>3)),6)];
                dayunmodall = [spatunmoddata(find((spatunmoddata(:,1)==1 & spatunmoddata(:,2)>5) | (spatunmoddata(:,1)>1 & spatunmoddata(:,2)>3)),1)...
                    spatunmoddata(find((spatunmoddata(:,1)==1 & spatunmoddata(:,2)>5) | (spatunmoddata(:,1)>1 & spatunmoddata(:,2)>3)),6)];
                for j = 1:ndays;
                    daymod{j} = daymodall(find(daymodall(:,1)==j),2);
                    dayunmod{j} = dayunmodall(find(dayunmodall(:,1)==j),2);
                    if j == 1;
                        meanday = [mean(daymodall(find(daymodall(:,1)==j),2)) mean(dayunmodall(find(dayunmodall(:,1)==j),2))]
                    else
                        meanday = [meanday; [mean(daymodall(find(daymodall(:,1)==j),2)) mean(dayunmodall(find(dayunmodall(:,1)==j),2))]];
                    end
                end
                
                
            elseif whichep == 1 %use first ep
                tit = 'First W Epoch';
                mod = spatmoddata(find((spatmoddata(:,1)==1 & spatmoddata(:,2)<5) | (spatmoddata(:,1)>1 & spatmoddata(:,2)<3)),6);
                unmod = spatunmoddata(find((spatunmoddata(:,1)==1 & spatunmoddata(:,2)>5) | (spatunmoddata(:,1)>1 & spatunmoddata(:,2)>3)),6);
               
                %day chunking
                daymodall = [spatmoddata(find((spatmoddata(:,1)==1 & spatmoddata(:,2)<5) | (spatmoddata(:,1)>1 & spatmoddata(:,2)<3)),1)...
                    spatmoddata(find((spatmoddata(:,1)==1 & spatmoddata(:,2)<5) | (spatmoddata(:,1)>1 & spatmoddata(:,2)<3)),6)];
                dayunmodall = [spatunmoddata(find((spatunmoddata(:,1)==1 & spatunmoddata(:,2)<5) | (spatunmoddata(:,1)>1 & spatunmoddata(:,2)<3)),1)...
                    spatunmoddata(find((spatunmoddata(:,1)==1 & spatunmoddata(:,2)<5) | (spatunmoddata(:,1)>1 & spatunmoddata(:,2)<3)),6)];
                for j = 1:ndays;
                    daymod{j} = daymodall(find(daymodall(:,1)==j),2); % i think I'll need these for errorbars                   
                    dayunmod{j} = dayunmodall(find(dayunmodall(:,1)==j),2);
                    if j == 1;
                        meanday = [mean(daymodall(find(daymodall(:,1)==j),2)) mean(dayunmodall(find(dayunmodall(:,1)==j),2))]
                    else
                        meanday = [meanday; [mean(daymodall(find(daymodall(:,1)==j),2)) mean(dayunmodall(find(dayunmodall(:,1)==j),2))]];
                    end
                end
                
            else %use all eps
                tit = 'All W Epochs (not combined)';
                mod = spatmoddata(:,6); %spatinfo{1,1}(1,i).output{1}(:,6);
                unmod = spatunmoddata(:,6); %spatinfo{1,2}(1,i).output{1}(:,6);
                
                %day chunking
                for j = 1:ndays;
                    daymod{j} = spatmoddata(find(spatmoddata(:,1)==j),6);
                    dayunmod{j} = spatunmoddata(find(spatunmoddata(:,1)==j),6);
                    if j == 1;
                        meanday = [mean(spatmoddata(find(spatmoddata(:,1)==j),2)) mean(spatunmoddata(find(spatunmoddata(:,1)==j),2))]
                    else
                        meanday = [meanday; [mean(spatmoddata(find(spatmoddata(:,1)==j),2)) mean(spatunmoddata(find(spatunmoddata(:,1)==j),2))]];
                    end
                end
                
            end
            edges = linspace(min([mod; unmod]),max([mod; unmod]),20); %%
            modcnt = histc(mod,edges);
            unmodcnt = histc(unmod,edges);
            figure(i)
            hold on
            subplot(1,2,1)
            plot(edges, modcnt/sum(modcnt), 'Color','b', 'linewidth', 3);
            hold on
            plot(edges, unmodcnt/sum(unmodcnt), 'Color','r', 'linewidth', 3);
            title({sprintf('%s %s Rip mod/unmod spatial info', animals{1,i}, area); tit; sprintf('n cells mod: %d unmod: %d',length(mod),length(unmod))})
            xlabel('info (bits/spike)')
            ylabel('fraction of total neurons')
            legend('mod', 'unmod')
            
            %stats test
            [h p] = kstest2(mod(~isnan(mod)),unmod(~isnan(unmod)));
            medmod = median(mod(~isnan(mod)));
            medunmod = median(unmod(~isnan(unmod)));
            sizes = [size(mod(~isnan(mod)),1) size(unmod(~isnan(unmod)),1) ];
            mod = mod(~isnan(mod));
            unmod = unmod(~isnan(unmod));
            subplot(1,2,2);
            hold on
            bar([1 2], [mean(mod) mean(unmod)], 0.8, 'EdgeColor', 'none')
            errorbar2([1 2], [mean(mod) mean(unmod)],  [stderr(mod) stderr(unmod)] , 0.3, 'k')
            xlim([0.3 2.7])
            set(gca, 'fontsize', 24)
            set(gca, 'xtick', [1 2], 'xticklabel', {'mod', 'unmod'})
            ylabel('info (bits/spike)')
            rp = ranksum(mod,unmod);
            title({animals{1,i}; sprintf('ranksum p = %s',rp); sprintf('K-Stest p = %s',p)})
            
            %plot info over days %%%%%%in progress
%             if plotoverdays == 1;
%                 subplot(2,2,3);
%                 hold on
%                 bar(1:ndays, meanday, 'EdgeColor','none')
%             end
            
            %plot info novel1(d1,2) vs familiar(d3-5) vs novel2(d6-..)
            
            
            %plot Wt1 vs Wt2
            
            
            %plot mosaic of rate field maps
        end
        
        
    end
    
    if plotcomban == 1
        for i = 1: length(animals);
            spatmoddata = (spatinfo{1,1}(1,i).output{1}); %houskeeping
            spatunmoddata = (spatinfo{1,2}(1,i).output{1});
            if i ==1
                if whichep == 2 %use last ep
                    tit = 'Last W Epoch';
                    mod = spatmoddata(find((spatmoddata(:,1)==1 & spatmoddata(:,2)>5) | (spatmoddata(:,1)>1 & spatmoddata(:,2)>3)),6);
                    unmod = spatunmoddata(find((spatunmoddata(:,1)==1 & spatunmoddata(:,2)>5) | (spatunmoddata(:,1)>1 & spatunmoddata(:,2)>3)),6);
                elseif whichep == 1 %use first ep
                    tit = 'First W Epoch';
                    mod = spatmoddata(find((spatmoddata(:,1)==1 & spatmoddata(:,2)<5) | (spatmoddata(:,1)>1 & spatmoddata(:,2)<3)),6);
                    unmod = spatunmoddata(find((spatunmoddata(:,1)==1 & spatunmoddata(:,2)>5) | (spatunmoddata(:,1)>1 & spatunmoddata(:,2)>3)),6);
                else %use all eps
                    tit = 'All W Epochs';
                    mod = spatmoddata(:,6); %spatinfo{1,1}(1,i).output{1}(:,6);
                    unmod = spatunmoddata(:,6); %spatinfo{1,2}(1,i).output{1}(:,6);
                end
            else
                if whichep == 2 %use last ep
                    mod = stack(mod,spatmoddata(find((spatmoddata(:,1)==1 & spatmoddata(:,2)>5) | (spatmoddata(:,1)>1 & spatmoddata(:,2)>3)),6));
                    unmod = stack(unmod,spatunmoddata(find((spatunmoddata(:,1)==1 & spatunmoddata(:,2)>5) | (spatunmoddata(:,1)>1 & spatunmoddata(:,2)>3)),6));
                elseif whichep == 1 %use first ep
                    mod = stack(mod,spatmoddata(find((spatmoddata(:,1)==1 & spatmoddata(:,2)<5) | (spatmoddata(:,1)>1 & spatmoddata(:,2)<3)),6));
                    unmod = stack(unmod,spatunmoddata(find((spatunmoddata(:,1)==1 & spatunmoddata(:,2)>5) | (spatunmoddata(:,1)>1 & spatunmoddata(:,2)>3)),6));
                else %use all eps
                    mod = stack(mod,spatmoddata(:,6)); %spatinfo{1,1}(1,i).output{1}(:,6);
                    unmod = stack(unmod,spatunmoddata(:,6)); %spatinfo{1,2}(1,i).output{1}(:,6);
                    %                     mod = stack(mod, spatinfo{1,1}(1,i).output{1}(:,6));
                    %                 unmod = stack(unmod, spatinfo{1,2}(1,i).output{1}(:,6));
                end
            end
        end
        [h p] = kstest2(mod(~isnan(mod)),unmod(~isnan(unmod)));
        medmod = median(mod(~isnan(mod)));
        medunmod = median(unmod(~isnan(unmod)));
        sizes = [size(mod(~isnan(mod)),1) size(unmod(~isnan(unmod)),1) ];
        mod = mod(~isnan(mod));
        unmod = unmod(~isnan(unmod));
        figure
        hold on
        bar([1 2], [mean(mod) mean(unmod)], 0.8,'EdgeColor','none')
        errorbar2([1 2], [mean(mod) mean(unmod)],  [stderr(mod) stderr(unmod)] , 0.3, 'k')
        xlim([0.3 2.7])
        set(gca, 'fontsize', 24)
        set(gca, 'xtick', [1 2], 'xticklabel', {'mod', 'unmod'})
        ylabel('info (bits/spike)')
        rp = ranksum(mod,unmod);
        title({sprintf('All Animals n(%d) %s',length(animals),tit); sprintf('ranksum p = %s',rp); sprintf('K-Stest p = %s',p)})
        
        %finish the save fig stuff
        %     figfile = [figdir,area,'_',statename,'_Ripple',kind,'_CorrCoeffHist']
        %     if savefig1==1,
        %         print('-depsc2', figfile); print('-djpeg', figfile); saveas(gcf,figfile,'fig');
        %     end
        
    end
end
return % for now






%dont do ttest bc not normal dist
% [r_realrdm p_realrdm] = ttest(allp_shuf<0.05,allp_rdmshuf<0.05)
