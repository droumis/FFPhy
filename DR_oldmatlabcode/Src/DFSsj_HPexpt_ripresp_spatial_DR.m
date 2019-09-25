
warning('off','all');
clear;
runscript = 0;
savedata = 0; % save data option - only works if runscript is also on 
savefigs=1;
cyclemaps = 0; %1 to preview each map before it gets saved

%nadal
runnadal =0;
plotcomban = 1; %plot figs for all combined animals data. 0 for nadal or any individual HP animal
plotnovelyperiods = 1; % seperate d1-2, d3-5, d6-8...0 for nadal

plotstuff = 1;
ploteachan = 1;
plotoverdays = 1; 
plotmaps = 1;
plotPFCCA1maps = 1;
% plottypemaps =1;
plotmodmaps =1; 
plotunmodmaps =1;
modUnmod = [1 2]; % don't change. 
plotsizefield = 1;
plotcorrcoef =1;
plotrateXinfo = 1; %needs all plotmaps
collapseEPinday = 1; % don't change. 
whichep = 0; % don't change. use 0 for all epochs collapsed in day/cell

savedir = '/mnt/data25/sjadhav/HPExpt/ProcessedDataDR/';
savefilename = 'HPAllriplorentag2';
savefile = [savedir savefilename]; area = 'PFC'; %clrunmod = 'r'; clrmod = 'b'; % PFC  %%%% make this into a for loop to collect mod and unmod
figdir = '/mnt/data25/sjadhav/HPExpt/Figures_DR/';
HPdir = '/mnt/data25/sjadhav/HPExpt/';

Veqn = '>3';
minV=str2num(Veqn(end));
mintime = 3;
traj = [1:4] ;

% If runscript, run Datafilter and save data
if runscript == 1
    for i =  modUnmod;
        %Animal selection
        %-----------------------------------------------------
                                animals = {'HPa' 'HPb' 'HPc'};
%                         animals = {'nadal'};
        
        %Filter creation
        %-----------------------------------------------------
        
        % Epoch filter
        % -------------
        if runnadal == 1;
            dayfilter = '8:17';
        else
            dayfilter = '1:8';
        end

        if runnadal == 1;
            runepochfilter = 'isequal($type, ''run'')';
        else
            runepochfilter = 'isequal($environment, ''wtr1'') || isequal($environment, ''wtr2'')';
        end

        % %-----------
        if i < 2;
            placecellfilter = '(strcmp($area, ''PFC'') && strcmp($ripmodtag2, ''y'') && ($numspikes > 100))';   % Ripple mod tag2
%             placecellfilter = '(strcmp($area, ''PFC'') && strcmp($ripmodtag, ''y'') && ($numspikes > 100))';   % Ripple mod
%                         placecellfilter = '(strcmp($area, ''PFC'') && strcmp($thetamodtag, ''y'') && ($numspikes > 100))';   % theta mod
            
        else
            placecellfilter = '(strcmp($area, ''PFC'') && strcmp($ripmodtag2, ''n'') && ($numspikes > 100))'; % Ripple unmod tag2
%             placecellfilter = '(strcmp($area, ''PFC'') && strcmp($ripmodtag, ''n'') && ($numspikes > 100))'; % Ripple unmod
%                         placecellfilter = '(strcmp($area, ''PFC'') && strcmp($thetamodtag, ''n'') && ($numspikes > 100))'; % theta unmod
        end
        
        % Time filter -
        %%-----------

        riptetfilter = '(isequal($descrip, ''riptet''))';
        timefilter = { {'DFTFsj_getlinstate', '(($state ~= -1) & (abs($linearvel) >= 3))', 6}, {'DFTFsj_getriptimes','($nripples == 0)','tetfilter',riptetfilter,'minthresh',3} }; %DR added velocity filter.. trying to get ride of v high prococc in data

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
        %%%%spatf = setfilterfunction(spatf, 'getrawdata_DR', {'linpos', 'spikes'},6,  minV); % without any varargins. use for rawdata collect
        %      %%%%  spatfraw = setfilterfunction(spatf, 'getspatialinfo_DRraw', {'linpos', 'spikes'},6,  minV, 'peakrate', 3,'appendindex', 1, 'incltraj', traj); %I changed it to output all trajs' summed per cell, not combined across eps yet
        %         pof = setfilterfunction(spatf, 'DFAsj_twodoccupancy', {'spikes', 'linpos','pos'}, 'binsize', 1, 'std', 2);
        %use_____________________
        spatf = setfilterfunction(spatf, 'getspatialinfo_DR', {'linpos', 'spikes'},6,  minV, 'peakrate', 3,'appendindex', 1, 'incltraj', traj); %I changed it to output all trajs' summed per cell, not combined across eps yet
        psf = setfilterfunction(spatf, 'DFAsj_filtercalclinfields_tf',{'spikes', 'linpos'}, 'binsize', 2);
        pmf = setfilterfunction(spatf, 'DFAsj_openfieldrate_tf',{'spikes', 'linpos', 'pos'}, 'binsize', 1, 'std', 2);
        Expplaf = setfilterfunction(spatf, 'DFAsj_getplacefieldparams_DR', {'linpos', 'spikes'}, 'binsize', 2);
        
        % Run analysis-----------------------
        %         rawspatinfo{i} = runfilter(spatfraw);
        %         %use __________
        spatinfo{i} = runfilter(spatf); %spatial info
        pfm{i} = runfilter(pmf);  % Place Field Map
        pfs{i} = runfilter(psf);  % Place Field Stability.. trajectories
        plaf{i} = runfilter(Expplaf);  % proportion track fields
    end
    disp('Finished running filter script');
    %--------------------- Finished Filter Function Run -------------------
    
    if savedata == 1
        clear plotsizefield plotnovelyperiods runnadal figopt1 runscript plotrateXinfo    savedata plotstuff ploteachan plotcomban plotoverdays plotnoveltyperiods plotmaps plotmodmaps plotunmodmaps modunmod cyclemaps savefigs plottrajs
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
% PLOT!

% --------------------------------------------------------------------_____________________________________________________________________________________________________

% plot( each animal) do this first//
if plotstuff ==1
    close all
    %make directory..
    mkdir(figdir,[savefilename,'_Mod',num2str(plotmodmaps),'_Unmod',num2str(plotunmodmaps),'_Ep',num2str(whichep)])
    currfigdir = [savefilename,'_Mod',num2str(plotmodmaps),'_Unmod',num2str(plotunmodmaps),'_Ep',num2str(whichep),'/'];
    
    if ploteachan == 1
        
        for i = 1: length(animals);
            
            mkdir([figdir currfigdir],animals{i})
            anfigdir = [figdir currfigdir animals{i}];
            %_____________
            %combine the spat info of all trajs within each ep for each cell
            for mu = 1:2;
                clear dta tmpnum hu x y z c ia ic spatdata_alltraj gu
                dta = spatinfo{1,mu}(1,i).output{1};
                tmpnum = nan(length(dta(:,1)),1);
                for hu = 1:length(dta(:,1));
                    x = dta(hu,1:4);
                    y = sprintf('%d',x);
                    z = str2num(y);
                    tmpnum(hu) = z;
                end
                [c ia ic] = unique(tmpnum, 'stable');
                spatdata_alltraj = [];
                for gu = 1:length(ia);
                    spatdata_alltraj(gu,6) = nansum(dta(find(ic == gu),6));
                    spatdata_alltraj(gu,1:4) = dta(ia(gu),1:4); %take the index of first spt info number %leaving the trj column (5) empty
                end
                %_________________
                if collapseEPinday == 1;
                    clear  tmpnum x y z c ia ic spatdata_alltrajep gu
                    %combine the spat info of all trajs within each ep for each cell
                    tmpnum = nan(length(spatdata_alltraj(:,1)),1);
                    for hu = 1:length(spatdata_alltraj(:,1));
                        x = spatdata_alltraj(hu,[1 3 4]);
                        y = sprintf('%d',x);
                        z = str2num(y);
                        tmpnum(hu) = z;
                    end
                    [c ia ic] = unique(tmpnum, 'stable');
                    spatdata_alltrajep = [];
                    for gu = 1:length(ia);
                        spatdata_alltrajep(gu,6) = nanmean(spatdata_alltraj(find(ic == gu),6)); %TAKING THE MEAN OF EPS.. 
                        spatdata_alltrajep(gu,[1 3 4]) = spatdata_alltraj(ia(gu),[1 3 4]); %take the index of first spt info number %leaving the ep column (2) empty
                    end
                end
                
                
                %_________________
                if mu == 1;
                    spatmoddata = spatdata_alltrajep;
                    spatmoddata_alltraj = spatdata_alltraj; %use this for plotting maps with epochs separate but trajs combined
                else
                    spatunmoddata = spatdata_alltrajep;
                    spatunmoddata_alltraj = spatdata_alltraj;  %use this for plotting maps with epochs separate but trajs combined
                end
            end
            %___________________
           
            ndays = length(unique(spatinfo{1}(1,i).epochs{1}(:,1)));

            if whichep == 2 %use last ep... Don't use this anymore
                tit = 'Last W Epoch';
                mod = spatmoddata(find((spatmoddata(:,1)==1 & spatmoddata(:,2)>5) | (spatmoddata(:,1)>1 & spatmoddata(:,2)>3)),6); %the last run ep
                unmod = spatunmoddata(find((spatunmoddata(:,1)==1 & spatunmoddata(:,2)>5) | (spatunmoddata(:,1)>1 & spatunmoddata(:,2)>3)),6); %the last run ep
                
                %day chunking
                daymodall = [spatmoddata(find((spatmoddata(:,1)==1 & spatmoddata(:,2)>5) | (spatmoddata(:,1)>1 & spatmoddata(:,2)>3)),1)...
                    spatmoddata(find((spatmoddata(:,1)==1 & spatmoddata(:,2)>5) | (spatmoddata(:,1)>1 & spatmoddata(:,2)>3)),6)];
                dayunmodall = [spatunmoddata(find((spatunmoddata(:,1)==1 & spatunmoddata(:,2)>5) | (spatunmoddata(:,1)>1 & spatunmoddata(:,2)>3)),1)...
                    spatunmoddata(find((spatunmoddata(:,1)==1 & spatunmoddata(:,2)>5) | (spatunmoddata(:,1)>1 & spatunmoddata(:,2)>3)),6)];
                %over novelty periods:::: columns are mod and unmod. rows are
                %the means of the 3 periods. (d1-2; d3-5; d6-8)
                meannoveltyperiods = [[mean(daymodall(find(daymodall(:,1) == 1 | daymodall(:,1) == 2),2)) mean(dayunmodall(find(dayunmodall(:,1)==1 | dayunmodall(:,1)==2),2))];...
                    [mean(daymodall(find(daymodall(:,1) == 3 | daymodall(:,1) == 4 | daymodall(:,1) == 5 ),2)) mean(dayunmodall(find(dayunmodall(:,1)==3 | dayunmodall(:,1)==4 | dayunmodall(:,1)==5),2))];...
                    [mean(daymodall(find(daymodall(:,1) == 6 | daymodall(:,1) == 7 | daymodall(:,1) == 8 ),2)) mean(dayunmodall(find(dayunmodall(:,1)==6 | dayunmodall(:,1)==7 | dayunmodall(:,1)==8),2))]];
                %over days:::: columns are mod and unmod. rows are the
                %means for each day
                for j = 1:ndays;
                    daymod{j} = daymodall(find(daymodall(:,1)==j),2);
                    dayunmod{j} = dayunmodall(find(dayunmodall(:,1)==j),2);
                    if j == 1;
                        meanday = [mean(daymodall(find(daymodall(:,1)==j),2)) mean(dayunmodall(find(dayunmodall(:,1)==j),2))];
                    else
                        meanday = [meanday; [mean(daymodall(find(daymodall(:,1)==j),2)) mean(dayunmodall(find(dayunmodall(:,1)==j),2))]];
                    end
                end
                
                
                
            elseif whichep == 1 %use first ep. Don't use this anymore
                tit = 'First W Epoch';
                mod = spatmoddata(find((spatmoddata(:,1)==1 & spatmoddata(:,2)<5) | (spatmoddata(:,1)>1 & spatmoddata(:,2)<3)),6);
                unmod = spatunmoddata(find((spatunmoddata(:,1)==1 & spatunmoddata(:,2)>5) | (spatunmoddata(:,1)>1 & spatunmoddata(:,2)>3)),6);
                
                %day chunking
                daymodall = [spatmoddata(find((spatmoddata(:,1)==1 & spatmoddata(:,2)<5) | (spatmoddata(:,1)>1 & spatmoddata(:,2)<3)),1)...
                    spatmoddata(find((spatmoddata(:,1)==1 & spatmoddata(:,2)<5) | (spatmoddata(:,1)>1 & spatmoddata(:,2)<3)),6)];
                dayunmodall = [spatunmoddata(find((spatunmoddata(:,1)==1 & spatunmoddata(:,2)<5) | (spatunmoddata(:,1)>1 & spatunmoddata(:,2)<3)),1)...
                    spatunmoddata(find((spatunmoddata(:,1)==1 & spatunmoddata(:,2)<5) | (spatunmoddata(:,1)>1 & spatunmoddata(:,2)<3)),6)];
                %over novelty periods:::: columns are mod and unmod. rows are
                %the means of the 3 periods. (d1-2; d3-5; d6-8)
                meannoveltyperiods = [[mean(daymodall(find(daymodall(:,1) == 1 | daymodall(:,1) == 2),2)) mean(dayunmodall(find(dayunmodall(:,1)==1 | dayunmodall(:,1)==2),2))];...
                    [mean(daymodall(find(daymodall(:,1) == 3 | daymodall(:,1) == 4 | daymodall(:,1) == 5 ),2)) mean(dayunmodall(find(dayunmodall(:,1)==3 | dayunmodall(:,1)==4 | dayunmodall(:,1)==5),2))];...
                    [mean(daymodall(find(daymodall(:,1) == 6 | daymodall(:,1) == 7 | daymodall(:,1) == 8 ),2)) mean(dayunmodall(find(dayunmodall(:,1)==6 | dayunmodall(:,1)==7 | dayunmodall(:,1)==8),2))]];
                %over days:::: columns are mod and unmod. rows are the
                %means for each day
                for j = 1:ndays;
                    daymod{j} = daymodall(find(daymodall(:,1)==j),2); % i think I'll need these for errorbars
                    dayunmod{j} = dayunmodall(find(dayunmodall(:,1)==j),2);
                    if j == 1;
                        meanday = [mean(daymodall(find(daymodall(:,1)==j),2)) mean(dayunmodall(find(dayunmodall(:,1)==j),2))];
                    else
                        meanday = [meanday; [mean(daymodall(find(daymodall(:,1)==j),2)) mean(dayunmodall(find(dayunmodall(:,1)==j),2))]];
                    end
                end
                 
            else %use all eps.. USE THIS! (whichep = 0)
                tit = 'All W Epochs';
                mod = spatmoddata(:,6); %spatinfo{1,1}(1,i).output{1}(:,6);
                unmod = spatunmoddata(:,6); %spatinfo{1,2}(1,i).output{1}(:,6);
                %day chunking
                %over novelty periods:::: columns are mod and unmod. rows are
                %the means of the 3 periods. (d1-2; d3-5; d6-8)
                meannoveltyperiods = [[mean(spatmoddata(find(spatmoddata(:,1) == 1 | spatmoddata(:,1) == 2),6)) mean(spatunmoddata(find(spatunmoddata(:,1)==1 | spatunmoddata(:,1)==2),6))];...
                    [mean(spatmoddata(find(spatmoddata(:,1) == 3 | spatmoddata(:,1) == 4 | spatmoddata(:,1) == 5 ),6)) mean(spatunmoddata(find(spatunmoddata(:,1)==3 | spatunmoddata(:,1)==4 | spatunmoddata(:,1)==5),6))];...
                    [mean(spatmoddata(find(spatmoddata(:,1) == 6 | spatmoddata(:,1) == 7 | spatmoddata(:,1) == 8 ),6)) mean(spatunmoddata(find(spatunmoddata(:,1)==6 | spatunmoddata(:,1)==7 | spatunmoddata(:,1)==8),6))]];
                %over days:::: columns are mod and unmod. rows are the
                %means for each day
                
                %                 for j = 1:ndays;
                for j = [str2num(dayfilter)];
                    daymod{j} = spatmoddata(find(spatmoddata(:,1)==j),6);
                    dayunmod{j} = spatunmoddata(find(spatunmoddata(:,1)==j),6);
                    if runnadal == 1;
                        if j == 8; %nadal
                            meanday = [mean(spatmoddata(find(spatmoddata(:,1)==j),6)) mean(spatunmoddata(find(spatunmoddata(:,1)==j),6))];
                        else
                            meanday = [meanday; [mean(spatmoddata(find(spatmoddata(:,1)==j),6)) mean(spatunmoddata(find(spatunmoddata(:,1)==j),6))]];
                        end
                    else
                        if j == 1; %HP
                            meanday = [mean(spatmoddata(find(spatmoddata(:,1)==j),6)) mean(spatunmoddata(find(spatunmoddata(:,1)==j),6))];
                        else
                            meanday = [meanday; [mean(spatmoddata(find(spatmoddata(:,1)==j),6)) mean(spatunmoddata(find(spatunmoddata(:,1)==j),6))]];
                        end
                    end
                    
                end
                
            end
            % plotting each an
            edges = linspace(min([mod; unmod]),max([mod; unmod]),20); %%
            modcnt = histc(mod,edges);
            unmodcnt = histc(unmod,edges);
            figure(i)
            hold on
            subplot(2,2,1)
            bar(edges,modcnt/sum(modcnt),'b');
            alpha(0.65);
%             set(gcf,'FaceAlpha',.5);
%             plot(edges, modcnt/sum(modcnt), 'Color',[0 0 .8], 'linewidth', 3);
            hold on
%             plot(edges, unmodcnt/sum(unmodcnt), 'Color',[0.8 0 0], 'linewidth', 3);
            bar(edges,unmodcnt/sum(unmodcnt),'r');
            alpha(0.65);
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
            subplot(2,2,2);
            hold on
            bar(1, mean(mod), 'b', 'EdgeColor', 'none')
            bar(2, mean(unmod),'r', 'EdgeColor', 'none')
            errorbar2([1 2], [mean(mod) mean(unmod)],  [stderr(mod) stderr(unmod)] , 0.3, 'k')
            xlim([0.3 2.7])
            set(gca, 'fontsize', 24)
            set(gca, 'xtick', [1 2], 'xticklabel', {'mod', 'unmod'})
            ylabel('info (bits/spike)')
            rp = ranksum(mod,unmod);
            title({['Spat Info' animals{1,i}]; sprintf('ranksum p = %s',rp); sprintf('K-Stest p = %s',p)})
            
            %plot info over days
            if plotoverdays == 1;
                subplot(2,2,3);
                hold on
                bar([str2num(dayfilter)], meanday, 'EdgeColor','none') %hpa animals nadal
                %                 errorbar2(1:ndays, meanday,  [stderr(mod) stderr(unmod)] , 0.3, 'k') %add error bars later if something interesting
                
                if runnadal == 1;
                    xlim([7.5 17.5]) %nadal
                else
                    xlim([0 8.5]) %HP
                end
                %
                set(gca, 'fontsize', 24)
                
                set(gca, 'xtick',[str2num(dayfilter)]) %nadal
                %                 set(gca, 'xtick', 1:ndays)
                ylabel('info (bits/spike)')
                xlabel('days')
                title('spat info over days');
            end
            
            %plot info novel1(d1,2) vs familiar(d3-5) vs novel2(d6-..)
            if plotnovelyperiods ==1;
                subplot(2,2,4);
                hold on
                bar(1:length(meannoveltyperiods), meannoveltyperiods, 'EdgeColor','none')
                %                 errorbar2(1:ndays, meanday,  [stderr(mod) stderr(unmod)] , 0.3, 'k') %add error bars later if something interesting
                xlim([0 3.5])
                set(gca, 'fontsize', 15)
                set(gca,  'xtick', [1:3], 'xticklabel', {'1:2', '3:5', '6:8'})
                ylabel('info (bits/spike)')
                xlabel('day groups')
                title('spat info over novelty periods');
            end

            %plot excit vs inhib modulated using typetag in cellinfo
%             eval(load([HPdir,sprintf('%s_direct/%scellinfo.mat',animals{1,i},animals{1,i})]));

            
            
            
            
            
            
            %save the info figure and close it
            figfile = [anfigdir,'/','InfoFig' ];
            if savefigs==1
                print('-djpeg', figfile);
                %saveas(gcf,figfile,'fig');
                %print('-depsc2', figfile);
            end
            if~cyclemaps ==0
                keyboard
            end
            close
            
            
            %____________________________________________________________________________________________
            %plotting size of field

            if plotsizefield ==1;
                figure
                % total frac under 3
                for mu = 1:2;
                    clear  tmpnum x y z c ia ic spatdata_alltrajep gu indexdata meanday
                    currindexmod=[]; cntallcells = 0; total_fracunder3 = []; pu =0;
                    nidxs_an = length(plaf{mu}(i).output{1});
                    for wu=1:nidxs_an
                        indexdata(wu,:) = [plaf{mu}(i).output{1}(wu).index plaf{mu}(i).output{1}(wu).total_fracunder3];
                        currindexmod(wu,:) = plaf{mu}(i).output{1}(wu).index;
                    end
                    tmpnum = nan(length(currindexmod(:,1)),1);
                    for hu = 1:length(currindexmod(:,1));
                        x = currindexmod(hu,[1 3 4]);
                        y = sprintf('%d',x);
                        z = str2num(y);
                        tmpnum(hu) = z;
                    end
                    [c ia ic] = unique(tmpnum, 'stable');
                    sizedata_allep = [];
                    for gu = 1:length(ia);
                        sizedata_allep(gu,5) = nanmean(indexdata((find(ic == gu)),5));
                        sizedata_allep(gu,[1 3 4]) = indexdata((ia(gu)),[1 3 4]);
                    end
                    if mu == 1;
                        sizefrac3moddata = sizedata_allep;
                    else
                        sizefrac3unmoddata = sizedata_allep;
                    end
                end
                
                subplot(2,2,1);
                hold on
                bar(1,mean(sizefrac3moddata(:,5)),'b');
                bar(2,mean(sizefrac3unmoddata(:,5)),'r');
                errorbar(1,mean(sizefrac3moddata(:,5)),stderr(sizefrac3moddata(:,5)), 0.3, 'k','LineWidth',3);
                errorbar(2,mean(sizefrac3unmoddata(:,5)),stderr(sizefrac3unmoddata(:,5)), 0.3, 'k','LineWidth',3);
                xlim([0.3 2.7])
                set(gca, 'fontsize', 24)
                [h_frac3,p_frac3] = ttest2(sizefrac3moddata(:,5),sizefrac3unmoddata(:,5)); %h=0, p=0.77
                title({['Proportion track active above 3Hz'];[sprintf('p = %d',p_frac3)]});
                ylabel('Proportion track active');
                set(gca,'XTick',[1 2],'XTickLabel',{'Mod';'Unmod'});
                
                %plot size over days
                for j = [str2num(dayfilter)];
                    daymod{j} = sizefrac3moddata(find(sizefrac3moddata(:,1)==j),5);
                    dayunmod{j} = sizefrac3unmoddata(find(sizefrac3unmoddata(:,1)==j),5);
                    if runnadal == 1;
                        if j == 8; %nadal
                            meanday = [mean(sizefrac3moddata(find(sizefrac3moddata(:,1)==j),5)) mean(sizefrac3unmoddata(find(sizefrac3unmoddata(:,1)==j),5))];
                        else
                            meanday = [meanday; [mean(sizefrac3moddata(find(sizefrac3moddata(:,1)==j),5)) mean(sizefrac3unmoddata(find(sizefrac3unmoddata(:,1)==j),5))]];
                        end
                    else
                        if j == 1; %HP
                            meanday = [mean(sizefrac3moddata(find(sizefrac3moddata(:,1)==j),5)) mean(sizefrac3unmoddata(find(sizefrac3unmoddata(:,1)==j),5))];
                        else
                            meanday = [meanday; [mean(sizefrac3moddata(find(sizefrac3moddata(:,1)==j),5)) mean(sizefrac3unmoddata(find(sizefrac3unmoddata(:,1)==j),5))]];
                        end
                    end
                    
                end
                
                subplot(2,2,3);
                hold on
                bar([str2num(dayfilter)], meanday, 'EdgeColor','none') %hpa animals nadal
                %                 errorbar2(1:ndays, meanday,  [stderr(mod) stderr(unmod)] , 0.3, 'k') %add error bars later if something interesting
                if runnadal == 1;
                    xlim([7.5 17.5]) %nadal
                else
                    xlim([0 8.5]) %HP
                end
                %
                set(gca, 'fontsize', 24)
                
                set(gca, 'xtick',[str2num(dayfilter)]) %nadal
                %                 set(gca, 'xtick', 1:ndays)
                ylabel('Proportion field size')
                xlabel('days')
                title('Proportion field size >3Hz over days');

                
                %field size using mean threshold per
                %cell_____________________________________________
                
                for mu = 1:2;
                    clear  tmpnum x y z c ia ic spatdata_alltrajep gu indexdata meanday
                    currindexmod=[]; cntallcells = 0; total_fracundermean = []; pu =0;
                    nidxs_an = length(plaf{mu}(i).output{1});
                    for wu=1:nidxs_an
                        indexdata(wu,:) = [plaf{mu}(i).output{1}(wu).index plaf{mu}(i).output{1}(wu).total_fracundermean];
                        currindexmod(wu,:) = plaf{mu}(i).output{1}(wu).index;
                    end
                    tmpnum = nan(length(currindexmod(:,1)),1);
                    for hu = 1:length(currindexmod(:,1));
                        x = currindexmod(hu,[1 3 4]);
                        y = sprintf('%d',x);
                        z = str2num(y);
                        tmpnum(hu) = z;
                    end
                    [c ia ic] = unique(tmpnum, 'stable');
                    sizedata_allep = [];
                    for gu = 1:length(ia);
                        sizedata_allep(gu,5) = nanmean(indexdata((find(ic == gu)),5));
                        sizedata_allep(gu,[1 3 4]) = indexdata((ia(gu)),[1 3 4]);
                    end
                    if mu == 1;
                        sizefracmeanmoddata = sizedata_allep;
                    else
                        sizefracmeanunmoddata = sizedata_allep;
                    end
                end
                %over days
                subplot(2,2,2);
                hold on
                bar(1,mean(sizefracmeanmoddata(:,5)),'b');
                bar(2,mean(sizefracmeanunmoddata(:,5)),'r');
                errorbar(1,mean(sizefracmeanmoddata(:,5)),stderr(sizefracmeanmoddata(:,5)), 0.3, 'k','LineWidth',3);
                errorbar(2,mean(sizefracmeanunmoddata(:,5)),stderr(sizefracmeanunmoddata(:,5)), 0.3, 'k','LineWidth',3);
                xlim([0.3 2.7])
                set(gca, 'fontsize', 24)
                [h_frac3,p_frac3] = ttest2(sizefracmeanmoddata(:,5),sizefracmeanunmoddata(:,5)); %h=0, p=0.77
                title({['Proportion track active above mean'];[sprintf('p = %d',p_frac3)]});
                ylabel('Proportion track active');
                set(gca,'XTick',[1 2],'XTickLabel',{'Mod';'Unmod'});
                
                %plot size over days
                for j = [str2num(dayfilter)];
                    daymod{j} = sizefracmeanmoddata(find(sizefracmeanmoddata(:,1)==j),5);
                    dayunmod{j} = sizefracmeanunmoddata(find(sizefracmeanunmoddata(:,1)==j),5);
                    if runnadal == 1;
                        if j == 8; %nadal
                            meanday = [mean(sizefracmeanmoddata(find(sizefracmeanmoddata(:,1)==j),5)) mean(sizefracmeanunmoddata(find(sizefracmeanunmoddata(:,1)==j),5))];
                        else
                            meanday = [meanday; [mean(sizefracmeanmoddata(find(sizefracmeanmoddata(:,1)==j),5)) mean(sizefracmeanunmoddata(find(sizefracmeanunmoddata(:,1)==j),5))]];
                        end
                    else
                        if j == 1; %HP
                            meanday = [mean(sizefracmeanmoddata(find(sizefracmeanmoddata(:,1)==j),5)) mean(sizefracmeanunmoddata(find(sizefracmeanunmoddata(:,1)==j),5))];
                        else
                            meanday = [meanday; [mean(sizefracmeanmoddata(find(sizefracmeanmoddata(:,1)==j),5)) mean(sizefracmeanunmoddata(find(sizefracmeanunmoddata(:,1)==j),5))]];
                        end
                    end
                    
                end
                
                subplot(2,2,4);
                hold on
                bar([str2num(dayfilter)], meanday, 'EdgeColor','none') %hpa animals nadal
                %                 errorbar2(1:ndays, meanday,  [stderr(mod) stderr(unmod)] , 0.3, 'k') %add error bars later if something interesting
                if runnadal == 1;
                    xlim([7.5 17.5]) %nadal
                else
                    xlim([0 8.5]) %HP
                end
                %
                set(gca, 'fontsize', 24)
                
                set(gca, 'xtick',[str2num(dayfilter)]) %nadal
                %                 set(gca, 'xtick', 1:ndays)
                ylabel('Proportion field size')
                xlabel('days')
                title('Proportion field size >mean/cell over days');

                
                figfile = [anfigdir,'/','FieldSizeFig' ];
                if savefigs==1
                    print('-djpeg', figfile);
                    %saveas(gcf,figfile,'fig');
                    %print('-depsc2', figfile);
                end
                if~cyclemaps ==0
                    keyboard
                end
                close
            end
            
            %______________________________________________________________________________________________________________________________________
            %corr coef
            if plotcorrcoef ==1;
                
                for k = 1:length(pfs{1}(i).output{1}(1,:));
                    currind = [];
                    if length(pfs{1}(i).output{1}(1,k).trajdata) ==4;  % if all trajs exist
                        
                        %mod gather trajs
                        currind = pfs{1}(i).output{1}(1,k).index;
                        
                        outleftmod = pfs{1}(i).output{1}(1,k).trajdata{1,1}(:,5);
                        outrightmod= pfs{1}(i).output{1}(1,k).trajdata{1,3}(:,5);
                        
                        inleftmod = pfs{1}(i).output{1}(1,k).trajdata{1,2}(:,5);
                        inrightmod = pfs{1}(i).output{1}(1,k).trajdata{1,4}(:,5);
                        
                        %compute corrcoef
                        outmodmintrunc = min(length(outleftmod), length(outrightmod));
                        [outmodR outmodP] = corrcoef(outrightmod(1:outmodmintrunc),outleftmod(1:outmodmintrunc),'rows','pairwise'); % rows-pairwise ~ compute R(i,j) using rows with no NaN values in either column i or j.
                        AllR{1}(k,5)  = outmodR(2,1);
                        AllR{1}(k,[1 3 4]) = currind(1,[1 3 4]);
                        outmodAllP(k,5) = outmodP(2,1);
                        outmodAllP(k,[1 2 3 4]) = currind(1,[1 2 3 4]);
                        
                        inmodmintrunc = min(length(inleftmod), length(inrightmod));
                        [inmodR inmodP] = corrcoef(inrightmod(1:inmodmintrunc),inleftmod(1:inmodmintrunc),'rows','pairwise'); % rows-pairwise ~ compute R(i,j) using rows with no NaN values in either column i or j.
                        AllR{2}(k,5)  = inmodR(2,1);
                        AllR{2}(k,[1 3 4]) = currind(1,[1 3 4]);
                        inmodAllP(k,5) = inmodP(2,1);
                        inmodAllP(k,[1 2 3 4]) = currind(1,[1 2 3 4]);
                        
                        
                    else
                        pfs{1}(i).output{1}(1,k).index %print any cells skipped
                        
                    end
                end
                
                %UNMOD________________
                for k = 1:length(pfs{2}(i).output{1}(1,:));
                    currind = [];
                    if length(pfs{2}(i).output{1}(1,k).trajdata) ==4;  % if all trajs exist
                        
                        %gather trajs
                        currind = pfs{2}(i).output{1}(1,k).index;
                        
                        outleftunmod = pfs{2}(i).output{1}(1,k).trajdata{1,1}(:,5);
                        outrightunmod= pfs{2}(i).output{1}(1,k).trajdata{1,3}(:,5);
                        
                        inleftunmod = pfs{2}(i).output{1}(1,k).trajdata{1,2}(:,5);
                        inrightunmod = pfs{2}(i).output{1}(1,k).trajdata{1,4}(:,5);
                        
                        %compute corrcoef
                        outunmodmintrunc = min(length(outleftunmod), length(outrightunmod));
                        [outunmodR outunmodP] = corrcoef(outrightunmod(1:outunmodmintrunc),outleftunmod(1:outunmodmintrunc),'rows','pairwise'); % rows-pairwise ~ compute R(i,j) using rows with no NaN values in either column i or j.
                        AllR{3}(k,5)  = outunmodR(2,1);
                        AllR{3}(k,[1 3 4]) = currind(1,[1 3 4]);
                        outunmodAllP(k,5) = outunmodP(2,1);
                        outunmodAllP(k,[1 2 3 4]) = currind(1,[1 2 3 4]);
                        
                        inunmodmintrunc = min(length(inleftunmod), length(inrightunmod));
                        [inunmodR inunmodP] = corrcoef(inrightunmod(1:inunmodmintrunc),inleftunmod(1:inunmodmintrunc),'rows','pairwise'); % rows-pairwise ~ compute R(i,j) using rows with no NaN values in either column i or j.
                        AllR{4}(k,5)  = inunmodR(2,1);
                        AllR{4}(k,[1 3 4]) = currind(1,[1 3 4]);
                        inunmodAllP(k,5) = inunmodP(2,1);
                        inunmodAllP(k,[1 2 3 4]) = currind(1,[1 2 3 4]);
                        
                    else
                        pfs{2}(i).output{1}(1,k).index %print any cells skipped
                        
                    end
                end
                
                %day chunking
                
                for mu = 1:length(AllR);
                    clear  tmpnum x y z c ia ic gu
                    currindexmod=[]; cntallcells = 0; total_fracunder3 = []; pu =0;
                    nidxs_an = length(AllR{mu});
                    tmpnum = nan(length(AllR{mu}(:,1)),1);
                    for hu = 1:length(AllR{mu}(:,1));
                        x = AllR{mu}(hu,[1 3 4]);
                        y = sprintf('%d',x);
                        z = str2num(y);
                        tmpnum(hu) = z;
                    end
                    [c ia ic] = unique(tmpnum, 'stable');
                    for gu = 1:length(ia);
                        corrdata_allep{mu}(gu,5) = nanmean(AllR{mu}((find(ic == gu)),5));
                        corrdata_allep{mu}(gu,[1 3 4]) = AllR{mu}((ia(gu)),[1 3 4]);
                    end
                end
                
                
                corrdata_allep{1} = corrdata_allep{1}((~isnan(corrdata_allep{1}(:,5))),:);
                corrdata_allep{3} = corrdata_allep{3}((~isnan(corrdata_allep{3}(:,5))),:);
                corrdata_allep{2} = corrdata_allep{2}((~isnan(corrdata_allep{2}(:,5))),:);
                corrdata_allep{4} = corrdata_allep{4}((~isnan(corrdata_allep{4}(:,5))),:);
                
                figure
                subplot(2,2,1);
                hold on
                
                bar(1, mean(corrdata_allep{1}(:,5)), 'b', 'EdgeColor', 'none')
                bar(2, mean(corrdata_allep{3}(:,5)),'r', 'EdgeColor', 'none')
                errorbar2([1 2], [mean(corrdata_allep{1}(:,5)) mean(corrdata_allep{3}(:,5))],  [stderr(corrdata_allep{1}(:,5)) stderr(corrdata_allep{3}(:,5))] , 0.3, 'k')
                xlim([0.3 2.7])
                set(gca, 'fontsize', 24)
                set(gca, 'xtick', [1 2], 'xticklabel', {'mod', 'unmod'})
                ylabel('corrcoef R')
                [h p] = kstest2(corrdata_allep{1}(~isnan(corrdata_allep{1}(:,5)),5),corrdata_allep{3}(~isnan(corrdata_allep{3}(:,5)),5));
                rp = ranksum(corrdata_allep{1}(~isnan(corrdata_allep{1}(:,5)),5),corrdata_allep{3}(~isnan(corrdata_allep{3}(:,5)),5));
                title({['CorrCoef Outbound' animals{1,i}]; sprintf('ranksum p = %s',rp); sprintf('K-Stest p = %s',p); [sprintf('n sig mod= %d unmod= %d',length(find(outmodAllP(:,5) < 0.05)), length(find(outunmodAllP(:,5) < 0.05)))]})
                
                subplot(2,2,2);
                hold on
                bar(1, mean(corrdata_allep{2}(:,5)), 'b', 'EdgeColor', 'none')
                bar(2, mean(corrdata_allep{4}(:,5)),'r', 'EdgeColor', 'none')
                errorbar2([1 2], [mean(corrdata_allep{2}(:,5)) mean(corrdata_allep{4}(:,5))],  [stderr(corrdata_allep{2}(:,5)) stderr(corrdata_allep{4}(:,5))] , 0.3, 'k')
                xlim([0.3 2.7])
                set(gca, 'fontsize', 24)
                set(gca, 'xtick', [1 2], 'xticklabel', {'mod', 'unmod'})
                ylabel('corrcoef R')
                [h p] = kstest2(corrdata_allep{2}(~isnan(corrdata_allep{2}(:,5)),5),corrdata_allep{4}(~isnan(corrdata_allep{4}(:,5)),5));
                rp = ranksum(corrdata_allep{2}(~isnan(corrdata_allep{2}(:,5)),5),corrdata_allep{4}(~isnan(corrdata_allep{4}(:,5)),5));
                title({['CorrCoef Inbound__' animals{1,i}]; sprintf('ranksum p = %s',rp); sprintf('K-Stest p = %s',p); [sprintf('n sig mod= %d unmod= %d',length(find(inmodAllP(:,5) < 0.05)), length(find(inunmodAllP(:,5) < 0.05)))]})
            end
            
            
            %gather indices of all significant corrcoef
            SigCorrCoefList{1} = {outmodAllP(find(outmodAllP(:,5) < 0.05),:)};
            SigCorrCoefList{2} = {inmodAllP(find(inmodAllP(:,5) < 0.05),:)};
            SigCorrCoefList{3} = {outunmodAllP(find(outunmodAllP(:,5) < 0.05),:)};
            SigCorrCoefList{4} = {inunmodAllP(find(inunmodAllP(:,5) < 0.05),:)};
            
            corrfile = [anfigdir,'/','SigCorrCoefList' ];
            if savefigs==1
                save(corrfile, 'SigCorrCoefList');
                %                 save(corrfile, Allsigncorrcoef);
            end
            
            
            %plot corr coef over days
            outmoddayR = corrdata_allep{1};
            outunmoddayR = corrdata_allep{3};
            inmoddayR = corrdata_allep{2};
            inunmoddayR = corrdata_allep{4};
            
            for j = [str2num(dayfilter)];
                if runnadal == 1;
                    if j == 8; %nadal
                        outmeanday = [mean(outmoddayR(find(outmoddayR(:,1)==j),5)) mean(outunmoddayR(find(outunmoddayR(:,1)==j),5))];
                        inmeanday = [mean(inmoddayR(find(inmoddayR(:,1)==j),5)) mean(inunmoddayR(find(inunmoddayR(:,1)==j),5))];
                    else
                        outmeanday = [outmeanday; [mean(outmoddayR(find(outmoddayR(:,1)==j),5)) mean(outunmoddayR(find(outunmoddayR(:,1)==j),5))]];
                        inmeanday = [inmeanday; [mean(inmoddayR(find(inmoddayR(:,1)==j),5)) mean(inunmoddayR(find(inunmoddayR(:,1)==j),5))]];
                    end
                else
                    if j == 1; %HP
                        outmeanday = [mean(outmoddayR(find(outmoddayR(:,1)==j),5)) mean(outunmoddayR(find(outunmoddayR(:,1)==j),5))];
                        inmeanday = [mean(inmoddayR(find(inmoddayR(:,1)==j),5)) mean(inunmoddayR(find(inunmoddayR(:,1)==j),5))];
                    else
                        outmeanday = [outmeanday; [mean(outmoddayR(find(outmoddayR(:,1)==j),5)) mean(outunmoddayR(find(outunmoddayR(:,1)==j),5))]];
                        inmeanday = [inmeanday; [mean(inmoddayR(find(inmoddayR(:,1)==j),5)) mean(inunmoddayR(find(inunmoddayR(:,1)==j),5))]];
                    end
                end
            end
            
            subplot(2,2,3);
            hold on
            bar([str2num(dayfilter)], outmeanday, 'EdgeColor','none') %hpa animals nadal
            %                 errorbar2(1:ndays, meanday,  [stderr(mod) stderr(unmod)] , 0.3, 'k') %add error bars later if something interesting
            if runnadal == 1;
                xlim([7.5 17.5]) %nadal
            else
                xlim([0 8.5]) %HP
            end
            %
            set(gca, 'fontsize', 24)
            set(gca, 'xtick',[str2num(dayfilter)]) %nadal
            %                 set(gca, 'xtick', 1:ndays)
            ylabel('corrcoef R')
            xlabel('days')
            title('Corrcoef Outbound over days');
            
            subplot(2,2,4);
            hold on
            bar([str2num(dayfilter)], inmeanday, 'EdgeColor','none') %hpa animals nadal
            %                 errorbar2(1:ndays, meanday,  [stderr(mod) stderr(unmod)] , 0.3, 'k') %add error bars later if something interesting
            if runnadal == 1;
                xlim([7.5 17.5]) %nadal
            else
                xlim([0 8.5]) %HP
            end
            %
            set(gca, 'fontsize', 24)
            set(gca, 'xtick',[str2num(dayfilter)]) %nadal
            %                 set(gca, 'xtick', 1:ndays)
            ylabel('corrcoef R')
            xlabel('days')
            title('Corrcoef Inbound over days');
            
            
            
            figfile = [anfigdir,'/','CorrCoefFig' ];
            if savefigs==1
                print('-djpeg', figfile);
                %saveas(gcf,figfile,'fig');
                %print('-depsc2', figfile);
            end
            if~cyclemaps ==0
                keyboard
            end
            close
            
            
            %___________________________________________________________________________________________________________________________________________
            if plotmaps == 1;
                
                
%                 if plottypemaps ==1;
%                     load([HPdir,sprintf('%s_direct/%scellinfo.mat',animals{1,i},animals{1,i})]);
%                     
%                  trajdata = [];
%                     figcnt = 0; totalplots = 0; totalmapplots=0;%Count total plots across figures
%                     clr = {'b','r','g','m','c','y','k','r'};
%                     figcnt = 0; totalplots = 0; maxrates = [];%Count total plots across figures
%                     for k=1:length(pfm{1}(i).output{1}); %modulated
%                         typeind = pfm{1}(i).output{1}(k).index;
%                         typecell = cellinfo{typeind(1)}{typeind(2)}{typeind(3)}{typeind(4)}.typetag;
%                         %plot excitated cells
%                         %______________________
%                         
%                         
%                         alltrajdata{i}{k}=pfs{1}(i).output{1}(k);
%                         if rem(totalmapplots+1,20)==0; %  subplots finished
%                             figcnt=figcnt+1;
%                         end
%                         figure(figcnt+1);  hold on;
%                         subplot(4,5,rem(totalplots,20)+1); hold on;
%                         trajdata = alltrajdata{i}{k}.trajdata;
%                         for o=1:length(trajdata),
%                             plot(trajdata{o}(:,5),[clr{o} '.-'],'Linewidth',2);
% %                             maxrate(o) = max(trajdata{o}(5:end-5,5));
%                         end
%                         if rem(totalplots+1,20)==1;
%                             legend('OL','IL','OR','IR');
%                             xlabel ('Position along linear trajectory (cm)','FontSize',14,'Fontweight','bold');
%                             ylabel ('Firing Rate (Hz)','FontSize',14,'Fontweight','bold');
%                             titstring = [animals{i}, '________ModEXCITED______', tit, '_____cell group:', num2str(figcnt)];
%                         end
%                         set(gca,'xtick',[])
%                         set(gca,'xticklabel',[])
%                         set(gca,'ytick',[])
%                         set(gca,'yticklabel',[])
%                         
%                         %________________________ FIRING RATE MAPS
%                         allmapdata{i}{k} = pfm{1}(i).output{1}(k); %just mod for now
%                         totalmapplots = totalplots+5;
%                         subplot(4,5,rem(totalmapplots,20)+1); hold on;
%                         imagesc(allmapdata{i}{k}.smoothedspikerate')
%                         %                         set(gca,'YLim',[0 110]);
%                         %                         set(gca,'XLim',[0 117]);
%                         
%                         tmpspatind = spatmoddata_alltraj(:,[1:4]);
%                         if ~rowfind(allmapdata{i}{k}.index,tmpspatind) == 0;
%                             currmapspatinfo = spatmoddata_alltraj(rowfind(allmapdata{i}{k}.index,tmpspatind),6);
%                             if currmapspatinfo <1;
%                                 color = 'b';
%                             elseif currmapspatinfo <3;
%                                 color = [0 .5 0];
%                             elseif currmapspatinfo <10;
%                                 color = 'r';
%                             elseif currmapspatinfo >=10;
%                                 color = 'k';
%                                 
%                             end
%                             %                             title([num2str(currmapspatinfo)], 'FontSize',24,'Fontweight','bold', 'Color', color);
%                             title({['SI=',num2str(currmapspatinfo)]; ['maxF=',num2str(max(max(allmapdata{i}{k}.smoothedspikerate)))]; [sprintf('D%d E%d T%d C%d',pfm{1}(i).output{1}(k).index(1),pfm{1}(i).output{1}(k).index(2),pfm{1}(i).output{1}(k).index(3),pfm{1}(i).output{1}(k).index(4))]}, 'FontSize',24,'Fontweight','bold', 'Color', color);
%                             maxratesmod(k) = max(max(allmapdata{i}{k}.smoothedspikerate));
%                             spatinfoallmod(k) = currmapspatinfo;
%                         else
%                             title(['~exist'], 'FontSize',24,'Fontweight','bold');
%                         end
%                         set(gca,'xtick',[])
%                         set(gca,'xticklabel',[])
%                         set(gca,'ytick',[])
%                         set(gca,'yticklabel',[])
%                         % Update plot number
%                         totalplots=totalplots+1;
%                         if (rem(totalplots,20) >4 & rem(totalplots,20) <10) | (rem(totalplots,20) >14);
%                             totalplots = totalplots+5; %skip the firing rate maps lines
%                         end
%                         if rem(totalmapplots+1,20)==0
%                             
%                             %save the full maps figures
%                             figfile = [anfigdir,'/','EXCMaps',num2str(figcnt)];
%                             [ax,hU]=suplabel(titstring,'t', [.08 .08 .84 .84]);
%                             set(hU,'FontSize',200,'Fontweight','bold')  %font controls aren't working.. i think it's a bug with ubuntu matlab2013a...
%                             if savefigs==1
%                                 print('-djpeg', figfile);
%                                 %saveas(gcf,figfile,'fig');
%                                 %print('-depsc2', figfile);
%                             end
%                             
%                             if ~cyclemaps == 0
%                                 keyboard;
%                             end
%                             close
%                         end
%                     end
%                     
%                     %save the last maps figure as it's most likely
%                     %not full
%                     figfile = [anfigdir,'/','EXCMaps',num2str(figcnt)];
%                     [ax,hU]=suplabel(titstring,'t', [.08 .08 .84 .84]);
%                     set(hU,'FontSize',200,'Fontweight','bold')  %font controls aren't working.. i think it's a bug with ubuntu matlab2013a...
%                     if savefigs==1
%                         print('-djpeg', figfile);
%                         %saveas(gcf,figfile,'fig');
%                         %print('-depsc2', figfile);
%                     end
%                     if~cyclemaps ==0
%                         keyboard
%                     end
%                     close
%                 end
%                 
                %__________________________________________________________________________________________________________________

            
            
                            %_____________________________________________________________________________________________________________________________________________________________
                
                if plotmodmaps ==1;
                    load([HPdir,sprintf('%s_direct/%scellinfo.mat',animals{1,i},animals{1,i})]);
                    trajdata = [];
                    figcnt = 0; totalplots = 0; totalmapplots=0;%Count total plots across figures
                    clr = {'b','r','g','m','c','y','k','r'};
                    figcnt = 0; totalplots = 0; maxrates = [];%Count total plots across figures
                    for k=1:length(pfm{1}(i).output{1}); %modulated
                        %______________________
                        
                        alltrajdata{i}{k}=pfs{1}(i).output{1}(k);
                        if rem(totalmapplots+1,20)==0; %  subplots finished
                            figcnt=figcnt+1;
                        end
                        figure(figcnt+1);  hold on;
                        subplot(4,5,rem(totalplots,20)+1); hold on;
                        trajdata = alltrajdata{i}{k}.trajdata;
                        for o=1:length(trajdata),
                            plot(trajdata{o}(:,5),[clr{o} '.-'],'Linewidth',2);
%                             maxrate(o) = max(trajdata{o}(5:end-5,5));
                        end
                        typeind = pfm{1}(i).output{1}(k).index;
                        typecell = cellinfo{typeind(1)}{typeind(2)}{typeind(3)}{typeind(4)}.typetag;
                        
                        if typecell == 'exc';
                            typclr = [0 .7 .93];
                        elseif typecell == 'inh';
                            typclr = [.87 1 .65];
                        else
                            typclr = [.8 .8 .8];
                        end
                        title(typecell,'Color',typclr)
                        set(gca,'Color',typclr);
                        set(gcf, 'InvertHardCopy', 'off');
                        if rem(totalplots+1,20)==1;
                            
                            xlabel ('Position along linear trajectory (cm)','FontSize',14,'Fontweight','bold');
                            ylabel ('Firing Rate (Hz)','FontSize',14,'Fontweight','bold');
                            %                             legend('OL','IL','OR','IR');
                            legfix = legend('OL','IL','OR','IR');
                            set(legfix, 'Color', 'none');
                            hText = findobj(legfix, 'type', 'text');
                            set(hText(4),'color',  [0 0 1]);
                            set(hText(3),'color',  [1 0 0]);
                            set(hText(2),'color',  [0 1 0]);
                            set(hText(1),'color',  [1 0 1]);
                            linesInPlot = findobj(legfix, 'type', 'line');
                            set(linesInPlot(1:8),'XData',[0.5 0.7]);

                            
                            titstring = [animals{i}, '________Mod______', tit, '_____cell group:', num2str(figcnt)];
                        end
                        set(gca,'xtick',[])
                        set(gca,'xticklabel',[])
                        set(gca,'ytick',[])
                        set(gca,'yticklabel',[])
                        
                        %________________________ FIRING RATE MAPS
                        allmapdata{i}{k} = pfm{1}(i).output{1}(k); %just mod for now
                        totalmapplots = totalplots+5;
                        subplot(4,5,rem(totalmapplots,20)+1); hold on;
                        imagesc(allmapdata{i}{k}.smoothedspikerate')
                        %                         set(gca,'YLim',[0 110]);
                        %                         set(gca,'XLim',[0 117]);
                        
                        titstring = [animals{i}, '________Mod______', tit, '_____cell group:', num2str(figcnt)];
                        
                        
                        tmpspatind = spatmoddata_alltraj(:,[1:4]);
                        if ~rowfind(allmapdata{i}{k}.index,tmpspatind) == 0;
                            currmapspatinfo = spatmoddata_alltraj(rowfind(allmapdata{i}{k}.index,tmpspatind),6);
                            if currmapspatinfo <1;
                                color = 'b';
                            elseif currmapspatinfo <3;
                                color = [0 .5 0];
                            elseif currmapspatinfo <10;
                                color = 'r';
                            elseif currmapspatinfo >=10;
                                color = 'k';
                                
                            end
                            %                             title([num2str(currmapspatinfo)], 'FontSize',24,'Fontweight','bold', 'Color', color);
                            title({['SI=',num2str(currmapspatinfo)]; ['maxF=',num2str(max(max(allmapdata{i}{k}.smoothedspikerate)))]; [sprintf('D%d E%d T%d C%d',pfm{1}(i).output{1}(k).index(1),pfm{1}(i).output{1}(k).index(2),pfm{1}(i).output{1}(k).index(3),pfm{1}(i).output{1}(k).index(4))]}, 'FontSize',24,'Fontweight','bold', 'Color', color);
                            maxratesmod(k) = max(max(allmapdata{i}{k}.smoothedspikerate));
                            spatinfoallmod(k) = currmapspatinfo;
                        else
                            title(['~exist'], 'FontSize',24,'Fontweight','bold');
                        end

%                         xlabel(typecell, 'Color', typclr)
                        set(gca,'xtick',[])
                        set(gca,'xticklabel',[])
                        set(gca,'ytick',[])
                        set(gca,'yticklabel',[])
                        % Update plot number
                        totalplots=totalplots+1;
                        if (rem(totalplots,20) >4 & rem(totalplots,20) <10) | (rem(totalplots,20) >14);
                            totalplots = totalplots+5; %skip the firing rate maps lines
                        end
                        if rem(totalmapplots+1,20)==0

                            %save the full maps figures
                            figfile = [anfigdir,'/','ModMaps',num2str(figcnt)];
                            [ax,hU]=suplabel(titstring,'t', [.08 .08 .84 .84]);
                            set(hU,'FontSize',200,'Fontweight','bold')  %font controls aren't working.. i think it's a bug with ubuntu matlab2013a...
                            if savefigs==1
                                print('-djpeg', figfile);
%                  print('-dpdf', figfile);
                                %saveas(gcf,figfile,'fig');
%                                 print('-depsc2', figfile);
                            end
                            
                            if ~cyclemaps == 0
                                keyboard;
                            end
                            close
                        end
                    end
                    
                    %save the last maps figure as it's most likely
                    %not full
                    figfile = [anfigdir,'/','ModMaps',num2str(figcnt)];
                    [ax,hU]=suplabel(titstring,'t', [.08 .08 .84 .84]);
                    set(hU,'FontSize',200,'Fontweight','bold')  %font controls aren't working.. i think it's a bug with ubuntu matlab2013a...
                    if savefigs==1
                        print('-djpeg', figfile);
%                  print('-dpdf', figfile);
%                         saveas(gcf,figfile,'pdf')
%                         print('-depsc2', figfile);
                    end
                    if~cyclemaps ==0
                        keyboard
                    end
                    close
                end
                
                %___________________________________________
                if plotunmodmaps == 1;
                    trajdata = [];
                    figcnt = 0; totalplots = 0; totalmapplots=0;%Count total plots across figures
                    clr = {'b','r','g','m','c','y','k','r'};
                    for k=1:length(pfm{2}(i).output{1}); %unmodulated
                        %__________TRAJS
                        alltrajdata{i}{k}=pfs{2}(i).output{1}(k);
                        if rem(totalmapplots+1,20)==0; %  subplots finished or first plot
                            figcnt=figcnt+1;
                        end
                        figure(figcnt+1);  hold on;
                        
                        subplot(4,5,rem(totalplots,20)+1); hold on;
                        trajdata = alltrajdata{i}{k}.trajdata;
                        for o=1:length(trajdata),
                            plot(trajdata{o}(:,5),[clr{o} '.-'],'Linewidth',2);
                            %maxrate(i) = max(trajdata{i}(5:end-5,5));
                        end
                        if rem(totalplots+1,20)==1;
                            legend('OL','IL','OR','IR');
                            xlabel ('Position along linear trajectory (cm)','FontSize',14,'Fontweight','bold');
                            ylabel ('Firing Rate (Hz)','FontSize',14,'Fontweight','bold');
                            
                            titstring = [animals{i}, '________UnMod______', tit, '_____cell group:', num2str(figcnt)];
                            
                        end
                        
                        
                        set(gca,'xtick',[])
                        set(gca,'xticklabel',[])
                        set(gca,'ytick',[])
                        set(gca,'yticklabel',[])
                        
                        
                        %__________________________FIRING MAPS
                        allmapdata{i}{k} = pfm{2}(i).output{1}(k); %unmod for now
                        totalmapplots = totalplots+5;
                        subplot(4,5,rem(totalmapplots,20)+1); hold on;
                        imagesc(allmapdata{i}{k}.smoothedspikerate')
                        %                         set(gca,'YLim',[0 110]);
                        %                         set(gca,'XLim',[0 117]);
                        
                        titstring = [animals{i}, '________Unmod______', tit, '_____cell group:', num2str(figcnt)];
                        
                        
                        
                        %                         tmpspatind = spatunmoddata(:,[1 2 3 5]); % old used when spat info func was spitting out dat ep tet traj cell
                        tmpspatind = spatunmoddata_alltraj(:,[1:4]);
                        %find spatial info for map and make it title..
                        %color it according to intensity
                        if ~rowfind(allmapdata{i}{k}.index,tmpspatind) == 0;
                            currmapspatinfo = spatunmoddata_alltraj(rowfind(allmapdata{i}{k}.index,tmpspatind),6);
                            if currmapspatinfo <1;
                                color = 'b';
                            elseif currmapspatinfo <3;
                                color = [0 .5 0];
                            elseif currmapspatinfo <10;
                                color = 'r';
                            elseif currmapspatinfo >=10;
                                color = 'k';
                            end
                            %                             title({[num2str(pfm{1}(i).output{1}(k).index)]; ['SI=',num2str(currmapspatinfo)]; ['maxF=',num2str(max(max(allmapdata{i}{k}.smoothedspikerate)))]}, 'FontSize',24,'Fontweight','bold', 'Color', color);
                            title({['SI=',num2str(currmapspatinfo)]; ['maxF=',num2str(max(max(allmapdata{i}{k}.smoothedspikerate)))]; [sprintf('D%d E%d T%d C%d',pfm{2}(i).output{1}(k).index(1),pfm{2}(i).output{1}(k).index(2),pfm{2}(i).output{1}(k).index(3),pfm{2}(i).output{1}(k).index(4))]}, 'FontSize',24,'Fontweight','bold', 'Color', color);
                            maxratesunmod(k) = max(max(allmapdata{i}{k}.smoothedspikerate));
                            spatinfoallunmod(k) = currmapspatinfo;
                            
                        else
                            title(['~exist'], 'FontSize',24,'Fontweight','bold');
                        end
                        set(gca,'xtick',[])
                        set(gca,'xticklabel',[])
                        set(gca,'ytick',[])
                        set(gca,'yticklabel',[])
                        % Update plot number
                        totalplots=totalplots+1;
                        if (rem(totalplots,20) >4 & rem(totalplots,20) <10) | (rem(totalplots,20) >14);
                            totalplots = totalplots+5; %skip the firing rate maps lines
                        end
                        if rem(totalmapplots+1,20)==0
                            
                            
                            %save the full maps figures
                            figfile = [anfigdir,'/','UnmodMaps',num2str(figcnt)];
                            
                            [ax,hU]=suplabel(titstring,'t', [.08 .08 .84 .84]);
                            set(hU,'FontSize',200,'Fontweight','bold')  %font controls aren't working.. i think it's a bug with ubuntu matlab2013a...
                            if savefigs==1
                                print('-djpeg', figfile);
                                %saveas(gcf,figfile,'fig');
                                %print('-depsc2', figfile);
                            end
                            
                            if ~cyclemaps == 0
                                keyboard;
                            end
                            close
                        end
                    end
                    
                    %save the last map figure as it's most likely
                    %not full
                    figfile = [anfigdir,'/','UnmodMaps',num2str(figcnt)];
                    
                    [ax,hU]=suplabel(titstring,'t', [.08 .08 .84 .84]);
                    set(hU,'FontSize',200,'Fontweight','bold')  %font controls aren't working.. i think it's a bug with ubuntu matlab2013a...
                    if savefigs==1
                        print('-djpeg', figfile);
                        %saveas(gcf,figfile,'fig');
                        %print('-depsc2', figfile);
                    end
                    if ~cyclemaps == 0
                        keyboard;
                    end
                    close
                end
            end
            
            % plot max rates by spat info for each an
            if plotrateXinfo == 1;
            figure
%             subplot(1,2,1); 
            hold on;
            scatter(maxratesmod, spatinfoallmod, 'MarkerFaceColor', 'b');
            modfit = polyfit(maxratesmod, spatinfoallmod, 1);
            modfitoutput = linspace(min(maxratesmod), max(maxratesmod))*modfit(1,1) + modfit(1,2);
            plot(modfitoutput, 'b');
%             xlabel('firing rate');
%             ylabel('spatial info');
%             title({'modulated'; sprintf('slope = %d',modfit(1,1))});
            %                      subplot(1,2,2); 
            hold on;
            scatter(maxratesunmod, spatinfoallunmod, 'r', 'MarkerFaceColor', 'r');
            unmodfit = polyfit(maxratesunmod, spatinfoallunmod, 1);
            unmodfitoutput = linspace(min(maxratesunmod), max(maxratesunmod))*unmodfit(1,1) + unmodfit(1,2);
            plot(unmodfitoutput, 'r');
            xlabel('firing rate');
            ylabel('spatial info');
            title([animals{i}, '________max rate X Spatial Info'])
            %             title({'unmodulated'; sprintf('slope = %d',unmodfit(1,1))});
            legfix2 = legend(['mod',sprintf('slope = %d',modfit(1,1))], ['unmod',sprintf('slope = %d',unmodfit(1,1))]);
            
            set(legfix2, 'Color', 'none');
            
            linesInPlot = findobj(legfix2, 'type', 'line');
            set(linesInPlot(:),'XData',[0.7 0.7]);
            
            redInPlot = findobj(legfix2, 'Color', 'r');
            blueInPlot = findobj(legfix2, 'Color', 'b');
            if ~isempty(redInPlot)
            set(redInPlot(:),'color', [1 1 1]);
            end
            if ~isempty(blueInPlot)
            set(blueInPlot(:),'color', [1 1 1]);
            end

            hText = findobj(legfix2, 'type', 'text');
            set(hText(2),'color',  [0 0 1]);
            set(hText(1),'color',  [1 0 0]);
            %                             set(hText(2),'color',  [0 1 0]);
            %                             set(hText(1),'color',  [1 0 1]);

            
%             [ax,hU]=suplabel('max rate X Spatial Info','t', [.08 .08 .84 .9]);
%             set(hU,'FontSize',200,'Fontweight','bold')
            
            figfile = [anfigdir,'/','RateXInfo'];
            if savefigs==1
                print('-djpeg', figfile);
                %saveas(gcf,figfile,'fig');
                %print('-depsc2', figfile);
            end
            if ~cyclemaps == 0
                keyboard;
            end
            close
            end
            

            
            
            
        end
    end
    

    
    %recompute data across animals__________________________________________________________________________________________________________________________________
    
    if plotcomban == 1
        for i = 1: length(animals);
            %combine the spat info of all trajs within each ep for each cell
            for mu = 1:2;
                dta = spatinfo{1,mu}(1,i).output{1};
                tmpnum = nan(length(dta(:,1)),1);
                for hu = 1:length(dta(:,1));
                    x = dta(hu,1:4);
                    y = sprintf('%d',x);
                    z = str2num(y);
                    tmpnum(hu) = z;
                end
                [c ia ic] = unique(tmpnum, 'stable');
                spatdata_alltraj = [];
                for gu = 1:length(ia);
                    %                     if gu < length(ia); %if it's not the last one
                    spatdata_alltraj(gu,6) = nansum(dta(find(ic == gu),6));
                    %                     spatdata_alltraj(gu,6) = nansum(dta(ia(gu):(ia(gu+1)-1),6)); %TAKING THE SUM OF TRAJS.. the trajs of the same cell will be in order which is why i can just take the sum of the range
                    spatdata_alltraj(gu,1:4) = dta(ia(gu),1:4); %take the index of first spt info number %leaving the ep column (2) empty
                    %                     else
                    %                         spatdata_alltraj(gu,6) = nansum(dta(ia(gu):end,6));
                    %                         spatdata_alltraj(gu,1:4) = dta(ia(gu),1:4);
                    %                     end
                end
                %_________________
                if collapseEPinday == 1;
                    %combine the spat info of all trajs within each ep for each cell
                    tmpnum = nan(length(spatdata_alltraj(:,1)),1);
                    for hu = 1:length(spatdata_alltraj(:,1));
                        x = spatdata_alltraj(hu,[1 3 4]);
                        y = sprintf('%d',x);
                        z = str2num(y);
                        tmpnum(hu) = z;
                    end
                    [c ia ic] = unique(tmpnum, 'stable');
                    spatdata_alltrajep = [];
                    for gu = 1:length(ia);
                        %                         spatdata_alltrajep(gu,6) = nanmean(spatdata_alltraj(find(ic == gu),6));
                        spatdata_alltrajep(gu,6) = nanmean(spatdata_alltraj(find(ic == gu),6)); %TAKING THE MEAN OF EPS.. eps of same cell aren't in order so i need to do a find on ic
                        spatdata_alltrajep(gu,[1 3 4]) = spatdata_alltraj(ia(gu),[1 3 4]); %take the index of first spt info number %leaving the ep column (2) empty
                    end
                end
                
                
                %_________________
                if mu == 1;
                    spatmoddata = spatdata_alltrajep;
                    spatmoddata_alltraj = spatdata_alltraj; %use this for plotting maps with epochs separate but trajs combined
                else
                    spatunmoddata = spatdata_alltrajep;
                    spatunmoddata_alltraj = spatdata_alltraj;  %use this for plotting maps with epochs separate but trajs combined
                end
            end
            %                 spatmoddata = (spatinfo{1,1}(1,i).output{1}); %houskeeping
            %                 spatunmoddata = (spatinfo{1,2}(1,i).output{1});
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
        bar(1, mean(mod), 'b','EdgeColor','none')
        bar(2, mean(unmod), 'r','EdgeColor','none')
        errorbar2([1 2], [mean(mod) mean(unmod)],  [stderr(mod) stderr(unmod)] , 0.3, 'k')
        xlim([0.3 2.7])
        set(gca, 'fontsize', 24)
        set(gca, 'xtick', [1 2], 'xticklabel', {'mod', 'unmod'})
        ylabel('info (bits/spike)')
        rp = ranksum(mod,unmod);
        title({sprintf('Spatial Info All Animals n(%d) %s',length(animals),tit); sprintf('ranksum p = %s',rp); sprintf('K-Stest p = %s',p); sprintf('n mod= %d unmod= %d',length(mod),length(unmod))})
        
        
    
    figfile = [figdir,currfigdir, 'AllInfo'];
    if savefigs==1
        print('-djpeg', figfile);
        %saveas(gcf,figfile,'fig');
        %print('-depsc2', figfile);
    end
    if ~cyclemaps == 0
        keyboard;
    end
    close
    
    %combined field
    %sizes__________________________________________________________________________________________________________________________________
    for i = 1: length(animals);
        
        for mu = 1:2;
            clear  tmpnum x y z c ia ic spatdata_alltrajep gu indexdata meanday
            currindexmod=[]; cntallcells = 0; total_fracundermean = []; pu =0;
            nidxs_an = length(plaf{mu}(i).output{1});
            for wu=1:nidxs_an
                indexdata(wu,:) = [plaf{mu}(i).output{1}(wu).index plaf{mu}(i).output{1}(wu).total_fracundermean];
                currindexmod(wu,:) = plaf{mu}(i).output{1}(wu).index;
            end
            tmpnum = nan(length(currindexmod(:,1)),1);
            for hu = 1:length(currindexmod(:,1));
                x = currindexmod(hu,[1 3 4]);
                y = sprintf('%d',x);
                z = str2num(y);
                tmpnum(hu) = z;
            end
            [c ia ic] = unique(tmpnum, 'stable');
            sizedata_allep = [];
            for gu = 1:length(ia);
                sizedata_allep(gu,5) = nanmean(indexdata((find(ic == gu)),5));
                sizedata_allep(gu,[1 3 4]) = indexdata((ia(gu)),[1 3 4]);
            end
            if mu == 1;
                sizefracmeanmoddata = sizedata_allep;
            else
                sizefracmeanunmoddata = sizedata_allep;
            end
        end
        
        if i ==1
            tit = 'All W Epochs';
            mod = sizefracmeanmoddata(:,5); 
            unmod = sizefracmeanunmoddata(:,5); 
        else
            mod = stack(mod,sizefracmeanmoddata(:,5)); 
            unmod = stack(unmod,sizefracmeanunmoddata(:,5)); 

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
bar(1, mean(mod), 'b','EdgeColor','none')
bar(2, mean(unmod), 'r','EdgeColor','none')
errorbar2([1 2], [mean(mod) mean(unmod)],  [stderr(mod) stderr(unmod)] , 0.3, 'k')
xlim([0.3 2.7])
set(gca, 'fontsize', 24)
set(gca, 'xtick', [1 2], 'xticklabel', {'mod', 'unmod'})
ylabel('Proportion track active above mean')
rp = ranksum(mod,unmod);
title({sprintf('Field Size (Above Cell Mean rate)___All Animals n(%d) %s',length(animals),tit); sprintf('ranksum p = %s',rp); sprintf('K-Stest p = %s',p); sprintf('n mod= %d unmod= %d',length(mod),length(unmod))})


figfile = [figdir,currfigdir, 'Allsize'];
if savefigs==1
    print('-djpeg', figfile);
    %saveas(gcf,figfile,'fig');
    %print('-depsc2', figfile);
end
if ~cyclemaps == 0
    keyboard;
end
close


%_______________________________________________________________________________________________________________________________________________________
%corr coefficient combined
for i = 1: length(animals);
                for k = 1:length(pfs{1}(i).output{1}(1,:));
                    currind = [];
                    if length(pfs{1}(i).output{1}(1,k).trajdata) ==4;  % if all trajs exist
                        
                        %mod gather trajs
                        currind = pfs{1}(i).output{1}(1,k).index;
                        
                        outleftmod = pfs{1}(i).output{1}(1,k).trajdata{1,1}(:,5);
                        outrightmod= pfs{1}(i).output{1}(1,k).trajdata{1,3}(:,5);
                        
                        inleftmod = pfs{1}(i).output{1}(1,k).trajdata{1,2}(:,5);
                        inrightmod = pfs{1}(i).output{1}(1,k).trajdata{1,4}(:,5);
                        
                        %compute corrcoef
                        outmodmintrunc = min(length(outleftmod), length(outrightmod));
                        [outmodR outmodP] = corrcoef(outrightmod(1:outmodmintrunc),outleftmod(1:outmodmintrunc),'rows','pairwise'); % rows-pairwise ~ compute R(i,j) using rows with no NaN values in either column i or j.
                        %             outmodAllR(k,5)  = outmodR(2,1);
                        %             outmodAllR(k,[1 3 4]) = currind(1,[1 3 4]);
                        AllR{1}(k,5)  = outmodR(2,1);
                        AllR{1}(k,[1 3 4]) = currind(1,[1 3 4]);
                        outmodAllP(k,5) = outmodP(2,1);
                        outmodAllP(k,[1 2 3 4]) = currind(1,[1 2 3 4]);
                        
                        inmodmintrunc = min(length(inleftmod), length(inrightmod));
                        [inmodR inmodP] = corrcoef(inrightmod(1:inmodmintrunc),inleftmod(1:inmodmintrunc),'rows','pairwise'); % rows-pairwise ~ compute R(i,j) using rows with no NaN values in either column i or j.
                        %             inmodAllR(k,5)  = inmodR(2,1);
                        %             inmodAllR(k,[1 3 4]) = currind(1,[1 3 4]);
                        AllR{2}(k,5)  = inmodR(2,1);
                        AllR{2}(k,[1 3 4]) = currind(1,[1 3 4]);
                        inmodAllP(k,5) = inmodP(2,1);
                        inmodAllP(k,[1 2 3 4]) = currind(1,[1 2 3 4]);
                        
                        
                    else
                        pfs{1}(i).output{1}(1,k).index %print any cells skipped
                        
                    end
                end
                
                %UNMOD________________
                for k = 1:length(pfs{2}(i).output{1}(1,:));
                    currind = [];
                    if length(pfs{2}(i).output{1}(1,k).trajdata) ==4;  % if all trajs exist
                        
                        %gather trajs
                        currind = pfs{2}(i).output{1}(1,k).index;
                        
                        outleftunmod = pfs{2}(i).output{1}(1,k).trajdata{1,1}(:,5);
                        outrightunmod= pfs{2}(i).output{1}(1,k).trajdata{1,3}(:,5);
                        
                        inleftunmod = pfs{2}(i).output{1}(1,k).trajdata{1,2}(:,5);
                        inrightunmod = pfs{2}(i).output{1}(1,k).trajdata{1,4}(:,5);
                        
                        %compute corrcoef
                        outunmodmintrunc = min(length(outleftunmod), length(outrightunmod));
                        [outunmodR outunmodP] = corrcoef(outrightunmod(1:outunmodmintrunc),outleftunmod(1:outunmodmintrunc),'rows','pairwise'); % rows-pairwise ~ compute R(i,j) using rows with no NaN values in either column i or j.
                        %             outunmodAllR(k,5)  = outunmodR(2,1);
                        %             outunmodAllR(k,[1 3 4]) = currind(1,[1 3 4]);
                        AllR{3}(k,5)  = outunmodR(2,1);
                        AllR{3}(k,[1 3 4]) = currind(1,[1 3 4]);
                        outunmodAllP(k,5) = outunmodP(2,1);
                        outunmodAllP(k,[1 2 3 4]) = currind(1,[1 2 3 4]);
                        
                        inunmodmintrunc = min(length(inleftunmod), length(inrightunmod));
                        [inunmodR inunmodP] = corrcoef(inrightunmod(1:inunmodmintrunc),inleftunmod(1:inunmodmintrunc),'rows','pairwise'); % rows-pairwise ~ compute R(i,j) using rows with no NaN values in either column i or j.
                        %             inunmodAllR(k,5)  = inunmodR(2,1);
                        %             inunmodAllR(k,[1 3 4]) = currind(1,[1 3 4]);
                        AllR{4}(k,5)  = inunmodR(2,1);
                        AllR{4}(k,[1 3 4]) = currind(1,[1 3 4]);
                        inunmodAllP(k,5) = inunmodP(2,1);
                        inunmodAllP(k,[1 2 3 4]) = currind(1,[1 2 3 4]);
                        
                    else
                        pfs{2}(i).output{1}(1,k).index %print any cells skipped
                        
                    end
                end
                
                %day chunking
                
                for mu = 1:length(AllR);
                    clear  tmpnum x y z c ia ic gu
                    currindexmod=[]; cntallcells = 0; total_fracunder3 = []; pu =0;
                    nidxs_an = length(AllR{mu});
                    tmpnum = nan(length(AllR{mu}(:,1)),1);
                    for hu = 1:length(AllR{mu}(:,1));
                        x = AllR{mu}(hu,[1 3 4]);
                        y = sprintf('%d',x);
                        z = str2num(y);
                        tmpnum(hu) = z;
                    end
                    [c ia ic] = unique(tmpnum, 'stable');
                    for gu = 1:length(ia);
                        corrdata_allep{mu}(gu,5) = nanmean(AllR{mu}((find(ic == gu)),5));
                        corrdata_allep{mu}(gu,[1 3 4]) = AllR{mu}((ia(gu)),[1 3 4]);
                    end
                end
                
                if i ==1
                    tit = 'All W Epochs';
                    outmod = corrdata_allep{1}(:,5);
                    outunmod = corrdata_allep{3}(:,5);
                    inmod = corrdata_allep{2}(:,5);
                    inunmod = corrdata_allep{4}(:,5);
                else
                    outmod = stack(outmod, corrdata_allep{1}(:,5));
                    outunmod = stack(outunmod, corrdata_allep{3}(:,5));
                    inmod = stack(inmod, corrdata_allep{2}(:,5));
                    inunmod = stack(inunmod, corrdata_allep{4}(:,5));
                end

end

[h p] = kstest2(outmod(~isnan(outmod)),outunmod(~isnan(outunmod)));
medoutmod = median(outmod(~isnan(outmod)));
medoutunmod = median(outunmod(~isnan(outunmod)));
sizes = [size(outmod(~isnan(outmod)),1) size(outunmod(~isnan(outunmod)),1) ];
outmod = outmod(~isnan(outmod));
outunmod = outunmod(~isnan(outunmod));
figure
subplot(1,2,1)
hold on
bar(1, mean(outmod), 'b','EdgeColor','none')
bar(2, mean(outunmod), 'r','EdgeColor','none')
errorbar2([1 2], [mean(outmod) mean(outunmod)],  [stderr(outmod) stderr(outunmod)] , 0.3, 'k')
xlim([0.3 2.7])
set(gca, 'fontsize', 24)
set(gca, 'xtick', [1 2], 'xticklabel', {'mod', 'unmod'})
ylabel('Corr Coef R')
rp = ranksum(outmod,outunmod);
title({sprintf('CorrCoef Outbound___All Animals n(%d) %s',length(animals),tit); sprintf('ranksum p = %s',rp); sprintf('K-Stest p = %s',p); sprintf('n mod= %d unmod= %d',length(outmod),length(outunmod))})


[h p] = kstest2(inmod(~isnan(inmod)),inunmod(~isnan(inunmod)));
medinmod = median(inmod(~isnan(inmod)));
medinunmod = median(inunmod(~isnan(inunmod)));
sizes = [size(inmod(~isnan(inmod)),1) size(inunmod(~isnan(inunmod)),1) ];
inmod = inmod(~isnan(inmod));
inunmod = inunmod(~isnan(inunmod));
subplot(1,2,2)
hold on
bar(1, mean(inmod), 'b','EdgeColor','none')
bar(2, mean(inunmod), 'r','EdgeColor','none')
errorbar2([1 2], [mean(inmod) mean(inunmod)],  [stderr(inmod) stderr(inunmod)] , 0.3, 'k')
xlim([0.3 2.7])
set(gca, 'fontsize', 24)
set(gca, 'xtick', [1 2], 'xticklabel', {'mod', 'unmod'})
ylabel('Corr Coef R')
rp = ranksum(inmod,inunmod);
title({sprintf('CorrCoef Inbound___All Animals n(%d) %s',length(animals),tit); sprintf('ranksum p = %s',rp); sprintf('K-Stest p = %s',p); sprintf('n mod= %d unmod= %d',length(inmod),length(inunmod))})

figfile = [figdir,currfigdir, 'Allcorrcoef'];
if savefigs==1
    print('-djpeg', figfile);
    %saveas(gcf,figfile,'fig');
    %print('-depsc2', figfile);
end
if ~cyclemaps == 0
    keyboard;
end
close


end
end


