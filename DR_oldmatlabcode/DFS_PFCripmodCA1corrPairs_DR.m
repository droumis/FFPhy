
warning('off','all');
clear; close all;

runscript = 1;
savedata = 1; % save data option - only works if runscript is also on
savefigs=0;
plotstuff = 1; %
plotEachcell = 1;
cyclemaps = 1;
% runnadal =0;

plotPFCCA1maps =1;
savedir = '/mnt/data25/sjadhav/HPExpt/ProcessedDataDR/';
savefilename = 'Nadal_sigcorrPair';
savefile = [savedir savefilename]; area = 'PFC'; %clrunmod = 'r'; clrmod = 'b'; % PFC
figdir = '/mnt/data25/sjadhav/HPExpt/Figures_DR/';

Veqn = '>3';
minV=str2num(Veqn(end));
mintime = 3;
traj = [1:4] ;

% If runscript, run Datafilter and save data
if runscript == 1
    %     for i =  modUnmod;
    %Animal selection
    %-----------------------------------------------------
%                 animals = {'HPa' 'HPb' 'HPc' 'nadal'};
%             animals = {'HPa' 'HPb' 'HPc'};
    %         animals = {'HPc'};
%     animals = {'HPa'};
            animals = {'nadal'};
    
    %Filter creation
    %-----------------------------------------------------
    % Epoch filter
    % -------------
    %     if runnadal == 1;
            dayfilter = '8:17';
    %     else
%     dayfilter = '1:8';
    %     end
    
    %     if runnadal == 1;
            runepochfilter = 'isequal($type, ''run'')';
    %     else
%     runepochfilter = 'isequal($environment, ''wtr1'') || isequal($environment, ''wtr2'')';
    %     end
    
    % %Cell filter
    % %-----------
    %     placecellfilter = '(strcmp($area, ''PFC'') || (strcmp($area, ''CA1'') || && ($numspikes > 100))';  % not mod/unmod
    %     placecellfilter = '(strcmp($area, ''PFC'') && ($numspikes > 100) || strcmp($area, ''CA1'') && ($numspikes > 100))';
%     placecellfilter = '( strcmp($tag, ''CA1Pyr'') || strcmp($tag, ''iCA1Pyr'') || strcmp($tag, ''PFC'')) && ($numspikes > 100)';
    placecellfilter = '( strcmp($area, ''CA1'') || strcmp($area, ''iCA1'') || strcmp($area, ''PFC'')) && ($numspikes > 100)';
   
    %         if i < 2;
    %                         placecellfilter = '(strcmp($area, ''PFC'') && strcmp($ripmodtag, ''y'') && ($numspikes > 100))';   % Ripple mod
    % %             placecellfilter = '(strcmp($area, ''PFC'') && strcmp($thetamodtag, ''y'') && ($numspikes > 100))';   % theta mod
    %
    %         else
    %                         placecellfilter = '(strcmp($area, ''PFC'') && strcmp($ripmodtag, ''n'') && ($numspikes > 100))'; % Ripple unmod
    % %             placecellfilter = '(strcmp($area, ''PFC'') && strcmp($thetamodtag, ''n'') && ($numspikes > 100))'; % theta unmod
    %         end
    
    % Time filter -
    %%-----------
    
    riptetfilter = '(isequal($descrip, ''riptet''))';
    timefilter = { {'DFTFsj_getlinstate', '(($state ~= -1) & (abs($linearvel) >= 3))', 6}, {'DFTFsj_getriptimes','($nripples == 0)','tetfilter',riptetfilter,'minthresh',3} }; %DR added velocity filter.. trying to get ride of v high prococc in data
    
    % Iterator
    % --------
    iterator = 'singlecellanal';
    
    % Filter creation
    % ----------------
    %     spatf = createfilter('animal',animals,'epochs',runepochfilter,'cells',placecellfilter,'excludetime', timefilter,'iterator',iterator);
    spatf = createfilter('animal',animals,'days',dayfilter,'epochs',runepochfilter,'cells',placecellfilter,'excludetime', timefilter,'iterator',iterator);
    
    %do i need this?
    spatf = testexcludetimes(spatf, mintime); %removes epochs from analysis if all epoch excluded by excludetimes, mintime = 30
    
    disp('Done Filter Creation');
    
    % Set analysis function
    % ----------------------
    %use_____________________
    psf = setfilterfunction(spatf, 'DFAsj_filtercalclinfields_tf',{'spikes', 'linpos'}, 'binsize', 2);
    pmf = setfilterfunction(spatf, 'DFAsj_openfieldrate_tf',{'spikes', 'linpos', 'pos'}, 'binsize', 1, 'std', 2);
    
    % Run analysis-----------------------
    %         pfs{i} = runfilter(pfs);  % Place Field Stability.. trajectories
    pfs = runfilter(psf);  % Place Field Stability.. trajectories
    pfm = runfilter(pmf);  % Place Field Map
    
    %     end
    disp('Finished running filter script');
    %--------------------- Finished Filter Function Run -------------------
    
    if savedata == 1
        clear runscript  savedata plotstuff cyclemaps savefigs plottrajs runnadal plotEachcell plotPFCCA1maps
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

% --------------------------------------------------------------------

%make maps for PFC-CA1 cells sig predictive (glm) pairwise
%correlative
if savefigs == 1;
    mkdir(figdir,savefilename)
end
if plotPFCCA1maps ==1;
    load('/mnt/data25/sjadhav/HPExpt/ProcessedData/HP_allPFCCA1sigcorridxs');
    pairdata = allPFCCA1sigcorridxs;
    pfclistidx = []; ju = 0; pfsLIST = []; pu =0;
    
    for yu = 1:length(pairdata);
        pfclistidx(yu) = pairdata(1,yu).PFCidx(:);
    end

    clear dta tmpnum hu x y z c ia ic spatdata_alltraj gu
    dta = pfclistidx;
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
    
    
        
    for pu = 1:length(pfs.output{1,1});
    pfsLIST(pu, :) = pfs.output{1}(1,pu).index;
    end
    
    
    
    
    
    figcnt = 0; totalplots = 0; totalmapplots=1; maxrates = []; totaltrajplots = 1;
    for ju = 1:length(pfclistidx(:,1)); %for each cell
        currPFCtrajdata = [];  currCA1ind = [];
        clr = {'b','r','g','m','c','y','k','r'};
        allCA1ind = pairdata(1,ju).CA1sigidxs(:,:); % each CA1 cell

        currPFCindex = pfclistidx(ju,:);
        currPFCtrajdata = pfs(currPFCindex(1,1)).output{1}(ju).trajdata;
        %         if rem(totaltrajplots-1,length(allCA1ind)+1)==0; %  subplots finished
        %             figcnt=figcnt+1;
        %         end
        %         figure(figcnt+1);
        
        subplot(2,length(allCA1ind)+1,totaltrajplots); hold on; %alter to fit various ca1 data size
        for o=1:length(currPFCtrajdata),
            plot(currPFCtrajdata{o}(:,5),[clr{o} '.-'],'Linewidth',2);
        end
        %         if rem(totalplots+1,20)==1;
        %             xlabel ('Position along linear trajectory (cm)','FontSize',14,'Fontweight','bold');
        %             ylabel ('Firing Rate (Hz)','FontSize',14,'Fontweight','bold');
        %             titstring = [animals{i}, '________Mod______', tit, '_____cell group:', num2str(figcnt)];
        %         end
        set(gca,'xtick',[])
        set(gca,'xticklabel',[])
        set(gca,'ytick',[])
        set(gca,'yticklabel',[])
        
        %________________________ FIRING RATE MAP for the PFC cell
        currPFCmapdata = pfm(currPFCindex(1,1)).output{1}(ju).smoothedspikerate;
        totalmapplots = totaltrajplots+length(allCA1ind)+1;
        subplot(2,length(allCA1ind)+1,totalmapplots); hold on; %alter for CA1
        imagesc(currPFCmapdata')
        
        %          if (rem(totalplots,20) >4 & rem(totalplots,20) <10) | (rem(totalplots,20) >14);
        %             totalplots = totalplots+5; %skip the firing rate maps lines
        %         end
        
        for tu = 1:length(allCA1ind);  %traj map for the CA1 cell%_____________________________________________________________________________
            totaltrajplots = totaltrajplots +1;
            
            currCA1ind = pairdata(1,ju).CA1sigidxs(tu,:); % each CA1 cell
            currCA1trajdata = pfs.output{1}(ju).trajdata;
            
            subplot(2,length(currCA1ind)+1,totaltrajplots); hold on; %alter to fit various ca1 data size
            for o=1:length(currPFCtrajdata),
                plot(currPFCtrajdata{o}(:,5),[clr{o} '.-'],'Linewidth',2);
            end
            %         if rem(totalplots+1,20)==1;
            %             xlabel ('Position along linear trajectory (cm)','FontSize',14,'Fontweight','bold');
            %             ylabel ('Firing Rate (Hz)','FontSize',14,'Fontweight','bold');
            %             titstring = [animals{i}, '________Mod______', tit, '_____cell group:', num2str(figcnt)];
            %         end
            set(gca,'xtick',[])
            set(gca,'xticklabel',[])
            set(gca,'ytick',[])
            set(gca,'yticklabel',[])
            
            %________________________ FIRING RATE MAP for the PFC cell
            currPFCmapdata = pfm(currPFCindex(1,1)).output{1}(ju).smoothedspikerate;
            totalmapplots = totaltrajplots+length(currCA1ind)+1;
            subplot(2,length(currCA1ind)+1,totalmapplots); hold on; %alter for CA1
            imagesc(currPFCmapdata')
            
            %          if (rem(totalplots,20) >4 & rem(totalplots,20) <10) | (rem(totalplots,20) >14);
            %             totalplots = totalplots+5; %skip the firing rate maps lines
            
        end
        
        
        
        
        
        
        
        
        
        
        
        if rem(totalmapplots,(length(currCA1ind)+1)*2)==0
            totalmapplots = 0; totaltrajplots = 0;
            %save the full maps figures
            figfile = [figdir,savefilename,'/','SigCorrMaps',num2str(figcnt)];
            %             [ax,hU]=suplabel(titstring,'t', [.08 .08 .84 .84]);
            %             set(hU,'FontSize',200,'Fontweight','bold')  %font controls aren't working.. i think it's a bug with ubuntu matlab2013a...
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
    
    
    
    
    
    
    
end
%     load([HPdir,sprintf('%s_direct/%scellinfo.mat',animals{1,i},animals{1,i})]);
%     trajdata = [];
%     figcnt = 0; totalplots = 0; totalmapplots=0;%Count total plots across figures
%     clr = {'b','r','g','m','c','y','k','r'};
%     figcnt = 0; totalplots = 0; maxrates = [];%Count total plots across figures
%     for k=1:length(pfm{1}(i).output{1}); %modulated
%         %______________________
%
%         alltrajdata{i}{k}=pfs{1}(i).output{1}(k);
%         if rem(totalmapplots+1,20)==0; %  subplots finished
%             figcnt=figcnt+1;
%         end
%         figure(figcnt+1);  hold on;
%         subplot(4,5,rem(totalplots,20)+1); hold on;
%         trajdata = alltrajdata{i}{k}.trajdata;
%         for o=1:length(trajdata),
%             plot(trajdata{o}(:,5),[clr{o} '.-'],'Linewidth',2);
%             %                             maxrate(o) = max(trajdata{o}(5:end-5,5));
%         end
%         typeind = pfm{1}(i).output{1}(k).index;
%         typecell = cellinfo{typeind(1)}{typeind(2)}{typeind(3)}{typeind(4)}.typetag;
%
%         if typecell == 'exc';
%             typclr = [0 .7 .93];
%         elseif typecell == 'inh';
%             typclr = [.87 1 .65];
%         else
%             typclr = [.8 .8 .8];
%         end
%         title(typecell,'Color',typclr)
%         set(gca,'Color',typclr);
%         set(gcf, 'InvertHardCopy', 'off');
%         if rem(totalplots+1,20)==1;
%
%             xlabel ('Position along linear trajectory (cm)','FontSize',14,'Fontweight','bold');
%             ylabel ('Firing Rate (Hz)','FontSize',14,'Fontweight','bold');
%             %                             legend('OL','IL','OR','IR');
%             legfix = legend('OL','IL','OR','IR');
%             set(legfix, 'Color', 'none');
%             hText = findobj(legfix, 'type', 'text');
%             set(hText(4),'color',  [0 0 1]);
%             set(hText(3),'color',  [1 0 0]);
%             set(hText(2),'color',  [0 1 0]);
%             set(hText(1),'color',  [1 0 1]);
%             linesInPlot = findobj(legfix, 'type', 'line');
%             set(linesInPlot(1:8),'XData',[0.5 0.7]);
%
%
%             titstring = [animals{i}, '________Mod______', tit, '_____cell group:', num2str(figcnt)];
%         end
%         set(gca,'xtick',[])
%         set(gca,'xticklabel',[])
%         set(gca,'ytick',[])
%         set(gca,'yticklabel',[])
%
%         %________________________ FIRING RATE MAPS
%         allmapdata{i}{k} = pfm{1}(i).output{1}(k); %just mod for now
%         totalmapplots = totalplots+5;
%         subplot(4,5,rem(totalmapplots,20)+1); hold on;
%         imagesc(allmapdata{i}{k}.smoothedspikerate')
%         %                         set(gca,'YLim',[0 110]);
%         %                         set(gca,'XLim',[0 117]);
%
%         titstring = [animals{i}, '________Mod______', tit, '_____cell group:', num2str(figcnt)];
%
%
%         tmpspatind = spatmoddata_alltraj(:,[1:4]);
%         if ~rowfind(allmapdata{i}{k}.index,tmpspatind) == 0;
%             currmapspatinfo = spatmoddata_alltraj(rowfind(allmapdata{i}{k}.index,tmpspatind),6);
%             if currmapspatinfo <1;
%                 color = 'b';
%             elseif currmapspatinfo <3;
%                 color = [0 .5 0];
%             elseif currmapspatinfo <10;
%                 color = 'r';
%             elseif currmapspatinfo >=10;
%                 color = 'k';
%
%             end
%             %                             title([num2str(currmapspatinfo)], 'FontSize',24,'Fontweight','bold', 'Color', color);
%             title({['SI=',num2str(currmapspatinfo)]; ['maxF=',num2str(max(max(allmapdata{i}{k}.smoothedspikerate)))]; [sprintf('D%d E%d T%d C%d',pfm{1}(i).output{1}(k).index(1),pfm{1}(i).output{1}(k).index(2),pfm{1}(i).output{1}(k).index(3),pfm{1}(i).output{1}(k).index(4))]}, 'FontSize',24,'Fontweight','bold', 'Color', color);
%             maxratesmod(k) = max(max(allmapdata{i}{k}.smoothedspikerate));
%             spatinfoallmod(k) = currmapspatinfo;
%         else
%             title(['~exist'], 'FontSize',24,'Fontweight','bold');
%         end
%
%         %                         xlabel(typecell, 'Color', typclr)
%         set(gca,'xtick',[])
%         set(gca,'xticklabel',[])
%         set(gca,'ytick',[])
%         set(gca,'yticklabel',[])
%         % Update plot number
%         totalplots=totalplots+1;
%         if (rem(totalplots,20) >4 & rem(totalplots,20) <10) | (rem(totalplots,20) >14);
%             totalplots = totalplots+5; %skip the firing rate maps lines
%         end
%         if rem(totalmapplots+1,20)==0
%
%             %save the full maps figures
%             figfile = [anfigdir,'/','ModMaps',num2str(figcnt)];
%             [ax,hU]=suplabel(titstring,'t', [.08 .08 .84 .84]);
%             set(hU,'FontSize',200,'Fontweight','bold')  %font controls aren't working.. i think it's a bug with ubuntu matlab2013a...
%             if savefigs==1
%                 print('-djpeg', figfile);
%                 %                  print('-dpdf', figfile);
%                 %saveas(gcf,figfile,'fig');
%                 %                                 print('-depsc2', figfile);
%             end
%
%             if ~cyclemaps == 0
%                 keyboard;
%             end
%             close
%         end
%     end
%
%     %save the last maps figure as it's most likely
%     %not full
%     figfile = [anfigdir,'/','ModMaps',num2str(figcnt)];
%     [ax,hU]=suplabel(titstring,'t', [.08 .08 .84 .84]);
%     set(hU,'FontSize',200,'Fontweight','bold')  %font controls aren't working.. i think it's a bug with ubuntu matlab2013a...
%     if savefigs==1
%         print('-djpeg', figfile);
%         %                  print('-dpdf', figfile);
%         %                         saveas(gcf,figfile,'pdf')
%         %                         print('-depsc2', figfile);
%     end
%     if~cyclemaps ==0
%         keyboard
%     end
%     close
% end
