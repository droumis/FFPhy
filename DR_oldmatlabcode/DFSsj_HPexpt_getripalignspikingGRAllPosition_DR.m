
% Ripple modulation of cells, especilally PFC cells. Time filter version of sj_HPexpt_ripalign_singlecell_getrip4.
% Will call DFAsj_getripalign.m
% Also see DFSsj_plotthetamod.m and DFSsj_HPexpt_xcorrmeasures2. Will gather data like these

clear; %close all;
runscript = 1;
savedata = 0; % save data option - only works if runscript is also on
figopt1 = 0; % Figure Options - Individual cells
cyclemaps =0; %pause plots.. DR not doing this to create a struct to overlay on firing
saverippos = 'rippos_Jan26';

savedir = '/data19/sjadhav/HPExpt/ProcessedDataDR/';
savefile = [savedir 'AllAn_PFCCA1_ripplepos_DR_Jan26']; area = 'PFC'; clr = 'b'; % PFC
%savefile = [savedir 'HP_ripplemod_CA1']; area = 'CA1';  clr = 'r';% CA1
%savefile = [savedir 'HP_ripplemod_PFC_speed']; area = 'PFC'; clr = 'b'; % PFC - low speed criterion

savefig1=0;


% Plot options
plotanimidx =  []; % To pick animals for plotting
plotdays = []; % If you only load data when runscript=0 and savedata=0, then this field will supplant days


%If runscript, run Datafilter and save data
if runscript == 1
    
    %Animal selection
    %-----------------------------------------------------
%     animals = {'HPa' 'HPb' 'HPc' 'nadal'};
    animals = {'HPa' 'HPb' 'HPc' 'Ndl' 'Rtl' 'Brg'};
    %Filter creation
    %-----------------------------------------------------
    
    % Epoch filter
    % -------------
      dayfilter = ''; % Shantanu - I am adding day filter to parse out epoch filter
    % Either Only do 1st w-track. 2 or 1 epochs per day
    % Or do Wtr1 and Wtr2, 2 epochs per day
%        runepochfilter = 'isequal($environment, ''wtr1'') || isequal($environment, ''wtr2'')';
    % GIDEON: THIS INSTEAD, CHECK WITH SHANTANU:
%     runepochfilter = 'isequal($type, ''run'')';
%DR
runepochfilter = 'isequal($type, ''run'') && ~isequal($environment, ''lin'')';
    
    
    %sleepepochfilter = 'isequal($type, ''sleep'')'; % Only pre and post sleep marked as sleep
    
    % Cell filter
    % -----------
    %cellfilter = 'strcmp($area, ''PFC'')'; % This includes all, including silent cells
    %cellfilter = ' (strcmp($area, ''CA1'') || strcmp($area, ''iCA1'')) && ($numspikes > 100) && ~strcmp($tag2, ''CA1Int'') && ~strcmp($tag2, ''iCA1Int'')';
    
    % GIDEON: NO INTERNEURON EXCLUSION:
    cellfilter =  '( strcmp($area, ''CA1'') || strcmp($area, ''iCA1'') || strcmp($area, ''PFC'')) && ($numspikes > 100)'; %'strcmp($area, ''PFC'')'; 
    
    
    % cellfilter = '(strcmp($area, ''CA1'') || strcmp($area, ''iCA1'')) && ($numspikes > 100) && ($meanrate < 7)';
    
    % cellfilter = 'strcmp($tag2, ''PFC'')';
    %cellfilter = 'strcmp($tag2, ''CA1Pyr'') && ($numspikes > 100)'; % This includes all.
    % For more exclusive choosing, use $tag. Take care of number of spikes while gathering data
    
    % Time filter
    % -----------
    riptetfilter = '(isequal($descrip, ''riptet''))';
    
    % The following are ripple time filter options. Instead of using time filters, within the function,
    % call getripples using the riptetfilter option to get the riptimes. Since you want a vector of riptimes.
    % You can also use inter-ripple-interval of 1 sec within the function, etc.
    
    % a) Using getriptimes: generates 1 vector with 1's for ripple times for all tetrodes.
    % Should be similar to getripples/ getripplees_direct
    
       timefilterrun_rip = {{'DFTFsj_getvelpos', '(($absvel <= 5))'}, {'DFTFsj_getriptimes','($nripples > 1)','tetfilter',riptetfilter,'minthresh',3}};
    
    % b) getripples. This is similar to getripltimes, but gives start and end of each ripple time
    % instead of a time filter. Also has an additional condition of at least 50 ms ripple.
    
    % Iterator
    % --------
    iterator = 'singlecellanal';
    
    % Filter creation
    % ----------------
        modf = createfilter('animal',animals,'days',dayfilter,'epochs',runepochfilter, 'cells', cellfilter, 'excludetime', timefilterrun_rip, 'iterator', iterator);
%     modf = createfilter('animal',animals,'epochs',runepochfilter, 'cells',...
%         cellfilter, 'iterator', iterator);
    disp('Done Filter Creation');
    
    % Set analysis function
    % ----------------------
    
    modf = setfilterfunction(modf,'DFAsj_getripalignspiking_position',{'spikes', 'ripples', 'tetinfo', 'pos'}); %
    %modf = setfilterfunction(modf,'DFAsj_getripalignspiking',{'spikes', 'ripples', 'tetinfo', 'pos'},'dospeed',1); %
    % Going to call getripples_tetinfo within function
    
    % Run analysis
    % ------------
    modf = runfilter(modf);
    disp('Finished running filter script');
    %--------------------- Finished Filter Function Run -------------------
    
    if savedata == 1
        clear figopt1 runscript plotdays plotanimidx savedata
        save([savefile saverippos]);
    end
    
else
    
    load([savefile saverippos]);
    
end % end runscript

if ~exist('savedata')
    return
end


% -------------------------  Filter Format Done -------------------------



% ----------------------------------
% Whether to gather data or to load previously gathered data
% -------------------------------------------------------------------
gatherdata = 1; savegatherdata = 1;
gatherdatafile = [savedir 'Ndl_ripplemod_PFC_pos_gather']; % PFC cells to Hipp ripples
%gatherdatafile = [savedir 'HP_ripplemod_CA1_gather']; % CA1 cells to Hipp ripples
%gatherdatafile = [savedir 'HP_ripplemodspeed_PFC_gather']; % PFC cells to Hipp ripples - low speed criterion


if gatherdata
    
    % Parameters if any
    % -----------------
    
    % -------------------------------------------------------------
    
    cnt=0; % Count how many cells will be kept based on nspikes in output: >0
    allanimindex=[]; alldataraster=[]; alldatahist = []; all_Nspk=[];
    
    for an = 1:length(modf)
        for i=1:length(modf(an).output{1})
            % Check for empty output - If Cell defined in rpoch and Nspks in ripple response wndow > 0
            if (modf(an).output{1}(i).Nspikes > 0)
                cnt=cnt+1;
                anim_index{an}(cnt,:) = modf(an).output{1}(i).index;
                % Only indexes
                animindex=[an modf(an).output{1}(i).index]; % Put animal index in front
                allanimindex = [allanimindex; animindex]; % Collect all Anim Day Epoch Tet Cell Index
                % Data
                alldataraster{cnt} = modf(an).output{1}(i).rip_spks_cell; % Only get raster and histogram response
                alldatahist{cnt} = modf(an).output{1}(i).rip_spkshist_cell;
                all_Nspk(cnt) = modf(an).output{1}(i).Nspikes;
                alldataraster_rdm{cnt} = modf(an).output{1}(i).rdm_spks_cell; % Only get raster and histogram response
                alldatahist_rdm{cnt} = modf(an).output{1}(i).rdm_spkshist_cell;
                % trialResps: Summed Nspks/trial in respective window
                alldatatrialResps{cnt} = modf(an).output{1}(i).trialResps;
                alldatatrialResps_bck{cnt} = modf(an).output{1}(i).trialResps_bck;
                alldatatrialResps_rdm{cnt} = modf(an).output{1}(i).trialResps_rdm;
                %New
                allrippos{cnt}=modf(an).output{1}(i).ripPos;
                
                %end
                if cnt==1
                    pret =  modf(an).output{1}(i).pret;
                    postt = modf(an).output{1}(i).postt;
                    binsize = modf(an).output{1}(i).binsize;
                    rwin = modf(an).output{1}(i).rwin;
                    bckwin = modf(an).output{1}(i).bckwin;
                    bins_resp  = modf(an).output{1}(i).bins_resp;
                    bins_bck = modf(an).output{1}(i).bins_bck;
                    timeaxis = modf(an).output{1}(i).timeaxis;
                end
            end
        end
        
    end
    
    % Consolidate single cells across epochs. Multiple methods: see also DFSsj_getcellinfo and DFSsj_xcorrmeasures2
    % ----------------------------------------------------------------------------
    
    allripplemod = struct;
    
    % Method 1
    % ---------------------------------------------
    %     cntcells=0;
    %     animdaytetcell = unique(allanimindex(:,[1 2 4 5]),'rows'); % Collapse across epochs
    %     for ind = 1:size(animdaytetcell,1)
    %         currhist=[]; currraster = [];
    %         for a = 1:size(allanimindex,1)
    %             if animdaytetcell(ind,:)==allanimindex(a,[1 2 4 5]) % Epoch match
    %                 currhist = [currhist; alldatahist{a}];
    %                 currraster = [currraster; alldataraster{a}];
    %             end
    %         end
    %         cntcells = cntcells + 1;
    %         allripplemod_idx(cntcells,:)=animdaytetcell(ind,:);
    %         allripplemod(cntcells).index=animdaytetcell(ind,:);
    %         allripplemod(cntcells).hist=currhist;
    %         allripplemod(cntcells).raster=currraster;
    %         allripplemod(cntcells).Nspk=Nspk;
    %     end
    
    
    % Method 2
    % ---------
    dummyindex=allanimindex;  % all anim-day-epoch-tet-cell indices
    cntcells=0;
    % all unique indices of anim-day-tet-cell, collapsing across epochs
    uniqueIndices=unique(allanimindex(:,[1,2,4,5]),'rows');
    % iterating only over the unique indices and finding matches in
    % allanimindex
    for i=1:length(uniqueIndices)
        animdaytetcell=uniqueIndices(i,:);
        %         ind=[];
        %         while rowfind(animdaytetcell,dummyindex(:,[1 2 4 5]))~=0          % collect all rows (epochs)
        %             ind = [ind rowfind(animdaytetcell,dummyindex(:,[1 2 4 5]))];        % finds the first matching row
        %             dummyindex(rowfind(animdaytetcell,dummyindex(:,[1 2 4 5])),:)=[0 0 0 0 0]; % after adding index, remove the corresponding row
        %             % so you could find the next one
        %         end
        
        %Gideon: simpler way of finding indices of same cell across epochs:
        ind=find(ismember(allanimindex(:,[1 2 4 5]),animdaytetcell,'rows'))';
        
        
        % Gather everything for the current cell across epochs
        currhist=[]; currraster=[]; currNspk=0;
        currhist_rdm=[]; currraster_rdm=[];
        currtrialResps=[]; currtrialResps_rdm=[]; currtrialResps_bck=[];
        currippos=[];
        for r=ind
            currNspk = currNspk + all_Nspk(r);
            currhist = [currhist; alldatahist{r}];
            currraster = [currraster, alldataraster{r}];
            currhist_rdm = [currhist_rdm; alldatahist_rdm{r}];
            currraster_rdm = [currraster_rdm, alldataraster_rdm{r}];
            currtrialResps = [currtrialResps, alldatatrialResps{r}];
            currtrialResps_rdm = [currtrialResps_rdm, alldatatrialResps_rdm{r}];
            currtrialResps_bck = [currtrialResps_bck, alldatatrialResps_bck{r}];
            currippos=[currippos;allrippos{r}];
        end
        if currNspk >= 50
            cntcells = cntcells + 1;
            allripplemod_idx(cntcells,:)=animdaytetcell;
            allripplemod(cntcells).index=animdaytetcell;
            allripplemod(cntcells).hist=currhist*(1000/binsize); % Convert to firing rate in Hz
            allripplemod(cntcells).raster=currraster;
            allripplemod(cntcells).Nspk=currNspk;
            allripplemod(cntcells).hist_rdm=currhist_rdm*(1000/binsize); % Convert to firing rate in Hz
            allripplemod(cntcells).raster_rdm=currraster_rdm;
            % Trial Resps
            allripplemod(cntcells).trialResps = currtrialResps';
            allripplemod(cntcells).trialResps_rdm = currtrialResps_rdm';
            allripplemod(cntcells).trialResps_bck = currtrialResps_bck';
            allripplemod(cntcells).rippos=currippos;
        end
    end
    
    % GIDEON
    respRange=[-100:200]+550;
    
    for ii=1:cntcells
        ii;
        curRast=allripplemod(ii).raster;
        numRips=length(curRast);
        curRastMat=zeros(numRips,1101);
        curRastMatS=zeros(numRips,1101);
        curRastMatF=zeros(numRips,1101);
        
        for i=1:numRips
            curRastMat(i,round(curRast{i}+551))=1;
            curRastMatS(i,:)=smooth(curRastMat(i,:),50);
            curRastMatF(i,:)=(curRastMatS(i,:)-mean(curRastMatS(i,:)));
            
        end
        %figure;
        meanResp=mean(curRastMatS);
        meanRespF=meanResp-mean(meanResp);
        
        meanRespFwin=meanRespF(respRange);
        curRastMatFwin=curRastMatF(:,respRange);
        
        respAmps=curRastMatFwin*meanRespFwin';
        allripplemod(ii).respAmps=respAmps;
        
    end
    %--------------------------------------------------
    indexALLlist = [];%DR
    clear ripplepos;
    numSTDs=0;
    for i=1:cntcells;
        indexALLlist = [indexALLlist; allripplemod(i).index];
        c1=allripplemod(i).rippos;
        if mean(allripplemod(i).trialResps)>mean(allripplemod(i).trialResps_bck);
            highresps2=find(allripplemod(i).trialResps>(mean(allripplemod(i).trialResps)+numSTDs*std((allripplemod(i).trialResps))));
        else
            highresps2=find(allripplemod(i).trialResps<(mean(allripplemod(i).trialResps)-1*std((allripplemod(i).trialResps))));
        end
        if ~isempty(highresps2)
            ripplepos{i} = c1(highresps2,:);
        else
            ripplepos{i} = [];
        end
    end
%     save(sprintf('%s_%s', savefile, saverippos), 'indexALLlist', 'ripplepos','allripplemod')
    save(sprintf('%s_%s', savedir, saverippos), 'indexALLlist', 'ripplepos','allripplemod')
    'saved'
   %--------------------------------------------------- 
%     numSTDs=2;
%     for i=1:cntcells
%         curridx = allripplemod(i).index;
%         day = curridx(2); tet = curridx(3); cell = curridx(4);
%         c1=allripplemod(i).rippos;
%         figure;
%          set(gcf, 'Position',[205 136 723 946]);
%         subplot(2,1,1)
%         plot(c1(:,1),c1(:,2),'bx')
%         highresps=find(allripplemod(i).respAmps>(mean(allripplemod(i).respAmps)+numSTDs*std((allripplemod(i).respAmps))))
%         hold on
%         if ~isempty(highresps)
%             plot(c1(highresps,1),c1(highresps,2),'ro','markerfacecolor','r')
%         end
%         subplot(2,1,2)
%         plot(c1(:,1),c1(:,2),'bx')
%         
%         if mean(allripplemod(i).trialResps)>mean(allripplemod(i).trialResps_bck)
%             highresps2=find(allripplemod(i).trialResps>(mean(allripplemod(i).trialResps)+numSTDs*std((allripplemod(i).trialResps))))
%         else
%             highresps2=find(allripplemod(i).trialResps<(mean(allripplemod(i).trialResps)-1*std((allripplemod(i).trialResps))))
%         end
%         
%         hold on
%         if ~isempty(highresps2)
%             plot(c1(highresps2,1),c1(highresps2,2),'ro','markerfacecolor','r')
%         end
%         title(['Day: ' num2str(day) ' tet: ' num2str(tet) ' cell: ' num2str(cell)],'FontSize',40,'Fontweight','normal');
% 
%         if cyclemaps == 1;
%             keyboard
%         end
%     end
    % Calculations/ Stats. Stats between response and b
    % Similar to ...getrip4
    % -----------------------------------------------------------
    
    % Save
    % -----
    if savegatherdata == 1
        save(gatherdatafile);
    end
    
else % gatherdata=0
    
    load(gatherdatafile);
    
end % end gather data


