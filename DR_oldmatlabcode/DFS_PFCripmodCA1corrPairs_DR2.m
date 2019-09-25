
warning('off','all');
clear; close all;

runscript = 0;
savedata = 0; % save data option - only works if runscript is also on
savefigs=0;
plotstuff = 1; %
plotEachcell = 1;
cyclemaps = 1;
peakthresh = 3;
% runnadal =0;

plotCorrCoef =1;
runPFCCA1maps =1;
plotPFCCA1maps =0;
savedir = '/mnt/data25/sjadhav/HPExpt/ProcessedDataDR/';
savefilename = 'AllAn_sigcorrPair_CA1structs_eachAn';
savefile = [savedir savefilename]; area = 'PFC'; %clrunmod = 'r'; clrmod = 'b'; % PFC
figdir = '/mnt/data25/sjadhav/HPExpt/Figures_DR/';
combineHPNdl = 1;
Veqn = '>3';
minV=str2num(Veqn(end));
mintime = 3;
traj = [1:4] ;

% If runscript, run Datafilter and save data
if runscript == 1
    %     for i =  modUnmod;
    %Animal selection
    %-----------------------------------------------------
    animals = {'HPa' 'HPb' 'HPc' 'nadal'};
    %             animals = {'HPa' 'HPb' 'HPc'};
    %         animals = {'HPc'};
    %     animals = {'HPa'};
    %         animals = {'nadal'};
    
    %Filter creation
    %-----------------------------------------------------
    % Epoch filter
    % -------------
    %     if runnadal == 1;
    %         dayfilter = '8:17';
    %     else
    dayfilter = '';
    %     end
    
    %     if runnadal == 1;
    runepochfilter = 'isequal($type, ''run'') && ~isequal($environment, ''lin'')';
    %     else
    %         runepochfilter = 'isequal($environment, ''wtr1'') || isequal($environment, ''wtr2'')';
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
    %     psf = setfilterfunction(spatf, 'DFAsj_filtercalclinfields_tf',{'spikes', 'linpos'}, 'binsize', 2);
    %     pmf = setfilterfunction(spatf, 'DFAsj_openfieldrate_tf',{'spikes', 'linpos', 'pos'}, 'binsize', 1, 'std', 2);
    fields = setfilterfunction(spatf, 'DFA_loadfields_DR', {'linfields', 'mapfields', 'cellinfo'});
    
    % Run analysis-----------------------
    %     flds = runfilter(psf);  % Place Field Stability.. trajectories
    %     pfm = runfilter(pmf);  % Place Field Map
    flds = runfilter(fields);  % Place Field Map
    
    %     end
    disp('Finished running filter script');
    %--------------------- Finished Filter Function Run -------------------
    
    %Create CA1 struct for ALL HP animals_____________________________________________
    ca1cnt = 0; ca1ndl = 0; ca1hpa =0; ca1hpb = 0; ca1hpc = 0;
    for ani = 1:3; % loop thru HP anims
        for ai = 1:length(flds(1,ani).output{1}); %loop over cells in anim
            if strcmp(flds(1,ani).output{1}(1,ai).area,'iCA1') || strcmp(flds(1,ani).output{1}(1,ai).area,'CA1'); %if ca1
                if length(flds(1,ani).output{1}(1,ai).trajdata) == 4; %if it has all 4 trajs
                    ca1cnt = ca1cnt +1;
                    HPCA1struct{ca1cnt,1}= flds(1,ani).output{1}(1,ai); %store ca1 data
                end
            end
        end
    end
    
    %create ca1 struct for nadal
    for ani = 4; % loop thru anims
        for ai = 1:length(flds(1,ani).output{1}); %loop over cells in anim
            if strcmp(flds(1,ani).output{1}(1,ai).area,'iCA1') || strcmp(flds(1,ani).output{1}(1,ai).area,'CA1');
                if length(flds(1,ani).output{1}(1,ai).trajdata) == 4; %if it has all 4 trajs
                    ca1ndl = ca1ndl +1;
                    NdlCA1struct{ca1ndl,1}= flds(1,ani).output{1}(1,ai); %store ca1 data
                end
            end
        end
    end
    
    %create ca1 struct for HPa
    for ani = 1; % loop thru anims
        for ai = 1:length(flds(1,ani).output{1}); %loop over cells in anim
            if strcmp(flds(1,ani).output{1}(1,ai).area,'iCA1') || strcmp(flds(1,ani).output{1}(1,ai).area,'CA1');
                if length(flds(1,ani).output{1}(1,ai).trajdata) == 4; %if it has all 4 trajs
                    ca1hpa = ca1hpa +1;
                    HPaCA1struct{ca1hpa,1}= flds(1,ani).output{1}(1,ai); %store ca1 data
                end
            end
        end
    end
    
    %create ca1 struct for HPb
    for ani = 2; % loop thru anims
        for ai = 1:length(flds(1,ani).output{1}); %loop over cells in anim
            if strcmp(flds(1,ani).output{1}(1,ai).area,'iCA1') || strcmp(flds(1,ani).output{1}(1,ai).area,'CA1');
                if length(flds(1,ani).output{1}(1,ai).trajdata) == 4; %if it has all 4 trajs
                    ca1hpb = ca1hpb +1;
                    HPbCA1struct{ca1hpb,1}= flds(1,ani).output{1}(1,ai); %store ca1 data
                end
            end
        end
    end
    
    %create ca1 struct for HPc
    for ani = 3; % loop thru anims
        for ai = 1:length(flds(1,ani).output{1}); %loop over cells in anim
            if strcmp(flds(1,ani).output{1}(1,ai).area,'iCA1') || strcmp(flds(1,ani).output{1}(1,ai).area,'CA1');
                if length(flds(1,ani).output{1}(1,ai).trajdata) == 4; %if it has all 4 trajs
                    ca1hpc = ca1hpc +1;
                    HPcCA1struct{ca1hpc,1}= flds(1,ani).output{1}(1,ai); %store ca1 data
                end
            end
        end
    end
    
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

if savefigs == 1;
    mkdir(figdir,savefilename)
end

%prepping data structs
%________________________________________________________________________
if runPFCCA1maps ==1;
    pfclistidx = []; ju = 0; fldsLIST = []; pu =0;  cnt = 0; trajdataall =[]; CorrCoefstruct = []; clear fldsdata.output fldsLIST fldsY fldsZ fldsShortList mfldsdata.output mfldsLIST mfldsY mfldsZ mfldsShortList
    clr = {'b','r','g','m','c','y','k','r'};
    load('/mnt/data25/sjadhav/HPExpt/ProcessedData/HP_allPFCCA1sigcorridxs'); %combined epochs
    %     load('/mnt/data25/sjadhav/HPExpt/ProcessedData/HP_allPFCCA1sigidxs'); %individual epochs... doubles the cells
    
    pairdata = allPFCCA1sigcorridxs;
    for yu = 1:length(pairdata);
        PFClistindex(yu,:) = pairdata(1,yu).PFCidx(:);
        PFCY = sprintf('%d', PFClistindex(yu,:)); %make string out of all columns an day ep tet cell
        PFCZ = str2num(PFCY); %convert string to num
        PFCShortList(yu,1) =  PFCZ; %save for later
    end
    
    %transform the flds traj data into a workable and searchable format
    for anims = 1:length(animals);
        for bu = 1:length(flds(anims).output{1});
            cnt = cnt+1;
            fldsdata.output{cnt} = flds(anims).output{1}(1,bu); %collect all data across animals
            fldsLIST(cnt,:) = [anims , fldsdata.output{cnt}.index([1 3 4])]; %collect indices W/O EP and add animal num in first column to match pairdata format
            fldsY = sprintf('%d', fldsLIST(cnt,:)); %make string out of all columns an day ep tet cell
            fldsZ = str2num(fldsY); %convert string to num
            fldsShortList(cnt,1) =  fldsZ; %save for later
        end
    end
    
    %transform the flds map data into a workable and searchable format
    cnt = 0;
    for anims = 1:length(animals);
        for bu = 1:length(flds(anims).output{1});
            cnt = cnt+1;
            mfldsdata.output{cnt} = flds(anims).output{1}(1,bu); %collect all data across animals
            mfldsLIST(cnt,:) = [anims, mfldsdata.output{cnt}.index([1 3 4 ])]; %collect indices W/O EP and add animal num in first column to match pairdata format
            mfldsY = sprintf('%d', mfldsLIST(cnt,:)); %make string out of all columns an day ep tet cell
            mfldsZ = str2num(mfldsY); %convert string to num
            mfldsShortList(cnt,1) =  mfldsZ; %save for later
        end
    end
    
    %____________________________________________________________________________________________________________________
    %loop through the PFCShortList. for each i pfc, find PFCShortList(i) in
    %fldsShortList then take that row num and plot the traj from
    %fldsdata.output and the map from pfmdata.output.. Then for the length of
    %the CA1sigind of the pairdata(1,i), create short list of CA1 indices, then loop through each row to find match in fldsShortList then take that row num and store the traj from
    %fldsdata.output and the map from pfmdata.output
    
    
    paircounter = 0; pospaircounter = 0; negpaircounter = 0; mocnt = 0; hpacnt = 0; hpbcnt = 0; hpccnt = 0; ndlcnt = 0;
    for i = 1:length(PFCShortList);
        clear pfcmatchtraj pfcmatchmap trajplot mapplot CA1X  CA1Y CA1Z CA1ShortList CA1matchtraj CA1matchmap CA1bX CA1bY CA1bZ CA1ShortListb CA1matchtrajb CA1matchmapb
        pfcmatchtraj = find(PFCShortList(i) == fldsShortList(:), 1, 'last'); %find corresponding index from psfdata.need to do last bc day 1 has lin track first... update.. now i'm filtering out lin track so that i should only have w track run epochs.. but im still only using the last ep.. 
        pfcmatchmap = find(PFCShortList(i) == mfldsShortList(:), 1, 'last'); %find corresponding index from pfmdata
        alltraj = [];
        for tu = 1:length(fldsdata.output{pfcmatchtraj}.trajdata);
            alltraj = [alltraj; fldsdata.output{pfcmatchtraj}.trajdata{tu}(:,5)];
        end
        if max(alltraj) > peakthresh; %use cell if 5Hz peak threshold. if not, this whole pfc-ca1 set will be skipped.
            trajplot{1,1} = fldsdata.output{pfcmatchtraj}.trajdata; %store pfc traj data
            mapplot{1,1} = mfldsdata.output{pfcmatchmap}.mapdata.smoothedspikerate; %store pfc map data.
            
            
            
            %loop over the list of CA1 cells for the current pfc cell and make a
            %short list. then loop through the shortlist, finding traj, maps,
            %and then store or plot data for each.
            %updated.. if pos idx isn't empty, do all this for those
            %then, if neg idx isn't empty, do all this for those concatenating
            %the pos, if any
            %update2: loop using sigidxs but add to data pos, then neg
            
            %         for cai = 1:length(pairdata(1,i).CA1sigidxs{1}); %length of pos +neg
            %pos CA1 data
            if ~isempty(pairdata(1,i).CA1posidxs);
                for posi = 1:length(pairdata(1,i).CA1posidxs(:,1)); %is this the correct tag?
                    CA1X = pairdata(1,i).CA1posidxs(posi,:); %is this the correct tag?
                    CA1Y = sprintf('%d', CA1X);
                    CA1Z = str2num(CA1Y);
                    CA1ShortList(posi,1) = CA1Z;
                    CA1matchtraj = find(CA1ShortList(posi) == fldsShortList(:), 1, 'last'); %find last corresponding index from psfdata. .need to do last bc day 1 has lin track first
                    CA1matchmap = find(CA1ShortList(posi) == mfldsShortList(:), 1, 'last'); %find corresponding index from pfmdata
                    alltraj = [];
                    for tu = 1:length(fldsdata.output{CA1matchtraj}.trajdata);
                        alltraj = [alltraj; fldsdata.output{CA1matchtraj}.trajdata{tu}(:,5)];
                    end
                    if max(alltraj) > peakthresh; %use cell if 5Hz peak threshold
                        trajplot{1, length(trajplot)+1} = fldsdata.output{CA1matchtraj}.trajdata; %store CA1 traj data next to pfc data
                        mapplot{1,length(mapplot)+1} = mfldsdata.output{CA1matchmap}.mapdata.smoothedspikerate; %store CA1 map data next to pfc data
                    end
                end %CA1 pos short list and store
            else
                CA1ShortList =[];
            end %if pos not empty
            
            %neg CA1 data
            if ~isempty(pairdata(1,i).CA1negidxs);
                for negi = 1:length(pairdata(1,i).CA1negidxs(:,1));
                    CA1bX = pairdata(1,i).CA1negidxs(negi,:);
                    CA1bY = sprintf('%d', CA1bX);
                    CA1bZ = str2num(CA1bY);
                    CA1ShortListb(negi,1) = CA1bZ;
                    CA1matchtrajb = find(CA1ShortListb(negi) == fldsShortList(:), 1, 'last'); %find corresponding index from psfdata.need to do last bc day 1 has lin track first
                    CA1matchmapb = find(CA1ShortListb(negi) == mfldsShortList(:), 1, 'last'); %find corresponding index from pfmdata
                    alltraj = [];
                    for tu = 1:length(fldsdata.output{CA1matchtrajb}.trajdata);
                        alltraj = [alltraj; fldsdata.output{CA1matchtrajb}.trajdata{tu}(:,5)];
                    end
                    if max(alltraj) > peakthresh; %use cell if 5Hz peak threshold
                        trajplot{1, length(trajplot)+1} = fldsdata.output{CA1matchtrajb}.trajdata; %store CA1 traj data next to pfc data
                        mapplot{1,length(mapplot)+1} = mfldsdata.output{CA1matchmapb}.mapdata.smoothedspikerate; %store CA1 map data next to pfc data
                    end
                end %CA1 neg short list and store
            else
                CA1ShortListb =[];
            end %if pos not empty
            
            %store all traj data
            trajdataall{length(trajdataall)+1,1} = trajplot; %store each pfc set for corrcoef
            
            %compute corr coef for current set.. also creating 3
            %different rand shuffled structs: all, pos, neg.. matching the
            %length of the corresponding non-shuff struct.
            if length(trajplot) >1; %if there are any ca1 cells
                for cu = 2:(length(trajplot)); %all ca1 cells ~ # of pairs w pfc
                    if length(trajplot{1}) ==4 && length(trajplot{cu}) ==4; %the pfc cell and ca1 cell must have all trajectories
                        paircounter = paircounter +1; %count each pair for storing
                        %find the index of a random ca1 cell,,update: i used
                        %this when only using 1 rand cell per pair for
                        %shuff struct..
                        %                             if pairdata(1,i).PFCidx(1,1) < 4; %HP
                        %                                 randca1 = HPCA1struct{randi(length(HPCA1struct)),1}; %take a random ca1 cell from one of the HP animals
                        %                             else %Ndl
                        %                                 randca1 = NdlCA1struct{randi(length(NdlCA1struct)),1}; %take a random ca1 cell from nadal
                        %                             end
                        for tru = 1:length(trajplot{1}); %for each trajectory
                            clear r p minlen;
                            minlen = min(length(trajplot{1}{tru}(:,5)), length(trajplot{cu}{tru}(:,5)));
                            [r, p] = corrcoef(trajplot{1}{tru}(1:minlen,5), trajplot{cu}{tru}(1:minlen,5), 'rows', 'pairwise'); %don't use rows if either has a nan
                            CorrCoefstruct{paircounter,1}.r(tru,1) = r(2,1); %store corr coefs
                            CorrCoefstruct{paircounter,1}.p(tru,1) = p(2,1); %store corr coefs p vals
                            %
                            %                                 for tru = 1:length(trajplot{1}); %for each trajectory
                            %                                     clear r p minlen;
                            %                                     minlen = min(length(trajplot{1}{tru}(:,5)), length(randca1.trajdata{tru}(:,5)));
                            %                                     [r, p] = corrcoef(trajplot{1}{tru}(1:minlen,5), randca1.trajdata{tru}(1:minlen,5), 'rows', 'pairwise'); %don't use rows if either has a nan
                            %                                     CorrCoefstructShuff{paircounter,1}.r(tru,1) = r(2,1); %store corr coefs
                            %                                     CorrCoefstructShuff{paircounter,1}.p(tru,1) = p(2,1); %store corr coefs p vals
                            %                                 end
                        end %each traj
                        
                        %run all that for the random cell.. update:
                        %moving the shuff struct outside of the traj
                        %loop so that i can loop thru all ca1 cells for
                        %each pfc cell. update2: moving the shuffled
                        %struct out of the pair loop entirely as it's
                        %now dependent on # of pfc cells.
                        %
                        %                             for mo = 1:length(HPCA1struct); %loop over all the HP CA1 cells for the shuffled struct
                        %                                 mocnt = mocnt +1;
                        %                                 for tru = 1:length(trajplot{1}); %for each trajectory
                        %                                     clear r p minlen;
                        %                                     minlen = min(length(trajplot{1}{tru}(:,5)), length(HPCA1struct{mo}.trajdata{tru}(:,5)));
                        %                                     [r, p] = corrcoef(trajplot{1}{tru}(1:minlen,5), HPCA1struct{mo}.trajdata{tru}(1:minlen,5), 'rows', 'pairwise'); %don't use rows if either has a nan
                        %                                     CorrCoefstructShuff{mocnt,1}.r(tru,1) = r(2,1); %store corr coefs  %update using a different counter bc now they'll be exp2 as more
                        %                                     CorrCoefstructShuff{mocnt,1}.p(tru,1) = p(2,1); %store corr coefs p vals
                        %                                 end
                        %                             end
                        
                        if cu-1 <= length(pairdata(1,i).CA1posidxs(:,1)); %if pos corr ca1 cell
                            CorrCoefstruct{paircounter, 1}.posneg = 'pos';
                            CorrCoefstruct{paircounter, 1}.pfcindex(1,:) = pairdata(1,i).PFCidx(1,:);
                            CorrCoefstruct{paircounter, 1}.ca1index(1,:) = pairdata(1,i).CA1sigidxs(cu-1,:);
                            pospaircounter = pospaircounter +1;
                            CorrCoefstructPOS{pospaircounter, 1} = CorrCoefstruct{paircounter, 1}; %make copy of pos struct for pos store
                            
                            %make struct of pos shuffled
                            %                                 CorrCoefstructShuff{paircounter, 1}.posneg = 'pos';
                            %                                 CorrCoefstructShuff{paircounter, 1}.pfcindex(1,:) = pairdata(1,i).PFCidx(1,:);
                            %                                 CorrCoefstructShuff{paircounter, 1}.ca1index(1,:) = randca1.index;
                            %                                 CorrCoefstructPOSShuff{pospaircounter, 1} = CorrCoefstructShuff{paircounter, 1}; %make copy of pos struct for pos store
                            
                        else %if neg corr ca1 cell, add -1 to (5,1), (5,2)
                            CorrCoefstruct{paircounter, 1}.posneg = 'neg';
                            CorrCoefstruct{paircounter, 1}.pfcindex(1,:) = pairdata(1,i).PFCidx(1,:);
                            CorrCoefstruct{paircounter, 1}.ca1index(1,:) = pairdata(1,i).CA1sigidxs(cu-1,:);
                            negpaircounter = negpaircounter +1;
                            CorrCoefstructNEG{negpaircounter, 1} = CorrCoefstruct{paircounter, 1}; %make copy of neg struct for neg store
                            
                            %                                 %make struct of neg shuffled
                            %                                 CorrCoefstructShuff{paircounter, 1}.posneg = 'neg';
                            %                                 CorrCoefstructShuff{paircounter, 1}.pfcindex(1,:) = pairdata(1,i).PFCidx(1,:);
                            %                                 CorrCoefstructShuff{paircounter, 1}.ca1index(1,:) = randca1.index;
                            %                                 CorrCoefstructNEGShuff{negpaircounter, 1} = CorrCoefstructShuff{paircounter, 1}; %make copy of neg struct for neg store
                        end
                    end %if cells have all trajs
                end %eac ca1 cell in set
            end %if any ca1 cells
            
            %                 %make massive shuff struct__________________
            %                 %for each 1 (pfc in short list)
            %                 if length(trajplot) >1; %if there are any ca1 cells
            %                     if pairdata(1,i).PFCidx(1,1) < 4%HPanimal
            %                         for mo = 1:length(HPCA1struct); %loop over all the HP CA1 cells for the shuffled struct
            %                             %loop over CA1 cells that are no significantly correlated with the current PFC cell
            %                             if isempty(find(str2num(sprintf('%d', HPCA1struct{mo,1}.index)) == CA1ShortList)); %if this ca1 cell is not on the pos ca1 list
            %                                 if isempty(find(str2num(sprintf('%d', HPCA1struct{mo,1}.index)) == CA1ShortListb)) ; %and if this ca1 cell is not on the neg ca1 list
            %                                     mocnt = mocnt +1;
            %                                     for tru = 1:length(trajplot{1}); %for each trajectory
            %                                         clear r p minlen;
            %                                         minlen = min(length(trajplot{1}{tru}(:,5)), length(HPCA1struct{mo}.trajdata{tru}(:,5)));
            %                                         [r, p] = corrcoef(trajplot{1}{tru}(1:minlen,5), HPCA1struct{mo}.trajdata{tru}(1:minlen,5), 'rows', 'pairwise'); %don't use rows if either has a nan
            %                                         CorrCoefstructShuff{mocnt,1}.r(tru,1) = r(2,1); %store corr coefs
            %                                         CorrCoefstructShuff{mocnt,1}.p(tru,1) = p(2,1); %store corr coefs p vals
            %                                     end
            %                                     CorrCoefstructShuff{mocnt, 1}.pfcindex(1,:) = pairdata(1,i).PFCidx(1,:);
            %                                     CorrCoefstructShuff{mocnt, 1}.ca1index(1,:) = HPCA1struct{mo}.index;
            %                                 end
            %                             end
            %                         end
            %                     else %nadal
            %                         for mo = 1:length(NdlCA1struct); %loop over all the HP CA1 cells for the shuffled struct
            %                             if isempty(find(str2num(sprintf('%d', NdlCA1struct{mo,1}.index)) == CA1ShortList));
            %                                 if isempty(find(str2num(sprintf('%d', NdlCA1struct{mo,1}.index)) == CA1ShortListb)); %loop over CA1 cells that are no significantly correlated with the current PFC cell
            %                                     mocnt = mocnt +1;
            %                                     for tru = 1:length(trajplot{1}); %for each trajectory
            %                                         clear r p minlen;
            %                                         minlen = min(length(trajplot{1}{tru}(:,5)), length(NdlCA1struct{mo}.trajdata{tru}(:,5)));
            %                                         [r, p] = corrcoef(trajplot{1}{tru}(1:minlen,5), NdlCA1struct{mo}.trajdata{tru}(1:minlen,5), 'rows', 'pairwise'); %don't use rows if either has a nan
            %                                         CorrCoefstructShuff{mocnt,1}.r(tru,1) = r(2,1); %store corr coefs
            %                                         CorrCoefstructShuff{mocnt,1}.p(tru,1) = p(2,1); %store corr coefs p vals
            %                                     end
            %                                     CorrCoefstructShuff{mocnt, 1}.pfcindex(1,:) = pairdata(1,i).PFCidx(1,:);
            %                                     CorrCoefstructShuff{mocnt, 1}.ca1index(1,:) = NdlCA1struct{mo}.index;
            %                                 end
            %                             end
            %                         end
            %                     end
            %                 end
            %
            
            % %make massive shuff struct for EACH animal as well as a huge one that
            % loops within each an but combines it to the full count.. this is
            % different than the structure above with corrcoefs each pfc cell with all
            % ca1 cells from across hp animals..
            %for each 1 (pfc in short list)
            if length(trajplot) >1; %if there are any ca1 cells
                if pairdata(1,i).PFCidx(1,1) == 1%HPa
                    for mo = 1:length(HPaCA1struct); %loop over all the HP CA1 cells for the shuffled struct
                        %loop over CA1 cells that are no significantly correlated with the current PFC cell
                        if isempty(find(str2num(sprintf('%d', HPaCA1struct{mo,1}.index)) == CA1ShortList)); %if this ca1 cell is not on the pos ca1 list
                            if isempty(find(str2num(sprintf('%d', HPaCA1struct{mo,1}.index)) == CA1ShortListb)) ; %and if this ca1 cell is not on the neg ca1 list
                                mocnt = mocnt +1; %across all anims
                                hpacnt = hpacnt +1; %within an
                                for tru = 1:length(trajplot{1}); %for each trajectory
                                    clear r p minlen;
                                    minlen = min(length(trajplot{1}{tru}(:,5)), length(HPaCA1struct{mo}.trajdata{tru}(:,5)));
                                    [r, p] = corrcoef(trajplot{1}{tru}(1:minlen,5), HPaCA1struct{mo}.trajdata{tru}(1:minlen,5), 'rows', 'pairwise'); %don't use rows if either has a nan
                                    CorrCoefstructShuff{mocnt,1}.r(tru,1) = r(2,1); %store corr coefs
                                    CorrCoefstructShuff{mocnt,1}.p(tru,1) = p(2,1); %store corr coefs p vals
                                    CorrCoefstructShuffHPa{hpacnt,1}.r(tru,1) = r(2,1); %store corr coefs
                                    CorrCoefstructShuffHPa{hpacnt,1}.p(tru,1) = p(2,1); %store corr coefs p vals
                                end
                                CorrCoefstructShuff{mocnt, 1}.pfcindex(1,:) = pairdata(1,i).PFCidx(1,:);
                                CorrCoefstructShuff{mocnt, 1}.ca1index(1,:) = HPaCA1struct{mo}.index;
                                CorrCoefstructShuffHPa{hpacnt, 1}.pfcindex(1,:) = pairdata(1,i).PFCidx(1,:);
                                CorrCoefstructShuffHPa{hpacnt, 1}.ca1index(1,:) = HPaCA1struct{mo}.index;
                            end
                        end
                    end
                    
                elseif pairdata(1,i).PFCidx(1,1) == 2%HPb
                    for mo = 1:length(HPbCA1struct); %loop over all the HP CA1 cells for the shuffled struct
                        %loop over CA1 cells that are no significantly correlated with the current PFC cell
                        if isempty(find(str2num(sprintf('%d', HPbCA1struct{mo,1}.index)) == CA1ShortList)); %if this ca1 cell is not on the pos ca1 list
                            if isempty(find(str2num(sprintf('%d', HPbCA1struct{mo,1}.index)) == CA1ShortListb)) ; %and if this ca1 cell is not on the neg ca1 list
                                mocnt = mocnt +1; %across all anims
                                hpbcnt = hpbcnt +1; %within an
                                for tru = 1:length(trajplot{1}); %for each trajectory
                                    clear r p minlen;
                                    minlen = min(length(trajplot{1}{tru}(:,5)), length(HPbCA1struct{mo}.trajdata{tru}(:,5)));
                                    [r, p] = corrcoef(trajplot{1}{tru}(1:minlen,5), HPbCA1struct{mo}.trajdata{tru}(1:minlen,5), 'rows', 'pairwise'); %don't use rows if either has a nan
                                    CorrCoefstructShuff{mocnt,1}.r(tru,1) = r(2,1); %store corr coefs
                                    CorrCoefstructShuff{mocnt,1}.p(tru,1) = p(2,1); %store corr coefs p vals
                                    CorrCoefstructShuffHPb{hpbcnt,1}.r(tru,1) = r(2,1); %store corr coefs
                                    CorrCoefstructShuffHPb{hpbcnt,1}.p(tru,1) = p(2,1); %store corr coefs p vals
                                end
                                CorrCoefstructShuff{mocnt, 1}.pfcindex(1,:) = pairdata(1,i).PFCidx(1,:);
                                CorrCoefstructShuff{mocnt, 1}.ca1index(1,:) = HPbCA1struct{mo}.index;
                                CorrCoefstructShuffHPb{hpbcnt, 1}.pfcindex(1,:) = pairdata(1,i).PFCidx(1,:);
                                CorrCoefstructShuffHPb{hpbcnt, 1}.ca1index(1,:) = HPbCA1struct{mo}.index;
                            end
                        end
                    end
                    
                elseif pairdata(1,i).PFCidx(1,1) == 3%HPc
                    for mo = 1:length(HPcCA1struct); %loop over all the HP CA1 cells for the shuffled struct
                        %loop over CA1 cells that are no significantly correlated with the current PFC cell
                        if isempty(find(str2num(sprintf('%d', HPcCA1struct{mo,1}.index)) == CA1ShortList)); %if this ca1 cell is not on the pos ca1 list
                            if isempty(find(str2num(sprintf('%d', HPcCA1struct{mo,1}.index)) == CA1ShortListb)) ; %and if this ca1 cell is not on the neg ca1 list
                                mocnt = mocnt +1; %across all anims
                                hpccnt = hpccnt +1; %within an
                                for tru = 1:length(trajplot{1}); %for each trajectory
                                    clear r p minlen;
                                    minlen = min(length(trajplot{1}{tru}(:,5)), length(HPcCA1struct{mo}.trajdata{tru}(:,5)));
                                    [r, p] = corrcoef(trajplot{1}{tru}(1:minlen,5), HPcCA1struct{mo}.trajdata{tru}(1:minlen,5), 'rows', 'pairwise'); %don't use rows if either has a nan
                                    CorrCoefstructShuff{mocnt,1}.r(tru,1) = r(2,1); %store corr coefs
                                    CorrCoefstructShuff{mocnt,1}.p(tru,1) = p(2,1); %store corr coefs p vals
                                    CorrCoefstructShuffHPc{hpccnt,1}.r(tru,1) = r(2,1); %store corr coefs
                                    CorrCoefstructShuffHPc{hpccnt,1}.p(tru,1) = p(2,1); %store corr coefs p vals
                                end
                                CorrCoefstructShuff{mocnt, 1}.pfcindex(1,:) = pairdata(1,i).PFCidx(1,:);
                                CorrCoefstructShuff{mocnt, 1}.ca1index(1,:) = HPcCA1struct{mo}.index;
                                CorrCoefstructShuffHPc{hpccnt, 1}.pfcindex(1,:) = pairdata(1,i).PFCidx(1,:);
                                CorrCoefstructShuffHPc{hpccnt, 1}.ca1index(1,:) = HPcCA1struct{mo}.index;
                            end
                        end
                    end
                    
                elseif pairdata(1,i).PFCidx(1,1) == 4%Nadal
                    for mo = 1:length(NdlCA1struct); %loop over all the HP CA1 cells for the shuffled struct
                        if isempty(find(str2num(sprintf('%d', NdlCA1struct{mo,1}.index)) == CA1ShortList));
                            if isempty(find(str2num(sprintf('%d', NdlCA1struct{mo,1}.index)) == CA1ShortListb)); %loop over CA1 cells that are no significantly correlated with the current PFC cell
                                mocnt = mocnt +1;
                                ndlcnt = ndlcnt +1;
                                for tru = 1:length(trajplot{1}); %for each trajectory
                                    clear r p minlen;
                                    minlen = min(length(trajplot{1}{tru}(:,5)), length(NdlCA1struct{mo}.trajdata{tru}(:,5)));
                                    [r, p] = corrcoef(trajplot{1}{tru}(1:minlen,5), NdlCA1struct{mo}.trajdata{tru}(1:minlen,5), 'rows', 'pairwise'); %don't use rows if either has a nan
                                    CorrCoefstructShuff{mocnt,1}.r(tru,1) = r(2,1); %store corr coefs
                                    CorrCoefstructShuff{mocnt,1}.p(tru,1) = p(2,1); %store corr coefs p vals
                                    CorrCoefstructShuffNdl{ndlcnt,1}.r(tru,1) = r(2,1); %store corr coefs
                                    CorrCoefstructShuffNdl{ndlcnt,1}.p(tru,1) = p(2,1); %store corr coefs p vals
                                end
                                CorrCoefstructShuff{mocnt, 1}.pfcindex(1,:) = pairdata(1,i).PFCidx(1,:);
                                CorrCoefstructShuff{mocnt, 1}.ca1index(1,:) = NdlCA1struct{mo}.index;
                                CorrCoefstructShuffNdl{ndlcnt, 1}.pfcindex(1,:) = pairdata(1,i).PFCidx(1,:);
                                CorrCoefstructShuffNdl{ndlcnt, 1}.ca1index(1,:) = NdlCA1struct{mo}.index;
                            end
                        end
                    end
                end
            end
            
            
            %-----------------------------------------------------------------------------------
            if plotPFCCA1maps ==1;
                for ru = 1:length(trajplot);
                    subplot(2,length(trajplot),ru); hold on;
                    for o = 1:length(trajplot{ru}); %plot each trajectory
                        plot(trajplot{ru}{o}(:,5),[clr{o} '.-'],'Linewidth',2);
                    end
                    if ru ==1; %add legend to pfc plot
                        legfix = legend('OL','IL','OR','IR');
                        legend(legfix, 'boxoff');
                        hText = findobj(legfix, 'type', 'text');
                        set(hText(4),'color',  [0 0 1]);
                        set(hText(3),'color',  [1 0 0]);
                        set(hText(2),'color',  [0 1 0]);
                        set(hText(1),'color',  [1 0 1]);
                        linesInPlot = findobj(legfix, 'type', 'line');
                        set(linesInPlot(1:8),'XData',[0 0]);
                    elseif ru-1 <= length(pairdata(1,i).CA1posidxs(:,1)); %blue the pos indxs
                        typclr = [0 .7 .93];
                        set(gca,'Color',typclr);
                        set(gcf, 'InvertHardCopy', 'off');
                    else %yellow the neg indxs
                        typclr = [1 1 .67];
                        set(gca,'Color',typclr);
                        set(gcf, 'InvertHardCopy', 'off');
                    end
                    set(gca,'xtick',[])
                    set(gca,'xticklabel',[])
                    set(gca,'ytick',[])
                    set(gca,'yticklabel',[])
                    %firing map
                    subplot(2, length(mapplot), ru+ length(trajplot)); hold on;
                    imagesc(mapplot{ru});
                    if ru ==1;
                        tit = {'PFC'; sprintf('%s %d %d %d', animals{pairdata(1,i).PFCidx(1,1)},pairdata(1,i).PFCidx(1,2), pairdata(1,i).PFCidx(1,3), pairdata(1,i).PFCidx(1,4))};  %label w index
                    else
                        tit = {'ca1'; sprintf('%s %d %d %d', animals{pairdata(1,i).CA1sigidxs(ru-1,1)},pairdata(1,i).CA1sigidxs(ru-1,2), pairdata(1,i).CA1sigidxs(ru-1,3), pairdata(1,i).CA1sigidxs(ru-1,4))}; %label w index
                    end
                    set(gca,'xtick',[])
                    set(gca,'xticklabel',[])
                    set(gca,'ytick',[])
                    set(gca,'yticklabel',[])
                    title(tit)
                end
                figfile = [figdir,savefilename,'/','SigCorrMaps',animals{pairdata(1,i).PFCidx(1,1)}, num2str(i)];
                if savefigs==1
                    print('-djpeg', figfile);
                end
                if ~cyclemaps == 0
                    keyboard;
                end
                close
            end
        end %if pfc cell didn't have peakthresh rate
    end %for each pfc cell
    %_____________________________________________________________________________________
    %gather means of trajs  and run stats
    
    if plotCorrCoef == 1;
        for yu = 1:length(CorrCoefstructPOS); %collect mean corrcoef for each pos pair into list
            posR(yu,1) = mean(CorrCoefstructPOS{yu,1}.r(:));
        end
        
        for qu = 1:length(CorrCoefstructNEG); %collect mean corrcoef for each neg pair into list
            negR(qu,1) = mean(CorrCoefstructNEG{qu,1}.r(:));
        end
        
        for wu = 1:length(CorrCoefstructShuff); %collect mean corrcoef for each pair into list
            allRShuff(wu,1) = mean(CorrCoefstructShuff{wu,1}.r(:));
        end
        
        
        %gather means of trajs for individual animals
        hpa = 0; hpb = 0; hpc = 0; ndl = 0;
        for fu = 1:length(CorrCoefstructPOS); %POS
            if CorrCoefstructPOS{fu,1}.pfcindex(1,1) == 1;
                hpa = hpa+1;
                CorrCoefstructPOSHPa{hpa,1} = CorrCoefstructPOS{fu,1};
                posRHPa(hpa, 1) = mean(CorrCoefstructPOS{fu,1}.r(:));
            elseif CorrCoefstructPOS{fu,1}.pfcindex(1,1) == 2;
                hpb = hpb+1;
                CorrCoefstructPOSHPb{hpb,1} = CorrCoefstructPOS{fu,1};
                posRHPb(hpb, 1) = mean(CorrCoefstructPOS{fu,1}.r(:));
            elseif CorrCoefstructPOS{fu,1}.pfcindex(1,1) == 3;
                hpc = hpc+1;
                CorrCoefstructPOSHPc{hpc,1} = CorrCoefstructPOS{fu,1};
                posRHPc(hpc, 1) = mean(CorrCoefstructPOS{fu,1}.r(:));
            elseif CorrCoefstructPOS{fu,1}.pfcindex(1,1) == 4;
                ndl = ndl+1;
                CorrCoefstructPOSNdl{ndl,1} = CorrCoefstructPOS{fu,1};
                posRNdl(ndl, 1) = mean(CorrCoefstructPOS{fu,1}.r(:));
            end
        end
        
        hpa = 0; hpb = 0; hpc = 0; ndl = 0; fu = 0;
        for fu = 1:length(CorrCoefstructNEG); %NEG
            if CorrCoefstructNEG{fu,1}.pfcindex(1,1) == 1;
                hpa = hpa+1;
                CorrCoefstructNEGHPa{hpa,1} = CorrCoefstructNEG{fu,1};
                negRHPa(hpa, 1) = mean(CorrCoefstructNEG{fu,1}.r(:));
            elseif CorrCoefstructNEG{fu,1}.pfcindex(1,1) == 2;
                hpb = hpb+1;
                CorrCoefstructNEGHPb{hpb,1} = CorrCoefstructNEG{fu,1};
                negRHPb(hpb, 1) = mean(CorrCoefstructNEG{fu,1}.r(:));
            elseif CorrCoefstructNEG{fu,1}.pfcindex(1,1) == 3;
                hpc = hpc+1;
                CorrCoefstructNEGHPc{hpc,1} = CorrCoefstructNEG{fu,1};
                negRHPc(hpc, 1) = mean(CorrCoefstructNEG{fu,1}.r(:));
            elseif CorrCoefstructNEG{fu,1}.pfcindex(1,1) == 4;
                ndl = ndl+1;
                CorrCoefstructNEGNdl{ndl,1} = CorrCoefstructNEG{fu,1};
                negRNdl(ndl, 1) = mean(CorrCoefstructNEG{fu,1}.r(:));
            end
        end
        
        wu = 0;
        for wu = 1:length(CorrCoefstructShuffHPa); %collect mean corrcoef for each pair into list
            allRShuffHPa(wu,1) = mean(CorrCoefstructShuffHPa{wu,1}.r(:));
        end
        
        wu = 0;
        for wu = 1:length(CorrCoefstructShuffHPb); %collect mean corrcoef for each pair into list
            allRShuffHPb(wu,1) = mean(CorrCoefstructShuffHPb{wu,1}.r(:));
        end
        
        wu = 0;
        for wu = 1:length(CorrCoefstructShuffHPc); %collect mean corrcoef for each pair into list
            allRShuffHPc(wu,1) = mean(CorrCoefstructShuffHPc{wu,1}.r(:));
        end
        
        wu = 0;
        for wu = 1:length(CorrCoefstructShuffNdl); %collect mean corrcoef for each pair into list
            allRShuffNdl(wu,1) = mean(CorrCoefstructShuffNdl{wu,1}.r(:));
        end
        
        posR = posR(~isnan(posR));
        negR = negR(~isnan(negR));
        
        posRHPa = posRHPa(~isnan(posRHPa));
        negRHPa = negRHPa(~isnan(negRHPa));
        
        posRHPb = posRHPb(~isnan(posRHPb));
        negRHPb = negRHPb(~isnan(negRHPb));
        
        posRHPc = posRHPc(~isnan(posRHPc));
        negRHPc = negRHPc(~isnan(negRHPc));
        
        posRNdl = posRNdl(~isnan(posRNdl));
        negRNdl = negRNdl(~isnan(negRNdl));
        
        allRShuff = allRShuff(~isnan(allRShuff));
        allRShuffHPa = allRShuffHPa(~isnan(allRShuffHPa));
        allRShuffHPb = allRShuffHPb(~isnan(allRShuffHPb));
        allRShuffHPc = allRShuffHPc(~isnan(allRShuffHPc));
        allRShuffNdl = allRShuffNdl(~isnan(allRShuffNdl));
        
        
        [hpos ppos] = kstest2(posR,allRShuff);
        [hneg pneg] = kstest2(negR,allRShuff);
        rspos = ranksum(posR,allRShuff);
        rsneg = ranksum(negR,allRShuff);
        
        
        [hposHPa pposHPa] = kstest2(posRHPa,allRShuffHPa);
        [hnegHPa pnegHPa] = kstest2(negRHPa,allRShuffHPa);
        rsposHPa = ranksum(posRHPa,allRShuffHPa);
        rsnegHPa = ranksum(negRHPa,allRShuffHPa);
        
        [hposHPb pposHPb] = kstest2(posRHPb,allRShuffHPb);
        [hnegHPb pnegHPb] = kstest2(negRHPb,allRShuffHPb);
        rsposHPb = ranksum(posRHPb,allRShuffHPb);
        rsnegHPb = ranksum(negRHPb,allRShuffHPb);
        
        
        [hposHPc pposHPc] = kstest2(posRHPc,allRShuffHPc);
        [hnegHPc pnegHPc] = kstest2(negRHPc,allRShuffHPc);
        rsposHPc = ranksum(posRHPc,allRShuffHPc);
        rsnegHPc = ranksum(negRHPc,allRShuffHPc);
        
        
        [hposNdl pposNdl] = kstest2(posRNdl,allRShuffNdl);
        [hnegNdl pnegNdl] = kstest2(negRNdl,allRShuffNdl);
        rsposNdl = ranksum(posRNdl,allRShuffNdl);
        rsnegNdl = ranksum(negRNdl,allRShuffNdl);
        
        
        %______________________________________________________________________________________
        %plot stats
        %             figure(1)
        %             hold on
        %             bar(1, mean(posR), 'FaceColor', [0 .7 .93], 'EdgeColor', 'none')
        %             bar(2, mean(allRShuff),'k', 'EdgeColor', 'none')
        %             bar(3, mean(negR),'FaceColor', [1 1 .67], 'EdgeColor', 'none')
        %             errorbar2([1 2 3], [mean(posR) mean(allRShuff) mean(negR) ],  [std(posR) std(allRShuff) stderr(negR)] , 0.3, 'k')
        %             %         xlim([0.3 2.7])
        %             set(gca, 'fontsize', 24)
        %             set(gca, 'xtick', [1 2 3], 'xticklabel', {'POS', 'AllShuff', 'NEG'})
        %             ylabel('CorrCoef')
        %             title({['Rip Corr PFC-CA1 All Animals']; sprintf('ranksum posP = %d  negP = %d',rspos, rsneg); sprintf('KStest posP = %d negP = %d ',ppos,pneg)})
        
        %each animal
        subplot(2,2,1);
        hold on
        bar(1, mean(posRHPa), 'FaceColor', [0 .7 .93], 'EdgeColor', 'none')
        bar(2, mean(allRShuffHPa),'k', 'EdgeColor', 'none')
        bar(3, mean(negRHPa),'FaceColor', [1 1 .67], 'EdgeColor', 'none')
        errorbar2([1 2 3], [mean(posRHPa) mean(allRShuffHPa) mean(negRHPa) ],  [std(posRHPa) std(allRShuffHPa) stderr(negRHPa)] , 0.3, 'k')
        %         xlim([0.3 2.7])
        set(gca, 'fontsize', 24)
        set(gca, 'xtick', [1 2 3], 'xticklabel', {'POS', 'AllShuff', 'NEG'})
        ylabel('CorrCoef')
        title({['Rip Corr PFC-CA1 HPA']; sprintf('ranksum posP = %d  negP = %d',rsposHPa, rsnegHPa); sprintf('KStest posP = %d negP = %d ',pposHPa,pnegHPa)})
        
        subplot(2,2,2);
        hold on
        bar(1, mean(posRHPb), 'FaceColor', [0 .7 .93], 'EdgeColor', 'none')
        bar(2, mean(allRShuffHPb),'k', 'EdgeColor', 'none')
        bar(3, mean(negRHPb),'FaceColor', [1 1 .67], 'EdgeColor', 'none')
        errorbar2([1 2 3], [mean(posRHPb) mean(allRShuffHPb) mean(negRHPb) ],  [std(posRHPb) std(allRShuffHPb) stderr(negRHPb)] , 0.3, 'k')
        %         xlim([0.3 2.7])
        set(gca, 'fontsize', 24)
        set(gca, 'xtick', [1 2 3], 'xticklabel', {'POS', 'AllShuff', 'NEG'})
        ylabel('CorrCoef')
        title({['Rip Corr PFC-CA1 HPb']; sprintf('ranksum posP = %d  negP = %d',rsposHPb, rsnegHPb); sprintf('KStest posP = %d negP = %d ',pposHPb,pnegHPb)})
        
        subplot(2,2,3);
        hold on
        bar(1, mean(posRHPc), 'FaceColor', [0 .7 .93], 'EdgeColor', 'none')
        bar(2, mean(allRShuffHPc),'k', 'EdgeColor', 'none')
        bar(3, mean(negRHPc),'FaceColor', [1 1 .67], 'EdgeColor', 'none')
        errorbar2([1 2 3], [mean(posRHPc) mean(allRShuffHPc) mean(negRHPc) ],  [std(posRHPc) std(allRShuffHPc) stderr(negRHPc)] , 0.3, 'k')
        %         xlim([0.3 2.7])
        set(gca, 'fontsize', 24)
        set(gca, 'xtick', [1 2 3], 'xticklabel', {'POS', 'AllShuff', 'NEG'})
        ylabel('CorrCoef')
        title({['Rip Corr PFC-CA1 HPc']; sprintf('ranksum posP = %d  negP = %d',rsposHPc, rsnegHPc); sprintf('KStest posP = %d negP = %d ',pposHPc,pnegHPc)})
        
        subplot(2,2,4);
        hold on
        bar(1, mean(posRNdl), 'FaceColor', [0 .7 .93], 'EdgeColor', 'none')
        bar(2, mean(allRShuffNdl),'k', 'EdgeColor', 'none')
        bar(3, mean(negRNdl),'FaceColor', [1 1 .67], 'EdgeColor', 'none')
        errorbar2([1 2 3], [mean(posRNdl) mean(allRShuffNdl) mean(negRNdl) ],  [std(posRNdl) std(allRShuffNdl) stderr(negRNdl)] , 0.3, 'k')
        %         xlim([0.3 2.7])
        set(gca, 'fontsize', 24)
        set(gca, 'xtick', [1 2 3], 'xticklabel', {'POS', 'AllShuff', 'NEG'})
        ylabel('CorrCoef')
        title({['Rip Corr PFC-CA1 Ndl']; sprintf('ranksum posP = %d  negP = %d',rsposNdl, rsnegNdl); sprintf('KStest posP = %d negP = %d ',pposNdl,pnegNdl)})
        

        figfile = [figdir,savefilename,'/','SigCorrCoefsALLshuff', num2str(ee)];
        if savefigs==1
            print('-djpeg', figfile);
        end
        if ~cyclemaps == 0
            keyboard;
        end
        close
    end
    
end %run


