

%update: rewriting the script because it's horrible.. look at DFS_DR_spatialvsripcorr.m for the newest version

%05/06/14... I didn't want to rewrite a lot of the code so this is a bit of a hack.  first, run the code with UncorrHACK == 1... this will generate and save the correct uncorrelated pair structure. then, run as UncorrHACK == 0, which will run 
%the correct positive and negative, and then load in the correct uncorrelated pairs.. ask demetris if you have any questions... i apologize in advance.



warning('off','all');
clear; close all;

runscript = 0;
savedata = runscript; % save data option - only works if runscript is also on
savefigs=1;
% plotstuff = 1; %plot anything
% plotEachcell = 1;
UncorrHACK = 0; %0 for real pos and neg... 1 for uncorrelated pairs.. see note at top
saveuncorr = UncorrHACK; % this is horrible coding. i'll pay for this later
cyclemaps =1;
peakthresh = 3;
fonttype = 'Arial';
titlesize = 16;
axissize = 16;
trajline = 4; %traj line width
% runnadal =0;
runsampling = 0; %run the sampling of 'shuffled' data to derive p values for pos/neg distributions..not using this anymore
runPFCCA1maps =1;
plotCorrCoef =1;
plotPFCCA1Trajs = 0;
plotPFCCA1maps =0;
plot2dPFC = 0;
mod1unmod0 = 2; %2 to skip saving the mod and unmod R structures.
% savedir = '/mnt/data25/sjadhav/HPExpt/ProcessedDataDR/';
savedir = '/mnt/data25/sjadhav/HPExpt/ProcessedDataDR/';
savefilename = 'RipCorrPairs_AllAn_041114'; %has to match saved data name. specify GLM or Corr pairs and filters used: velocity filter <=.....nrip >= (#tetrodes ripples detected); std > of ripple detection power
savefigfilename = 'RipCorrPairs_AllAn_050614'; %fig name
savefile = [savedir savefilename]; %area = 'PFC'; %clrunmod = 'r'; clrmod = 'b'; % PFC
figdir = '/mnt/data25/sjadhav/HPExpt/Figures_DR/';
loadrippos = '/mnt/data25/sjadhav/HPExpt/ProcessedDataDR/AllAn_PFCCA1_ripplepos_DR_vel5tet1'; %load the rip positions from DFSsj_HPexpt_getripalignspikingGRAllPosition_DR.m
% loadpairindices = '/home/droumis/MATLAB/allPFCCA1sigidxs_withcorr'; %load the correct indices!! specificy GLM or Corr Pairs in savefilename above
% loadpairindices = '/mnt/data25/sjadhav/HPExpt/HP_ProcessedData/allPFCCA1sigidxs'; %load the correct indices!! specificy GLM or Corr Pairs in savefilename above
% loadpairindices = '/mnt/data25/sjadhav/HPExpt/HP_ProcessedData/allPFCCA1sigidxs_mod'; %load the correct indices!! specificy GLM or Corr Pairs in savefilename above.. SANITY CHECK
loadpairindices = '/mnt/data25/sjadhav/HPExpt/HP_ProcessedData/April19spatindices'; %load the correct indices!! specificy GLM or Corr Pairs in savefilename above.. 

% crossvalindices = '/home/droumis/MATLAB/crossValidated'; %load the cross val indices..
% combineHPNdl = 1;
% Veqn = '>3';
% minV=str2num(Veqn(end));
% mintime = 3;
% traj = [1:4] ;
mkdir([figdir savefigfilename]);

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
    dayfilter = ''; %leave blank to take all days from HP animals and Ndl
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
    timefilter = { {'DFTFsj_getlinstate', '(($state ~= -1) & (abs($linearvel) >= 3))', 6}, {'DFTFsj_getriptimes','($nripples == 0)','tetfilter',riptetfilter,'minthresh',3} }; %DR added velocity filter.. trying to get ride of v high prococc in data..update, i dont think this is necessary any more bc
    ... Im using the linfield and mapfield structures that should have been generated using a speed filter
        
% Iterator
% --------
iterator = 'singlecellanal';

% Filter creation
% ----------------
%     spatf = createfilter('animal',animals,'epochs',runepochfilter,'cells',placecellfilter,'excludetime', timefilter,'iterator',iterator);
spatf = createfilter('animal',animals,'days',dayfilter,'epochs',runepochfilter,'cells',placecellfilter,'excludetime', timefilter,'iterator',iterator);

%do i need this? DR commented out 2/26/14... GIdeon doesn't use this either..
%     spatf = testexcludetimes(spatf, mintime); %removes epochs from analysis if all epoch excluded by excludetimes, mintime = 30

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
%--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

%CREATE CA1 STRUCTS FOR EACH AN. ONE THAT HAS ALL THE DATA AND ONE THAT IS A SHORT LIST (AN, NO EP) OF INDICES..
%     These will be used for finding data for all eps of each pair, and seperately excluding cells from the shuffle if they are sig noise corr with curr pfc cell

ca1cnt = 0; ca1ndl = 0; ca1hpa =0; ca1hpb = 0; ca1hpc = 0;
%     for ani = 1:3; % loop thru HP anims
%         for ai = 1:length(flds(1,ani).output{1}); %loop over cells in anim
%             if strcmp(flds(1,ani).output{1}(1,ai).area,'iCA1') || strcmp(flds(1,ani).output{1}(1,ai).area,'CA1'); %if ca1
%                 if length(flds(1,ani).output{1}(1,ai).trajdata) == 4; %if it has all 4 trajs
%                     ca1cnt = ca1cnt +1;
%                     HPCA1struct{ca1cnt,1}= flds(1,ani).output{1}(1,ai); %store ca1 data
%                     HPCA1structShortList(ca1cnt, 1) = str2num(sprintf('%d',[flds(1,ani).output{1}(1,ai).index([1 3 4])])); %don't use this struct/ doesnt have an prefix
%                 end
%             end
%         end
%     end
%     HPCA1structShortListCombEp = unique(HPCA1structShortList(:),'stable'); %get list of cells across epochs

%create ca1 struct for nadal
ani = 4;
for ai = 1:length(flds(1,ani).output{1}); %loop over cells in anim
    if strcmp(flds(1,ani).output{1}(1,ai).area,'iCA1') || strcmp(flds(1,ani).output{1}(1,ai).area,'CA1');
        if length(flds(1,ani).output{1}(1,ai).trajdata) == 4; %if it has all 4 trajs
            ca1ndl = ca1ndl +1;
            NdlCA1struct{ca1ndl,1}= flds(1,ani).output{1}(1,ai); %store ca1 data
            NdlCA1structShortList(ca1ndl, 1) = str2num(sprintf('%d',[4 flds(1,ani).output{1}(1,ai).index([1 3 4])])); %[animal day tet cell]
        end
    end
end

NdlCA1structShortListCombEp = unique(NdlCA1structShortList(:),'stable'); %get list of cells across epochs

%create ca1 struct for HPa
ani = 1; %
% HPaCA1structShortListDay = cell(1,8);
for ai = 1:length(flds(1,ani).output{1}); %loop over cells in anim
    if strcmp(flds(1,ani).output{1}(1,ai).area,'iCA1') || strcmp(flds(1,ani).output{1}(1,ai).area,'CA1');
        if length(flds(1,ani).output{1}(1,ai).trajdata) == 4; %if it has all 4 trajs
            ca1hpa = ca1hpa +1;
            HPaCA1struct{ca1hpa,1}= flds(1,ani).output{1}(1,ai); %store ca1 data
            HPaCA1structShortList(ca1hpa, 1) = str2num(sprintf('%d',[1 flds(1,ani).output{1}(1,ai).index([1 3 4])]));
            %             HPaCA1structShortListDay{flds(1,ani).output{1}(1,ai).index(1)} = [HPaCA1structShortListDay{flds(1,ani).output{1}(1,ai).index(1)}; str2num(sprintf('%d',[1 flds(1,ani).output{1}(1,ai).index([1 3 4])]))];
        end
    end
end

HPaCA1structShortListCombEp = unique(HPaCA1structShortList(:),'stable'); %get list of cells across epochs

%create ca1 struct for HPb
ani = 2;
for ai = 1:length(flds(1,ani).output{1}); %loop over cells in anim
    if strcmp(flds(1,ani).output{1}(1,ai).area,'iCA1') || strcmp(flds(1,ani).output{1}(1,ai).area,'CA1');
        if length(flds(1,ani).output{1}(1,ai).trajdata) == 4; %if it has all 4 trajs
            ca1hpb = ca1hpb +1;
            HPbCA1struct{ca1hpb,1}= flds(1,ani).output{1}(1,ai); %store ca1 data
            HPbCA1structShortList(ca1hpb, 1) = str2num(sprintf('%d',[2 flds(1,ani).output{1}(1,ai).index([1 3 4])]));
        end
    end
end

HPbCA1structShortListCombEp = unique(HPbCA1structShortList(:),'stable'); %get list of cells across epochs

%create ca1 struct for HPc
ani = 3;
for ai = 1:length(flds(1,ani).output{1}); %loop over cells in anim
    if strcmp(flds(1,ani).output{1}(1,ai).area,'iCA1') || strcmp(flds(1,ani).output{1}(1,ai).area,'CA1');
        if length(flds(1,ani).output{1}(1,ai).trajdata) == 4; %if it has all 4 trajs
            ca1hpc = ca1hpc +1;
            HPcCA1struct{ca1hpc,1}= flds(1,ani).output{1}(1,ai); %store ca1 data
            HPcCA1structShortList(ca1hpc, 1) = str2num(sprintf('%d',[3 flds(1,ani).output{1}(1,ai).index([1 3 4])]));
        end
    end
end

HPcCA1structShortListCombEp = unique(HPcCA1structShortList(:),'stable'); %get list of cells across epochs
%--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

if savedata == 1
    clear saveuncorr UncorrHACK runsampling mod1unmod0 crossvalindices savefigfilename runscript  savedata plotstuff cyclemaps savefigs plottrajs runnadal plotEachcell plotPFCCA1maps plotPFCCA1Trajs plotCorrCoef ripposfile loadrippos figdir loadpairindices savefilename savedir runPFCCA1maps peakthresh  fonttype titlesize axissize trajline plot2dPFC
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


% --------------------------------------------------------------------

if savefigs == 1;
    mkdir(figdir,savefilename)
end

%prepping data structs
%________________________________________________________________________
if runPFCCA1maps ==1;
    pfclistidx = []; ju = 0; fldsLIST = []; pu =0;  cnt = 0; trajdataall =[]; CorrCoefstruct = []; clear fldsdata.output fldsLIST fldsY fldsZ fldsShortList mfldsdata.output mfldsLIST mfldsY mfldsZ mfldsShortList  allnegPFCCA1listcheck allposPFCCA1listcheck
    allposPFCCA1listcheck = []; allnegPFCCA1listcheck = []; pairdata = {}; pairdataall = {};
    clr = {'b','r','g','m','c','y','k','r'};
    load(loadpairindices); %combined epochs
       
    %     load(crossvalindices);
    
    %     allPFCCA1sigidxs = allPFCCA1sigidxs(crossValidated); %only use cross validated sets
    
    figdir = '/mnt/data25/sjadhav/HPExpt/Figures_DR/';
    countnow = 0;
    glmpaircnt = 0;
    corrsetcnt = 0;
    
                
%%% update April 14... GR provided new indices in another different format with each cell being a ca1pfc pair... need to change back to work with code, again
% __________________________________________________________________________________________________    
%convert the gr/sj struct to the older sj format....
%first collect all data in other format.. regardless of significance
%get number and indices of unique pfc cells
for i = 1:length(corrindsForSpatial);
    listedPFCfrompairs(i,:) = corrindsForSpatial(i).index([1 2 5 6]); %collect all pfc indices
    ShortlistedPFCfrompairs(i) = str2num(sprintf('%d',listedPFCfrompairs(i,:))); %shortlist indices for matching and unique
end
uniquePFCcells = unique(ShortlistedPFCfrompairs,'stable'); %list unique indices for pairdata # 
for i = 1:length(uniquePFCcells);
    pfcmatch = find(uniquePFCcells(i) == ShortlistedPFCfrompairs,1,'first'); %get corresponding unique index in shortlist
    pairdataall(i).PFCidx = listedPFCfrompairs(pfcmatch,:); %use the unconcatenated indices for pairdata
    pairdataall(i).CA1idx = [];
    pairdataall(i).pvals = [];
    pairdataall(i).betas = [];
end

%match pairs to pfc cell index in pairdata struct and inser ca1 indices, p vals, r vals
for i = 1:length(corrindsForSpatial);
    currpfcindpair = str2num(sprintf('%d',corrindsForSpatial(i).index([1 2 5 6])));
    currpfcmatch = find(currpfcindpair == uniquePFCcells,1,'first'); %find the corresponding pfc cell struct to put this data in
    pairdataall(currpfcmatch).CA1idx = [pairdataall(currpfcmatch).CA1idx; corrindsForSpatial(i).index([1 2 3 4])];
    pairdataall(currpfcmatch).pvals = [pairdataall(currpfcmatch).pvals; corrindsForSpatial(i).allp_epcomb];
    pairdataall(currpfcmatch).betas = [pairdataall(currpfcmatch).betas; corrindsForSpatial(i).allr_epcomb];
end %      

%Now weed out the significant p val pairs
allsigpaircount = 0;
for i = 1:length(pairdataall);
    if ~isempty(pairdataall(i).pvals);
        corrsetcnt = corrsetcnt +1;
        
        if UncorrHACK == 1;
                        %DR THIS IS IMPORTANT. writing this hack to spit out uncorrelated cells as the positive group, so that I don't have to rewrite a lot of the code. see note at top.
            pairdata(corrsetcnt).CA1sigidxs = pairdataall(i).CA1idx(pairdataall(i).pvals>0.05,:);
            pairdata(corrsetcnt).sigbetas = pairdataall(i).betas(pairdataall(i).pvals>0.05,1);
            allsigpaircount = allsigpaircount + length(pairdata(corrsetcnt).sigbetas);
            
            pairdata(corrsetcnt).CA1posidxs = pairdata(corrsetcnt).CA1sigidxs;
            pairdata(corrsetcnt).CA1negidxs = pairdata(corrsetcnt).CA1sigidxs;
            pairdata(corrsetcnt).sigbetaspos = pairdata(corrsetcnt).sigbetas;
            pairdata(corrsetcnt).sigbetasneg = pairdata(corrsetcnt).sigbetas;
            
            pairdata(corrsetcnt).PFCidx = pairdataall(i).PFCidx;
            pairdata(corrsetcnt).CA1idx = pairdataall(i).CA1idx;
            pairdata(corrsetcnt).pvals = pairdataall(i).pvals;
            pairdata(corrsetcnt).betas = pairdataall(i).betas;
            
        else
            
        %USE these to get actual negative and positive... 
            pairdata(corrsetcnt).CA1sigidxs = pairdataall(i).CA1idx(pairdataall(i).pvals<0.05,:);
            pairdata(corrsetcnt).sigbetas = pairdataall(i).betas(pairdataall(i).pvals<0.05,1);
            allsigpaircount = allsigpaircount + length(pairdata(corrsetcnt).sigbetas); %84 sig pairs in all
            
            pairdata(corrsetcnt).CA1posidxs = pairdata(corrsetcnt).CA1sigidxs(pairdata(corrsetcnt).sigbetas>0,:);
            pairdata(corrsetcnt).CA1negidxs = pairdata(corrsetcnt).CA1sigidxs(pairdata(corrsetcnt).sigbetas<0,:);
            pairdata(corrsetcnt).sigbetaspos = pairdata(corrsetcnt).sigbetas(pairdata(corrsetcnt).sigbetas>0,:);
            pairdata(corrsetcnt).sigbetasneg = pairdata(corrsetcnt).sigbetas(pairdata(corrsetcnt).sigbetas<0,:);
            
            pairdata(corrsetcnt).PFCidx = pairdataall(i).PFCidx;
            pairdata(corrsetcnt).CA1idx = pairdataall(i).CA1idx;
            pairdata(corrsetcnt).pvals = pairdataall(i).pvals;
            pairdata(corrsetcnt).betas = pairdataall(i).betas;

        end

    else
        'pfc cells without any ca1 cells from indices'
         [i pairdataall(i).PFCidx] %this shouldn't spit out any pfc cells now that we're only looking at corr pairs
    end
end
% __________________________________________________________________________________________________    

    %this is using the old index format
% __________________________________________________________________________________________________    
%     %convert the gr struct to sj format
%     for i = 1:length(allPFCCA1sigidxs)
%         countnow = countnow + 1;
%         
%         if ~isempty(allPFCCA1sigidxs(i).pvals);
%             glmpaircnt = glmpaircnt +1;
%             %
%             %GLM pairs_____________________________________
%             %         pairdata(glmpaircnt).CA1sigidxs = allPFCCA1sigidxs(i).CA1idx(allPFCCA1sigidxs(i).pvals<0.05,:);
%             %         pairdata(glmpaircnt).sigbetas = allPFCCA1sigidxs(i).betas(allPFCCA1sigidxs(i).pvals<0.05,1);
%             %
%             %         pairdata(glmpaircnt).CA1posidxs = pairdata(glmpaircnt).CA1sigidxs(pairdata(glmpaircnt).sigbetas>0,:);
%             %         pairdata(glmpaircnt).CA1negidxs = pairdata(glmpaircnt).CA1sigidxs(pairdata(glmpaircnt).sigbetas<0,:);
%             %         pairdata(glmpaircnt).sigbetaspos = pairdata(glmpaircnt).sigbetas(pairdata(glmpaircnt).sigbetas>0,:);
%             %         pairdata(glmpaircnt).sigbetasneg = pairdata(glmpaircnt).sigbetas(pairdata(glmpaircnt).sigbetas<0,:);
%             %
%             %         pairdata(glmpaircnt).PFCidx = allPFCCA1sigidxs(i).PFCidx;
%             %         pairdata(glmpaircnt).CA1idx = allPFCCA1sigidxs(i).CA1idx;
%             %         pairdata(glmpaircnt).pvals = allPFCCA1sigidxs(i).pvals;
%             %         pairdata(glmpaircnt).betas = allPFCCA1sigidxs(i).betas;
%             
%             %         %Corr pairs_____________________________________
%             
%             %to use all pairs,  not just sig_____
%             %             pairdata(glmpaircnt).CA1sigidxs = allPFCCA1sigidxs(i).CA1idx(:,:);
%             %             pairdata(glmpaircnt).sigbetas = allPFCCA1sigidxs(i).corrrvals(:,1);
%             
%             %sig pairs----
%             pairdata(glmpaircnt).CA1sigidxs = allPFCCA1sigidxs(i).CA1idx(allPFCCA1sigidxs(i).corrpvals<0.05,:);
%             pairdata(glmpaircnt).sigbetas = allPFCCA1sigidxs(i).corrrvals(allPFCCA1sigidxs(i).corrpvals<0.05,1);
%             
%             pairdata(glmpaircnt).CA1posidxs = pairdata(glmpaircnt).CA1sigidxs(pairdata(glmpaircnt).sigbetas>0,:);
%             pairdata(glmpaircnt).CA1negidxs = pairdata(glmpaircnt).CA1sigidxs(pairdata(glmpaircnt).sigbetas<0,:);
%             pairdata(glmpaircnt).sigbetaspos = pairdata(glmpaircnt).sigbetas(pairdata(glmpaircnt).sigbetas>0,:);
%             pairdata(glmpaircnt).sigbetasneg = pairdata(glmpaircnt).sigbetas(pairdata(glmpaircnt).sigbetas<0,:);
%             
%             pairdata(glmpaircnt).PFCidx = allPFCCA1sigidxs(i).PFCidx;
%             pairdata(glmpaircnt).CA1idx = allPFCCA1sigidxs(i).CA1idx;
%             pairdata(glmpaircnt).pvals = allPFCCA1sigidxs(i).corrpvals;
%             pairdata(glmpaircnt).betas = allPFCCA1sigidxs(i).corrrvals;
%             %
% 
%         else
%             [i allPFCCA1sigidxs(i).PFCidx]
%         end
%     end
% __________________________________________________________________________________________________    
    
    %     pairdata = allPFCCA1sigidxs;
    for yu = 1:length(pairdata);
        PFCShortList(yu,1) = str2num(sprintf('%d', pairdata(1,yu).PFCidx(:)));
        %         PFClistindex(yu,:) = pairdata(1,yu).PFCidx(:);
        %         PFCY = sprintf('%d', PFClistindex(yu,:)); %make string out of all columns AN DAY TET CELL
        %         PFCZ = str2num(PFCY); %convert string to num
        %         PFCShortList(yu,1) =  PFCZ; %save for later
    end
    
    %transform the flds traj data into a workable and searchable format
    cnt = 0;
    for anims = 1:length(animals);
        for bu = 1:length(flds(anims).output{1});
            cnt = cnt+1;
            fldsdata.output{cnt} = flds(anims).output{1}(1,bu); %collect all data across animals
            fldsShortList(cnt,1) = str2num(sprintf('%d', [anims, fldsdata.output{cnt}.index([1 3 4])]));
            
            %             fldsLIST(cnt,:) = [anims , fldsdata.output{cnt}.index([1 3 4])]; %collect indices W/O EP and add animal num in first column to match pairdata format
            %             fldsY = sprintf('%d', fldsLIST(cnt,:)); %make string out of all columns an day ep tet cell
            %             fldsZ = str2num(fldsY); %convert string to num
            %             fldsShortList(cnt,1) =  fldsZ; %save for later
        end
    end
    
    %transform the flds map data into a workable and searchable format
    cnt = 0;
    for anims = 1:length(animals);
        for bu = 1:length(flds(anims).output{1});
            cnt = cnt+1;
            mfldsdata.output{cnt} = flds(anims).output{1}(1,bu); %collect all data across animals
            mfldsShortList(cnt,1) = str2num(sprintf('%d',  [anims, mfldsdata.output{cnt}.index([1 3 4 ])]));
            
            %             mfldsLIST(cnt,:) = [anims, mfldsdata.output{cnt}.index([1 3 4 ])]; %collect indices W/O EP and add animal num in first column to match pairdata format
            %             mfldsY = sprintf('%d', mfldsLIST(cnt,:)); %make string out of all columns an day ep tet cell
            %             mfldsZ = str2num(mfldsY); %convert string to num
            %             mfldsShortList(cnt,1) =  mfldsZ; %save for later
        end
    end
    
    %FIND AND COLLECT TRAJ AND MAP DATA FOR EACH PFC CA1 PAIR IN CURRENT SET
    %--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    %____________________________________________________________________________________________________________________
    %loop through the PFCShortList. for each i pfc, find PFCShortList(i) in fldsShortList then use that index with fldsdata.output to get the
    %trajectory and map data from fldsdata.output and the map from pfmdata.output..Then for the length of the CA1sigind of the pairdata(1,i),
    %create short list of CA1 indices, then loop through each row to find match in fldsShortList then use that index and store the traj from
    %fldsdata.output and the map from pfmdata.output
    
    CorrCoefcombPOSCOUNT =0; CorrCoefcombPOSNEGCOUNT = 0; CorrCoefcombNEGCOUNT = 0;
    paircounter = 0; pospaircounter = 0; negpaircounter = 0; mocnt = 0; hpacnt = 0; hpbcnt = 0; hpccnt = 0; ndlcnt = 0; epochcounter = 0; AllCorrRip_Sp = []; pfcskipped3Hz =0; ca1skipped3Hz =0; ca1posgood = 0; ca1neggood = 0; pairsskippedPFC3Hz = 0;
    
    for i = 1:length(PFCShortList);
        poscelsinset = 0; negcelsinset = 0; CA1ShortList = []; CA1ShortListb =[];
        clear tmptrajplotcombtraj hpaca1shuffEP temptrajplot pfcmatchtraj pfcmatchmap trajplot mapplot CA1X  CA1Y CA1Z...
            CA1matchtraj CA1matchmap CA1bX CA1bY CA1bZ  CA1matchtrajb CA1matchmapb pfcAnDaEp CA1POSAnDaEp CA1NEGAnDaEp...
            trajplotallEps hpaca1shuffEP hpbca1shuffEP hpcca1shuffEP ndlca1shuffEP trajplotEpAvg mapplotEpAvg
        pfcmatchtraj = find(PFCShortList(i) == fldsShortList(:), 4); %find corresponding indices from psfdata. finds all epochs
        pfcmatchmap = find(PFCShortList(i) == mfldsShortList(:), 4); %find corresponding index from pfmdata. finds all epochs
        
        %>peakthreshHz Check prep struct__right now im using all the epochs of the pfc cell if it surpasses peakthreshHz in any epoch.
        alltraj = [];
        for ou = 1:length(pfcmatchtraj(:,1));
            for tu = 1:length(fldsdata.output{pfcmatchtraj(ou, 1)}.trajdata);
                alltraj = [alltraj; fldsdata.output{pfcmatchtraj(ou,1)}.trajdata{tu}(:,5)];
            end
        end
        
        if max(alltraj) > peakthresh; %use cell if >3Hz peak threshold. if not, this whole pfc-ca1 set will be skipped.
            for ou = 1:length(pfcmatchtraj(:,1)); %now storing all the epoch matches (max 4 run epochs for some of gideon's)
                trajplot{ou,1} = fldsdata.output{pfcmatchtraj(ou,1)}.trajdata; %store pfc traj data
                mapplot{ou,1} = mfldsdata.output{pfcmatchmap(ou,1)}.mapdata.smoothedspikerate; %store pfc map data.
                %                 pfcAnDaEp(ou,1) =  str2num(sprintf('%d',[fldsLIST(pfcmatchtraj(ou,1),1) fldsdata.output{pfcmatchtraj(ou,1)}.index([1 2])])); % make a struct of An Day Eps for current pfc cell for later compring epochs of paired ca1 cells so that only matching epochs get corrcoefs..not doing this anymore
            end
            
            %loop over the list of CA1 cells for the current pfc cell and make a short list. then loop through the shortlist, finding traj, maps, and then store
            %or plot data for each. Seperately finds/stores pos noise corr then neg noise corr pairs
            
            paircnt = 1; %start at 1 bc pfc cell is already in first column. first ca1 cell below at index 2
            if ~isempty(pairdata(1,i).CA1posidxs);
                for posi = 1:length(pairdata(1,i).CA1posidxs(:,1));
                    CA1ShortList(posi,1) = str2num(sprintf('%d', pairdata(1,i).CA1posidxs(posi,:)));
                    %
                    %                     CA1X = pairdata(1,i).CA1posidxs(posi,:);
                    %                     CA1Y = sprintf('%d', CA1X);
                    %                     CA1Z = str2num(CA1Y);
                    %                     CA1ShortList(posi,1) = CA1Z;
                    CA1matchtraj = find(CA1ShortList(posi) == fldsShortList(:),4); %taking all epochs.. 4 bc gideon max runs ep =4
                    %                     CA1matchmap = find(CA1ShortList(posi) == mfldsShortList(:),4);
                    %find corresponding index frompfmdata... update.. not using map index anymore as it seems to be equal to the traj index
                    alltraj = [];
                    for ju = 1:length(CA1matchtraj);
                        for tu = 1:length(fldsdata.output{CA1matchtraj(ju,1)}.trajdata);
                            alltraj = [alltraj; fldsdata.output{CA1matchtraj(ju,1)}.trajdata{tu}(:,5)];
                        end
                    end
                    if max(alltraj) > peakthresh; %use cell if >Hz peak threshold
                        %check to make sure CA1matchtraj and CA1matchmap only contain the same epochs as the paired pfc cell
                        %version DR5 or greater.. I'm no longer limiting the epochs to those that match the pfc cell. I decided this after looking at rate maps
                        %over epochs, which seem consistent, and the fact that the shuffled versions didn't filter for epoch matches
                        CA1matchtrajEP = [];
                        for ru = 1:length(CA1matchtraj(:,1));
                            %                             CA1POSAnDaEp(ru,1) = str2num(sprintf('%d',[fldsLIST(CA1matchtraj(ru,1),1)  fldsdata.output{CA1matchtraj(ru,1)}.index([1 2])])); %make a struct of An Day Eps for current ca1 cell for compring epochs of paired ca1 cells so that only matching epochs get corrcoefs..update.. not doing this anymore
                            %                             if ~isempty(find(CA1POSAnDaEp(ru,1) == pfcAnDaEp(:))); %if there is a matching epoch between the ca1 cell and pfc cell
                            CA1matchtrajEP = [CA1matchtrajEP; CA1matchtraj(ru,1)]; %add the matching ep to the new match list for data gather. now adds all eps
                            %                             end
                        end
                        paircnt = paircnt +1;
                        poscelsinset = poscelsinset +1;
                        for ru = 1:length(CA1matchtrajEP(:,1));
                            trajplot{ru, paircnt} = fldsdata.output{CA1matchtrajEP(ru,1)}.trajdata; %store CA1 traj data next to pfc data
                            mapplot{ru,paircnt} = mfldsdata.output{CA1matchtrajEP(ru,1)}.mapdata.smoothedspikerate; %store CA1 map data next to pfc data..replaced CA1matchtrajmap with CA1matchtrajEP
                        end
                        ca1posgood =  ca1posgood + 1;
                    else
                        ca1skipped3Hz = ca1skipped3Hz+1;
                    end
                end %CA1 pos short list and store
            else
                CA1ShortList =[];
            end %if pos not empty
            
            %neg CA1 data
            if ~isempty(pairdata(1,i).CA1negidxs);
                for negi = 1:length(pairdata(1,i).CA1negidxs(:,1));
                    CA1ShortListb(negi,1) =  str2num(sprintf('%d', pairdata(1,i).CA1negidxs(negi,:)));
                    %                     CA1bX = pairdata(1,i).CA1negidxs(negi,:);
                    %                     CA1bY = sprintf('%d', CA1bX);
                    %                     CA1bZ = str2num(CA1bY);
                    %                     CA1ShortListb(negi,1) = CA1bZ;
                    CA1matchtrajb = find(CA1ShortListb(negi) == fldsShortList(:), 4); %find corresponding index from psfdata. 4 bc gideon max runs ep =4.. finds ALL epochs
                    %                     CA1matchmapb = find(CA1ShortListb(negi) == mfldsShortList(:), 4); %find corresponding index from pfmdata
                    alltraj = [];
                    for juu = 1:length(CA1matchtrajb);
                        for tu = 1:length(fldsdata.output{CA1matchtrajb(juu,1)}.trajdata);
                            alltraj = [alltraj; fldsdata.output{CA1matchtrajb(juu,1)}.trajdata{tu}(:,5)];
                        end
                    end
                    if max(alltraj) > peakthresh; %use cell if 3Hz peak threshold
                        CA1matchtrajEPb = [];
                        for ruu = 1:length(CA1matchtrajb(:,1));
                            %                             CA1NEGAnDaEp(ruu,1) = str2num(sprintf('%d',[fldsLIST(CA1matchtrajb(ruu,1),1)  fldsdata.output{CA1matchtrajb(ruu,1)}.index([1 2])])); %make a struct of An Day Eps for current ca1 cell for compring epochs of paired ca1 cells so that only matching epochs get corrcoefs. update.. not doing this anymore
                            %                             if ~isempty(find(CA1NEGAnDaEp(ruu,1) == pfcAnDaEp(:))); %if there is a matching epoch between the ca1 cell and pfc cell
                            CA1matchtrajEPb = [CA1matchtrajEPb; CA1matchtrajb(ruu,1)]; %add the matching ep to the new match list for data gather
                            %                             end
                        end
                        paircnt = paircnt +1;
                        negcelsinset = negcelsinset +1;
                        for ruu = 1:length(CA1matchtrajEPb(:,1));
                            trajplot{ruu, paircnt} = fldsdata.output{CA1matchtrajEPb(ruu,1)}.trajdata; %store CA1 traj data next to pfc data
                            mapplot{ruu,paircnt} = mfldsdata.output{CA1matchtrajEPb(ruu,1)}.mapdata.smoothedspikerate; %store CA1 map data next to pfc data.. replaced map index with traj index
                        end
                        ca1neggood = ca1neggood + 1;
                    else
                        ca1skipped3Hz = ca1skipped3Hz +1;
                    end
                end %CA1 neg short list and store
            else
                CA1ShortListb =[];
            end %if pos not empty
            
            %store all traj data
            %             trajdataall{length(trajdataall)+1,1} = trajplot; %store each pfc set for corrcoef
            
            %nanmean the spatial firing across epochs__________________________________________________________________________________________________________________________________________________________________________________
            
            for trajcell = 1:length(trajplot(1,:)); %for each cell in traj plot
                celltraj1 = []; celltraj2 = []; celltraj3 = []; celltraj4 = []; mapcell =[];
                for epcell = 1:length(trajplot(~cellfun('isempty',trajplot(:,trajcell)),trajcell)); %takes all the epochs that exist for the  cell
                    %find the min length of the trajectories to be averaged
                    if isempty(celltraj1); %if it's the first epoch
                        %                         minlen1 = length(trajplot{epcell,trajcell}{1,1}(:,5));
                        %                         minlen2 = length(trajplot{epcell,trajcell}{1,2}(:,5));
                        %                         minlen3 = length(trajplot{epcell,trajcell}{1,3}(:,5));
                        %                         minlen4 = length(trajplot{epcell,trajcell}{1,4}(:,5));
                        %                         celltraj1 = [trajplot{epcell,trajcell}{1,1}(1:minlen1,5)]; % trajectory 1
                        %                         celltraj2 = [trajplot{epcell,trajcell}{1,2}(1:minlen2,5)]; % trajectory 2
                        %                         celltraj3 = [trajplot{epcell,trajcell}{1,3}(1:minlen3,5)]; %trajectory 3
                        %                         celltraj4 = [trajplot{epcell,trajcell}{1,4}(1:minlen4,5)]; % trajectory 4
                        %                         %Trajs
                        celltraj1 = [trajplot{epcell,trajcell}{1,1}(:,5)]; % trajectory 1
                        celltraj2 = [trajplot{epcell,trajcell}{1,2}(:,5)]; % trajectory 2
                        celltraj3 = [trajplot{epcell,trajcell}{1,3}(:,5)]; %trajectory 3
                        celltraj4 = [trajplot{epcell,trajcell}{1,4}(:,5)]; % trajectory 4
                        
                        %MAP DATA
                        mapminlenrow = length(mapplot{epcell,trajcell}(:,1));
                        mapminlencol = length(mapplot{epcell,trajcell}(1,:));
                        mapcell = mapplot{epcell,trajcell};
                        
                        
                    else%all the rest of the epochs
                        %MAPs.. since maps are 2Dimensional, I need to find the mean every epoch
                        mapminlenrow = min(length(mapcell(:,1)), length(mapplot{epcell,trajcell}(:,1))); %use the length from the shortest row
                        mapminlencol = min(length(mapcell(1,:)), length(mapplot{epcell,trajcell}(1,:))); %use the length from the shortest col
                        mapcell = (mapcell(1:mapminlenrow,1:mapminlencol) + mapplot{epcell,trajcell}(1:mapminlenrow,1:mapminlencol))./2; %take the elementwise mean of the same sized spatial matrices
                        
                        
                        %Trajs
                        minlen1 = min(length(celltraj1), length(trajplot{epcell,trajcell}{1,1}(:,5))); %use the length from the shortest trajectory across epochs for each trajectory
                        minlen2 = min(length(celltraj2), length(trajplot{epcell,trajcell}{1,2}(:,5)));
                        minlen3 = min(length(celltraj3), length(trajplot{epcell,trajcell}{1,3}(:,5)));
                        minlen4 = min(length(celltraj4), length(trajplot{epcell,trajcell}{1,4}(:,5)));
                        celltraj1 = [celltraj1(1:minlen1,:) trajplot{epcell,trajcell}{1,1}(1:minlen1,5)]; %horizontally concatenate trajectory 1
                        celltraj2 = [celltraj2(1:minlen2,:) trajplot{epcell,trajcell}{1,2}(1:minlen2,5)]; %horizontally concatenate trajectory 2
                        celltraj3 = [celltraj3(1:minlen3,:) trajplot{epcell,trajcell}{1,3}(1:minlen3,5)]; %horizontally concatenate trajectory 3
                        celltraj4 = [celltraj4(1:minlen4,:) trajplot{epcell,trajcell}{1,4}(1:minlen4,5)]; %horizontally concatenate trajectory 4
                        
                    end %first epoch or not
                end %all the epochs for each cell
                trajplotEpAvg{1, trajcell}{1,1}(:,5) = nanmean(celltraj1,2); %row wise (spatial bin) nanmean across epochs for each trajectory
                trajplotEpAvg{1, trajcell}{1,2}(:,5) = nanmean(celltraj2,2);
                trajplotEpAvg{1, trajcell}{1,3}(:,5) = nanmean(celltraj3,2);
                trajplotEpAvg{1, trajcell}{1,4}(:,5) = nanmean(celltraj4,2);
                
                mapplotEpAvg{1,trajcell} = mapcell; %collect the epoch averages maps for each cell in this set
                
            end %each cell
            
            trajplotallEps = trajplot; %saving old traj plot (not combined across epochs) for nostalgia
            
            %concatenating the trajectories as to get a single corrcoef later
            for eachcell = 1:length(trajplotEpAvg);
                combtraj = [];
                for eachtraj = 1:length(trajplotEpAvg{eachcell});
                    if eachtraj == 2 || eachtraj ==4; %inbound trajs
                        trajplotEpAvg{eachcell}{eachtraj} = flipud(trajplotEpAvg{eachcell}{eachtraj}); %flip the vector for the inbound trajectories so that the concatenated trajectories are continuous with how the animal ran.
                    end
                    combtraj = [combtraj; trajplotEpAvg{eachcell}{eachtraj}];
                end
                tmptrajplotcombtraj{eachcell}{1} = combtraj;
            end
            
            trajplot = tmptrajplotcombtraj; %redefining trajplot be averaged across epochs and concatenated trajectories.. update: now also with inbound trajs fliped for continuity with actual behavior
            
            %summing the pos and neg cells and doing  a corr coef for
            %pos~pfc neg~pfc and (pos-neg)~pfc
            %--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            normPOS=[]; normNEG = []; negcels = 0; poscels = 0; clear combnormPOS combnormNEG combnormPOSNEG combnormPOSarith combnormNEGarith combnormPOSNEGarith
            for cels = 2:length(trajplot);
                if cels-1 <= poscelsinset; %length(pairdata(1,i).CA1posidxs(:,1)); %pos cell
                    poscels = poscels +1;
                    %normalize by firing rate within cell
                    if poscels == 1; %if its the first cell.
                        normPOS(:,poscels) = trajplot{1,cels}{1}(:,5)./max(trajplot{1,cels}{1}(:,5));
                    else
                        minlen = min(length(normPOS(:,1)), length(trajplot{1,cels}{1}(:,5)));
                        normPOS = [normPOS(1:minlen,:) trajplot{1,cels}{1}(1:minlen,5)./max(trajplot{1,cels}{1}(1:minlen,5))];
                    end
                else cels-1 <= poscelsinset + negcelsinset;%neg cells
                    negcels = negcels +1;
                    if negcels ==1;
                        normNEG(:,negcels) = trajplot{1,cels}{1}(:,5)./max(trajplot{1,cels}{1}(:,5));
                    else
                        minlen = min(length(normNEG(:,1)), length(trajplot{1,cels}{1}(:,5)));
                        normNEG = [normNEG(1:minlen,:) trajplot{1,cels}{1}(1:minlen,5)./max(trajplot{1,cels}{1}(1:minlen,5))];
                    end
                    
                end
            end
            
            if poscelsinset >1; %length(pairdata(1,i).CA1posidxs(:,1)) > 1; %if there are more than 1 pos cells in this set
                CorrCoefcombPOSCOUNT = CorrCoefcombPOSCOUNT +1;
                combnormPOS = nansum(normPOS,2);
                minlen = min(length(trajplot{1,1}{1}(:,5)), length(combnormPOS));
                [r, p] = corrcoef(trajplot{1}{1}(1:minlen,5), combnormPOS(1:minlen,1), 'rows', 'pairwise'); %don't use rows if either has a nan
                CorrCoefcombPOS(CorrCoefcombPOSCOUNT,1) =  r(2,1);
            end
            
            if negcelsinset > 1; %length(pairdata(1,i).CA1negidxs(:,1)) >1; %if there are more than 1 neg cells in this set
                CorrCoefcombNEGCOUNT = CorrCoefcombNEGCOUNT +1;
                combnormNEG = nansum(normNEG,2);
                minlen = min(length(trajplot{1,1}{1}(:,5)), length(combnormNEG));
                [r, p] = corrcoef(trajplot{1}{1}(1:minlen,5), combnormNEG(1:minlen,1), 'rows', 'pairwise'); %don't use rows if either has a nan
                CorrCoefcombNEG(CorrCoefcombNEGCOUNT,1) =  r(2,1);
            end
            
            if poscelsinset > 0 && negcelsinset > 0%(poscelsinset > 1 & negcelsinset > 0) | (poscelsinset > 0 & negcelsinset > 1) ; %length(pairdata(1,i).CA1posidxs(:,1))>1 | length(pairdata(1,i).CA1negidxs(:,1))>1; %if there are more than 1 positive and negative cells for this set.
                CorrCoefcombPOSNEGCOUNT = CorrCoefcombPOSNEGCOUNT +1;
                combnormPOSarith = nansum(normPOS,2); %need to redefine in case only 1 cell
                combnormNEGarith = nansum(normNEG,2);
                minlen = min(length(combnormPOSarith), length(combnormNEGarith));
                combnormPOSNEGarith = nansum([combnormPOSarith(1:minlen) -combnormNEGarith(1:minlen)],2); %poscells - negcells element-wise
                minlen = min(length(trajplot{1,1}{1}(:,5)), length(combnormPOSNEGarith));
                [r, p] = corrcoef(trajplot{1}{1}(1:minlen,5), combnormPOSNEGarith(1:minlen,1), 'rows', 'pairwise'); %don't use rows if either has a nan
                CorrCoefcombPOSNEG(CorrCoefcombPOSNEGCOUNT,1) =  r(2,1);
            end
            
            
            
            %compute corr coef for each pair in current set.._____________________________________________________________________________________
            %--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            clear CorrCoefstructSET; pairinset = 0; %clear this for each set
            if length(trajplot(1,:)) >1; %if there are any ca1 cells
                for cu = 2:(length(trajplot(1,:))); %all ca1 cells ~ # of pairs w pfc
                    paircounter = paircounter +1; %count each pair for storing
                    pairinset = pairinset +1;
                    tru = 1; %for each trajectory.. update.. now that i've concatenated all trajs, there will only be one array (one long trajectory)
                    clear r p minlen;
                    minlen = min(length(trajplot{1,1}{tru}(:,5)), length(trajplot{1,cu}{tru}(:,5)));
                    [r, p] = corrcoef(trajplot{1, 1}{tru}(1:minlen,5), trajplot{1, cu}{tru}(1:minlen,5), 'rows', 'pairwise'); %don't use rows if either has a nan
                    CorrCoefstruct{paircounter,1}.r(tru,1) = r(2,1); %store corr coefs
                    CorrCoefstruct{paircounter,1}.p(tru,1) = p(2,1); %store corr coefs p vals
                    CorrCoefstructSET(pairinset,1) = r(2,1); %store the corr coefs of just this set for plotting purposes
                    
                    %                     end %each traj
                    
                    if cu-1 <= poscelsinset; %length(pairdata(1,i).CA1posidxs(:,1)); %if pos corr ca1 cell. cu-1 so that the count starts at the beginning of the pos indices list
                        CorrCoefstruct{paircounter, 1}.posneg = 'pos';
                        CorrCoefstruct{paircounter, 1}.pfcindex(1,:) = pairdata(1,i).PFCidx(1,:);
                        CorrCoefstruct{paircounter, 1}.ca1index(1,:) = pairdata(1,i).CA1posidxs(cu-1,:);
                        CorrCoefstruct{paircounter, 1}.sigbetas = pairdata(1,i).sigbetaspos(cu-1,:);
                        sigbetasscatter(paircounter, 1) = pairdata(1,i).sigbetaspos(cu-1,:);
                        spatcorrvalsscatter(paircounter, 1) = r(2,1);
                        pospaircounter = pospaircounter +1;
                        CorrCoefstructPOS{pospaircounter, 1} = CorrCoefstruct{paircounter, 1}; %make copy of pos struct for pos store
                        
                    else cu-1 <= poscelsinset + negcelsinset; %if neg corr ca1 cell
                        CorrCoefstruct{paircounter, 1}.posneg = 'neg';
                        CorrCoefstruct{paircounter, 1}.pfcindex(1,:) = pairdata(1,i).PFCidx(1,:);
                        CorrCoefstruct{paircounter, 1}.ca1index(1,:) = pairdata(1,i).CA1negidxs(cu-poscelsinset-1,:);
                        CorrCoefstruct{paircounter, 1}.sigbetas = pairdata(1,i).sigbetasneg(cu-poscelsinset-1,:);
                        sigbetasscatter(paircounter, 1) = pairdata(1,i).sigbetasneg(cu-poscelsinset-1,:);
                        spatcorrvalsscatter(paircounter, 1) = r(2,1);
                        negpaircounter = negpaircounter +1;
                        CorrCoefstructNEG{negpaircounter, 1} = CorrCoefstruct{paircounter, 1}; %make copy of neg struct for neg store
                        
                    end
                end %eac ca1 cell in set
            end %if any ca1 cells
            
            % %make massive shuff struct for EACH animal with corrcoefs for each pfc cell with all ca1 cells from within the same animal_____________________
            %--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            %             if length(trajplot(1,:)) >1; %if there are any ca1 cells
            if pairdata(1,i).PFCidx(1,1) == 1%HPa
                for mo = 1:length(HPaCA1structShortListCombEp); %loop over all the HP CA1 cells for the shuffled struct XXX loop over unique(HPaCA1structShortList) instead
                    clear hpaca1shuffEP
                    %loop over CA1 cells that are not significantly correlated with the current PFC cell
                    clear currca1cell
                    currca1cell = num2str(HPaCA1structShortListCombEp(mo)); %this line and the next are used to limit the shuff struct to ca1 cells from the same day of the pfc cell.
                    if str2num(currca1cell(2)) == pairdata(1,i).PFCidx(1,2);
                        if isempty(find(HPaCA1structShortListCombEp(mo) == CA1ShortList)); %if this ca1 cell is not on the pos ca1 list.
                            if isempty(find(HPaCA1structShortListCombEp(mo)  == CA1ShortListb)); %and if this ca1 cell is not on the neg ca1 list
                                mocnt = mocnt +1; %across all anims
                                hpacnt = hpacnt +1; %within an
                                %combine trajectories across epochs for
                                %current CA1 cell.
                                %                                     hpamatches = find(str2num(sprintf('%d',[HPaCA1struct{mo,1}.index([1 3 4])])) == HPaCA1structShortList); %find all epochs
                                hpamatches = find(HPaCA1structShortListCombEp(mo) == HPaCA1structShortList(:)); %find all epochs
                                
                                for hpau = 1:length(hpamatches); %collect traj data for each epoch
                                    hpaca1shufftrajs{hpau,1} = HPaCA1struct{hpamatches(hpau)}.trajdata;
                                end
                                celltraj1 = []; celltraj2 = []; celltraj3 = []; celltraj4 = [];
                                for epcell = 1:length(hpaca1shufftrajs(:,1));
                                    %find the min length of the trajectories to be averaged
                                    if isempty(celltraj1); %if it's the first epoch
                                        minlen1 = length(hpaca1shufftrajs{epcell,1}{1,1}(:,5));
                                        minlen2 = length(hpaca1shufftrajs{epcell,1}{1,2}(:,5));
                                        minlen3 = length(hpaca1shufftrajs{epcell,1}{1,3}(:,5));
                                        minlen4 = length(hpaca1shufftrajs{epcell,1}{1,4}(:,5));
                                        celltraj1 = [hpaca1shufftrajs{epcell,1}{1,1}(1:minlen1,5)]; %horizontally concatenate trajectory 1
                                        celltraj2 = [hpaca1shufftrajs{epcell,1}{1,2}(1:minlen2,5)]; %horizontally concatenate trajectory 2
                                        celltraj3 = [hpaca1shufftrajs{epcell,1}{1,3}(1:minlen3,5)]; %horizontally concatenate trajectory 3
                                        celltraj4 = [hpaca1shufftrajs{epcell,1}{1,4}(1:minlen4,5)]; %horizontally concatenate trajectory 4
                                    else%all the rest of the epochs
                                        minlen1 = min(length(celltraj1), length(hpaca1shufftrajs{epcell,1}{1,1}(:,5)));
                                        minlen2 = min(length(celltraj2), length(hpaca1shufftrajs{epcell,1}{1,2}(:,5)));
                                        minlen3 = min(length(celltraj3), length(hpaca1shufftrajs{epcell,1}{1,3}(:,5)));
                                        minlen4 = min(length(celltraj4), length(hpaca1shufftrajs{epcell,1}{1,4}(:,5)));
                                        celltraj1 = [celltraj1(1:minlen1,:) hpaca1shufftrajs{epcell,1}{1,1}(1:minlen1,5)]; %horizontally concatenate trajectory 1
                                        celltraj2 = [celltraj2(1:minlen2,:) hpaca1shufftrajs{epcell,1}{1,2}(1:minlen2,5)]; %horizontally concatenate trajectory 2
                                        celltraj3 = [celltraj3(1:minlen3,:) hpaca1shufftrajs{epcell,1}{1,3}(1:minlen3,5)]; %horizontally concatenate trajectory 3
                                        celltraj4 = [celltraj4(1:minlen4,:) hpaca1shufftrajs{epcell,1}{1,4}(1:minlen4,5)]; %horizontally concatenate trajectory 4
                                    end
                                end
                                
                                hpaca1shuffEP{1, 1}{1,1}(:,5) = nanmean(celltraj1,2);
                                hpaca1shuffEP{1, 1}{1,2}(:,5) = nanmean(celltraj2,2);
                                hpaca1shuffEP{1, 1}{1,3}(:,5) = nanmean(celltraj3,2);
                                hpaca1shuffEP{1, 1}{1,4}(:,5) = nanmean(celltraj4,2);
                                
                                %concatenate trajectories
                                combtraj = [];
                                clear tmptrajplotcombtraj;
                                for eachtraj = 1:length(hpaca1shuffEP{1, 1}); %all trajs
                                    combtraj = [combtraj; hpaca1shuffEP{1, 1}{eachtraj}];
                                end
                                tmptrajplotcombtraj{1}{1} = combtraj;
                                hpaca1shuffEP = tmptrajplotcombtraj;
                                
                                
                                for tru = 1; %for each trajectory. now just 1 concatenated
                                    clear r p minlen;
                                    minlen = min(length(trajplot{1}{tru}(:,5)), length(hpaca1shuffEP{1, 1}{1,tru}(:,5)));
                                    [r, p] = corrcoef(trajplot{1}{tru}(1:minlen,5), hpaca1shuffEP{1, 1}{1,tru}(1:minlen,5), 'rows', 'pairwise'); %don't use rows if either has a nan
                                    CorrCoefstructShuff{mocnt,1}.r(tru,1) = r(2,1); %store corr coefs
                                    CorrCoefstructShuff{mocnt,1}.p(tru,1) = p(2,1); %store corr coefs p vals
                                    CorrCoefstructShuffHPa{hpacnt,1}.r(tru,1) = r(2,1); %store corr coefs
                                    CorrCoefstructShuffHPa{hpacnt,1}.p(tru,1) = p(2,1); %store corr coefs p vals
                                end
                                CorrCoefstructShuff{mocnt, 1}.pfcindex(1,:) = pairdata(1,i).PFCidx(1,:);
                                CorrCoefstructShuff{mocnt, 1}.ca1index(1,:) = HPaCA1struct{hpamatches(1)}.index; %taking the index of the first epoch used
                                CorrCoefstructShuffHPa{hpacnt, 1}.pfcindex(1,:) = pairdata(1,i).PFCidx(1,:);
                                CorrCoefstructShuffHPa{hpacnt, 1}.ca1index(1,:) = HPaCA1struct{hpamatches(1)}.index;
                            end
                        end
                    end
                end
                
            elseif pairdata(1,i).PFCidx(1,1) == 2%HPb
                for mo = 1:length(HPbCA1structShortListCombEp); %loop over all the HP CA1 cells for the shuffled struct XXX loop over unique(HPbCA1structShortList) instead
                    clear hpbca1shuffEP
                    %loop over CA1 cells that are no significantly correlated with the current PFC cell
                    clear currca1cell
                    currca1cell = num2str(HPbCA1structShortListCombEp(mo)); %this line and the next are used to limit the shuff struct to ca1 cells from the same day of the pfc cell.
                    if str2num(currca1cell(2)) == pairdata(1,i).PFCidx(1,2);
                        if isempty(find(HPbCA1structShortListCombEp(mo) == CA1ShortList)); %if this ca1 cell is not on the pos ca1 list.
                            if isempty(find(HPbCA1structShortListCombEp(mo)  == CA1ShortListb)) ; %and if this ca1 cell is not on the neg ca1 list
                                mocnt = mocnt +1; %across all anims
                                hpbcnt = hpbcnt +1; %within an
                                %combine trajectories across epochs for
                                %current CA1 cell.
                                %                                     hpbmatches = find(str2num(sprintf('%d',[hpbCA1struct{mo,1}.index([1 3 4])])) == hpbCA1structShortList); %find all epochs
                                hpbmatches = find(HPbCA1structShortListCombEp(mo) == HPbCA1structShortList(:)); %find all epochs
                                
                                for hpbu = 1:length(hpbmatches); %collect traj data for each epoch
                                    hpbca1shufftrajs{hpbu,1} = HPbCA1struct{hpbmatches(hpbu)}.trajdata;
                                end
                                celltraj1 = []; celltraj2 = []; celltraj3 = []; celltraj4 = [];
                                for epcell = 1:length(hpbca1shufftrajs(:,1));
                                    %find the min length of the trajectories to be averaged
                                    if isempty(celltraj1); %if it's the first epoch
                                        minlen1 = length(hpbca1shufftrajs{epcell,1}{1,1}(:,5));
                                        minlen2 = length(hpbca1shufftrajs{epcell,1}{1,2}(:,5));
                                        minlen3 = length(hpbca1shufftrajs{epcell,1}{1,3}(:,5));
                                        minlen4 = length(hpbca1shufftrajs{epcell,1}{1,4}(:,5));
                                        celltraj1 = [hpbca1shufftrajs{epcell,1}{1,1}(1:minlen1,5)]; %horizontally concatenate trajectory 1
                                        celltraj2 = [hpbca1shufftrajs{epcell,1}{1,2}(1:minlen2,5)]; %horizontally concatenate trajectory 2
                                        celltraj3 = [hpbca1shufftrajs{epcell,1}{1,3}(1:minlen3,5)]; %horizontally concatenate trajectory 3
                                        celltraj4 = [hpbca1shufftrajs{epcell,1}{1,4}(1:minlen4,5)]; %horizontally concatenate trajectory 4
                                    else%all the rest of the epochs
                                        minlen1 = min(length(celltraj1), length(hpbca1shufftrajs{epcell,1}{1,1}(:,5)));
                                        minlen2 = min(length(celltraj2), length(hpbca1shufftrajs{epcell,1}{1,2}(:,5)));
                                        minlen3 = min(length(celltraj3), length(hpbca1shufftrajs{epcell,1}{1,3}(:,5)));
                                        minlen4 = min(length(celltraj4), length(hpbca1shufftrajs{epcell,1}{1,4}(:,5)));
                                        celltraj1 = [celltraj1(1:minlen1,:) hpbca1shufftrajs{epcell,1}{1,1}(1:minlen1,5)]; %horizontally concatenate trajectory 1
                                        celltraj2 = [celltraj2(1:minlen2,:) hpbca1shufftrajs{epcell,1}{1,2}(1:minlen2,5)]; %horizontally concatenate trajectory 2
                                        celltraj3 = [celltraj3(1:minlen3,:) hpbca1shufftrajs{epcell,1}{1,3}(1:minlen3,5)]; %horizontally concatenate trajectory 3
                                        celltraj4 = [celltraj4(1:minlen4,:) hpbca1shufftrajs{epcell,1}{1,4}(1:minlen4,5)]; %horizontally concatenate trajectory 4
                                    end
                                end
                                
                                hpbca1shuffEP{1, 1}{1,1}(:,5) = nanmean(celltraj1,2);
                                hpbca1shuffEP{1, 1}{1,2}(:,5) = nanmean(celltraj2,2);
                                hpbca1shuffEP{1, 1}{1,3}(:,5) = nanmean(celltraj3,2);
                                hpbca1shuffEP{1, 1}{1,4}(:,5) = nanmean(celltraj4,2);
                                
                                
                                %concatenate trajectories
                                combtraj = [];
                                clear tmptrajplotcombtraj;
                                for eachtraj = 1:length(hpbca1shuffEP{1, 1}); %all trajs
                                    combtraj = [combtraj; hpbca1shuffEP{1, 1}{eachtraj}];
                                end
                                tmptrajplotcombtraj{1}{1} = combtraj;
                                hpbca1shuffEP = tmptrajplotcombtraj;
                                
                                
                                for tru = 1; %for each trajectory
                                    clear r p minlen;
                                    minlen = min(length(trajplot{1}{tru}(:,5)), length(hpbca1shuffEP{1, 1}{1,tru}(:,5)));
                                    [r, p] = corrcoef(trajplot{1}{tru}(1:minlen,5), hpbca1shuffEP{1, 1}{1,tru}(1:minlen,5), 'rows', 'pairwise'); %don't use rows if either has a nan
                                    CorrCoefstructShuff{mocnt,1}.r(tru,1) = r(2,1); %store corr coefs
                                    CorrCoefstructShuff{mocnt,1}.p(tru,1) = p(2,1); %store corr coefs p vals
                                    CorrCoefstructShuffHPb{hpbcnt,1}.r(tru,1) = r(2,1); %store corr coefs
                                    CorrCoefstructShuffHPb{hpbcnt,1}.p(tru,1) = p(2,1); %store corr coefs p vals
                                end
                                CorrCoefstructShuff{mocnt, 1}.pfcindex(1,:) = pairdata(1,i).PFCidx(1,:);
                                CorrCoefstructShuff{mocnt, 1}.ca1index(1,:) = HPbCA1struct{hpbmatches(1)}.index; %taking the index of the first epoch used
                                CorrCoefstructShuffHPb{hpbcnt, 1}.pfcindex(1,:) = pairdata(1,i).PFCidx(1,:);
                                CorrCoefstructShuffHPb{hpbcnt, 1}.ca1index(1,:) = HPbCA1struct{hpbmatches(1)}.index;
                            end
                        end
                    end
                end
                
                
            elseif pairdata(1,i).PFCidx(1,1) == 3%HPc
                for mo = 1:length(HPcCA1structShortListCombEp); %loop over all the HP CA1 cells for the shuffled struct XXX loop over unique(HPaCA1structShortList) instead
                    clear hpcca1shuffEP
                    %loop over CA1 cells that are no significantly correlated with the current PFC cell
                    clear currca1cell
                    currca1cell = num2str(HPcCA1structShortListCombEp(mo)); %this line and the next are used to limit the shuff struct to ca1 cells from the same day of the pfc cell.
                    if str2num(currca1cell(2)) == pairdata(1,i).PFCidx(1,2);
                        if isempty(find(HPcCA1structShortListCombEp(mo) == CA1ShortList)); %if this ca1 cell is not on the pos ca1 list.
                            if isempty(find(HPcCA1structShortListCombEp(mo)  == CA1ShortListb)) ; %and if this ca1 cell is not on the neg ca1 list
                                mocnt = mocnt +1; %across all anims
                                hpccnt = hpccnt +1; %within an
                                %combine trajectories across epochs for
                                %current CA1 cell.
                                %                                     hpamatches = find(str2num(sprintf('%d',[HPaCA1struct{mo,1}.index([1 3 4])])) == HPaCA1structShortList); %find all epochs
                                hpcmatches = find(HPcCA1structShortListCombEp(mo) == HPcCA1structShortList(:)); %find all epochs
                                
                                for hpcu = 1:length(hpcmatches); %collect traj data for each epoch
                                    hpcca1shufftrajs{hpcu,1} = HPcCA1struct{hpcmatches(hpcu)}.trajdata;
                                end
                                celltraj1 = []; celltraj2 = []; celltraj3 = []; celltraj4 = [];
                                for epcell = 1:length(hpcca1shufftrajs(:,1));
                                    %find the min length of the trajectories to be averaged
                                    if isempty(celltraj1); %if it's the first epoch
                                        minlen1 = length(hpcca1shufftrajs{epcell,1}{1,1}(:,5));
                                        minlen2 = length(hpcca1shufftrajs{epcell,1}{1,2}(:,5));
                                        minlen3 = length(hpcca1shufftrajs{epcell,1}{1,3}(:,5));
                                        minlen4 = length(hpcca1shufftrajs{epcell,1}{1,4}(:,5));
                                        celltraj1 = [hpcca1shufftrajs{epcell,1}{1,1}(1:minlen1,5)]; %horizontally concatenate trajectory 1
                                        celltraj2 = [hpcca1shufftrajs{epcell,1}{1,2}(1:minlen2,5)]; %horizontally concatenate trajectory 2
                                        celltraj3 = [hpcca1shufftrajs{epcell,1}{1,3}(1:minlen3,5)]; %horizontally concatenate trajectory 3
                                        celltraj4 = [hpcca1shufftrajs{epcell,1}{1,4}(1:minlen4,5)]; %horizontally concatenate trajectory 4
                                    else%all the rest of the epochs
                                        minlen1 = min(length(celltraj1), length(hpcca1shufftrajs{epcell,1}{1,1}(:,5)));
                                        minlen2 = min(length(celltraj2), length(hpcca1shufftrajs{epcell,1}{1,2}(:,5)));
                                        minlen3 = min(length(celltraj3), length(hpcca1shufftrajs{epcell,1}{1,3}(:,5)));
                                        minlen4 = min(length(celltraj4), length(hpcca1shufftrajs{epcell,1}{1,4}(:,5)));
                                        celltraj1 = [celltraj1(1:minlen1,:) hpcca1shufftrajs{epcell,1}{1,1}(1:minlen1,5)]; %horizontally concatenate trajectory 1
                                        celltraj2 = [celltraj2(1:minlen2,:) hpcca1shufftrajs{epcell,1}{1,2}(1:minlen2,5)]; %horizontally concatenate trajectory 2
                                        celltraj3 = [celltraj3(1:minlen3,:) hpcca1shufftrajs{epcell,1}{1,3}(1:minlen3,5)]; %horizontally concatenate trajectory 3
                                        celltraj4 = [celltraj4(1:minlen4,:) hpcca1shufftrajs{epcell,1}{1,4}(1:minlen4,5)]; %horizontally concatenate trajectory 4
                                    end
                                end
                                
                                hpcca1shuffEP{1, 1}{1,1}(:,5) = nanmean(celltraj1,2);
                                hpcca1shuffEP{1, 1}{1,2}(:,5) = nanmean(celltraj2,2);
                                hpcca1shuffEP{1, 1}{1,3}(:,5) = nanmean(celltraj3,2);
                                hpcca1shuffEP{1, 1}{1,4}(:,5) = nanmean(celltraj4,2);
                                
                                
                                %concatenate trajectories
                                combtraj = [];
                                clear tmptrajplotcombtraj;
                                for eachtraj = 1:length(hpcca1shuffEP{1, 1}); %all trajs
                                    combtraj = [combtraj; hpcca1shuffEP{1, 1}{eachtraj}];
                                end
                                tmptrajplotcombtraj{1}{1} = combtraj;
                                hpcca1shuffEP = tmptrajplotcombtraj;
                                
                                
                                for tru = 1; %for each trajectory
                                    clear r p minlen;
                                    minlen = min(length(trajplot{1}{tru}(:,5)), length(hpcca1shuffEP{1, 1}{1,tru}(:,5)));
                                    [r, p] = corrcoef(trajplot{1}{tru}(1:minlen,5), hpcca1shuffEP{1, 1}{1,tru}(1:minlen,5), 'rows', 'pairwise'); %don't use rows if either has a nan
                                    CorrCoefstructShuff{mocnt,1}.r(tru,1) = r(2,1); %store corr coefs
                                    CorrCoefstructShuff{mocnt,1}.p(tru,1) = p(2,1); %store corr coefs p vals
                                    CorrCoefstructShuffHPc{hpccnt,1}.r(tru,1) = r(2,1); %store corr coefs
                                    CorrCoefstructShuffHPc{hpccnt,1}.p(tru,1) = p(2,1); %store corr coefs p vals
                                end
                                CorrCoefstructShuff{mocnt, 1}.pfcindex(1,:) = pairdata(1,i).PFCidx(1,:);
                                CorrCoefstructShuff{mocnt, 1}.ca1index(1,:) = HPcCA1struct{hpcmatches(1)}.index; %taking the index of the first epoch used
                                CorrCoefstructShuffHPc{hpccnt, 1}.pfcindex(1,:) = pairdata(1,i).PFCidx(1,:);
                                CorrCoefstructShuffHPc{hpccnt, 1}.ca1index(1,:) = HPcCA1struct{hpcmatches(1)}.index;
                            end
                        end
                    end
                end
                
                
            elseif pairdata(1,i).PFCidx(1,1) == 4%Nadal
                for mo = 1:length(NdlCA1structShortListCombEp); %loop over all the HP CA1 cells for the shuffled struct XXX loop over unique(NdlCA1structShortList) instead
                    clear ndlca1shuffEP
                    %loop over CA1 cells that are no significantly correlated with the current PFC cell
                    clear currca1cell daymatchesndl
                    %bc nadal has 1 or 2 digit day vals.. need to adjust how i check the day match..
                    daymatchesndl = find(NdlCA1structShortListCombEp(mo) == NdlCA1structShortList,1); %find the first match of the shortlist in the long list to then use to get the unconcatented indices from the data struct below
                    currca1cell = num2str(NdlCA1struct{daymatchesndl}.index(1)); %this line and the next are used to limit the shuff struct to ca1 cells from the same day of the pfc cell.
                    if str2num(currca1cell) == pairdata(1,i).PFCidx(1,2);
                        if isempty(find(NdlCA1structShortListCombEp(mo) == CA1ShortList)); %if this ca1 cell is not on the pos ca1 list.
                            if isempty(find(NdlCA1structShortListCombEp(mo)  == CA1ShortListb)) ; %and if this ca1 cell is not on the neg ca1 list
                                mocnt = mocnt +1; %across all anims
                                ndlcnt = ndlcnt +1; %within an
                                %combine trajectories across epochs for
                                %current CA1 cell.
                                %                                     hpamatches = find(str2num(sprintf('%d',[HPaCA1struct{mo,1}.index([1 3 4])])) == HPaCA1structShortList); %find all epochs
                                ndlmatches = find(NdlCA1structShortListCombEp(mo) == NdlCA1structShortList(:)); %find all epochs
                                
                                for ndlu = 1:length(ndlmatches); %collect traj data for each epoch
                                    ndlca1shufftrajs{ndlu,1} = NdlCA1struct{ndlmatches(ndlu)}.trajdata;
                                end
                                celltraj1 = []; celltraj2 = []; celltraj3 = []; celltraj4 = [];
                                for epcell = 1:length(ndlca1shufftrajs(:,1));
                                    %find the min length of the trajectories to be averaged
                                    if isempty(celltraj1); %if it's the first epoch
                                        minlen1 = length(ndlca1shufftrajs{epcell,1}{1,1}(:,5));
                                        minlen2 = length(ndlca1shufftrajs{epcell,1}{1,2}(:,5));
                                        minlen3 = length(ndlca1shufftrajs{epcell,1}{1,3}(:,5));
                                        minlen4 = length(ndlca1shufftrajs{epcell,1}{1,4}(:,5));
                                        celltraj1 = [ndlca1shufftrajs{epcell,1}{1,1}(1:minlen1,5)]; %horizontally concatenate trajectory 1
                                        celltraj2 = [ndlca1shufftrajs{epcell,1}{1,2}(1:minlen2,5)]; %horizontally concatenate trajectory 2
                                        celltraj3 = [ndlca1shufftrajs{epcell,1}{1,3}(1:minlen3,5)]; %horizontally concatenate trajectory 3
                                        celltraj4 = [ndlca1shufftrajs{epcell,1}{1,4}(1:minlen4,5)]; %horizontally concatenate trajectory 4
                                    else%all the rest of the epochs
                                        minlen1 = min(length(celltraj1), length(ndlca1shufftrajs{epcell,1}{1,1}(:,5)));
                                        minlen2 = min(length(celltraj2), length(ndlca1shufftrajs{epcell,1}{1,2}(:,5)));
                                        minlen3 = min(length(celltraj3), length(ndlca1shufftrajs{epcell,1}{1,3}(:,5)));
                                        minlen4 = min(length(celltraj4), length(ndlca1shufftrajs{epcell,1}{1,4}(:,5)));
                                        celltraj1 = [celltraj1(1:minlen1,:) ndlca1shufftrajs{epcell,1}{1,1}(1:minlen1,5)]; %horizontally concatenate trajectory 1
                                        celltraj2 = [celltraj2(1:minlen2,:) ndlca1shufftrajs{epcell,1}{1,2}(1:minlen2,5)]; %horizontally concatenate trajectory 2
                                        celltraj3 = [celltraj3(1:minlen3,:) ndlca1shufftrajs{epcell,1}{1,3}(1:minlen3,5)]; %horizontally concatenate trajectory 3
                                        celltraj4 = [celltraj4(1:minlen4,:) ndlca1shufftrajs{epcell,1}{1,4}(1:minlen4,5)]; %horizontally concatenate trajectory 4
                                    end
                                end
                                
                                ndlca1shuffEP{1, 1}{1,1}(:,5) = nanmean(celltraj1,2);
                                ndlca1shuffEP{1, 1}{1,2}(:,5) = nanmean(celltraj2,2);
                                ndlca1shuffEP{1, 1}{1,3}(:,5) = nanmean(celltraj3,2);
                                ndlca1shuffEP{1, 1}{1,4}(:,5) = nanmean(celltraj4,2);
                                
                                
                                %concatenate trajectories
                                combtraj = [];
                                clear tmptrajplotcombtraj;
                                for eachtraj = 1:length(ndlca1shuffEP{1, 1}); %all trajs
                                    combtraj = [combtraj; ndlca1shuffEP{1, 1}{eachtraj}];
                                end
                                tmptrajplotcombtraj{1}{1} = combtraj;
                                ndlca1shuffEP = tmptrajplotcombtraj;
                                
                                
                                for tru = 1; %for each trajectory
                                    clear r p minlen;
                                    minlen = min(length(trajplot{1}{tru}(:,5)), length(ndlca1shuffEP{1, 1}{1,tru}(:,5)));
                                    [r, p] = corrcoef(trajplot{1}{tru}(1:minlen,5), ndlca1shuffEP{1, 1}{1,tru}(1:minlen,5), 'rows', 'pairwise'); %don't use rows if either has a nan
                                    CorrCoefstructShuff{mocnt,1}.r(tru,1) = r(2,1); %store corr coefs
                                    CorrCoefstructShuff{mocnt,1}.p(tru,1) = p(2,1); %store corr coefs p vals
                                    CorrCoefstructShuffNdl{ndlcnt,1}.r(tru,1) = r(2,1); %store corr coefs
                                    CorrCoefstructShuffNdl{ndlcnt,1}.p(tru,1) = p(2,1); %store corr coefs p vals
                                end
                                CorrCoefstructShuff{mocnt, 1}.pfcindex(1,:) = pairdata(1,i).PFCidx(1,:);
                                CorrCoefstructShuff{mocnt, 1}.ca1index(1,:) = NdlCA1struct{ndlmatches(1)}.index; %taking the index of the first epoch used
                                CorrCoefstructShuffNdl{ndlcnt, 1}.pfcindex(1,:) = pairdata(1,i).PFCidx(1,:);
                                CorrCoefstructShuffNdl{ndlcnt, 1}.ca1index(1,:) = NdlCA1struct{ndlmatches(1)}.index;
                            end
                        end
                    end
                end
            end
            %             end
            
            
            %--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            %plot concatenated trajectories of each ca1 on the pfc with the
            %combined cells on the bottom
            
            if plotPFCCA1Trajs ==1;
                if length(trajplot(1,:)) >1; %if there are any ca1 cells
                    
                    %but first collect the noise corr vals.. need to reorder them
                    %to fit the pos then neg order of the trajdata
                    clear  sigrset
                    %                     sigrpos = (pairdata(1,i).rsig([find(pairdata(1,i).rsig(:,1)>0,length(pairdata(1,i).rsig))],1))'; %collect the positive noise corr vals for this set
                    %                     sigrneg = (pairdata(1,i).rsig([find(pairdata(1,i).rsig(:,1)<0,length(pairdata(1,i).rsig))],1))';  %collect the negative noise corr vals for this set
                    sigrset = [pairdata(1,i).sigbetaspos; pairdata(1,i).sigbetasneg]; %vert stack so that pos first, as ive collected the traj data
                    complot = 0;
                    if exist('combnormPOS') | exist('combnormNEG') | exist('combnormPOSNEGarith'); %if there's a combined struct
                        complot = 1;
                    end
                    for ru = 2:length(trajplot)+complot; %all the ca1 cells in this set;; USING data averages across epochs AND concatenated trajectories
                        subplot(length(trajplot)-1+complot,1,ru-1); hold on; %make fig that fits all cells in set
                        if ru-1 <= poscelsinset; %length(pairdata(1,i).CA1posidxs(:,1)); %pos cell
                            plot(trajplot{1}{1}(:,5)./max(trajplot{1}{1}(:,5)), '.-','Color', 'k','Linewidth',trajline);
                            plot(trajplot{ru}{1}(:,5)./max(trajplot{ru}{1}(:,5)), '.-','Color',[0 .7 .93],'Linewidth',trajline);
                            tit = {sprintf('%s PFC(%d %d %d) CA1(%d %d %d)  RipCor(%0.5f) SpCor(%0.5f)', animals{pairdata(1,i).CA1sigidxs(ru-1,1)}, pairdata(1,i).PFCidx(1,2), pairdata(1,i).PFCidx(1,3), pairdata(1,i).PFCidx(1,4),...
                                pairdata(1,i).CA1sigidxs(ru-1,2), pairdata(1,i).CA1sigidxs(ru-1,3), pairdata(1,i).CA1sigidxs(ru-1,4), sigrset(ru-1,1),CorrCoefstructSET(ru-1,1))}; %label w pfc indx, ca1 idx, spcorr
                        elseif ru-1 <= poscelsinset + negcelsinset; %length(pairdata(1,i).CA1sigidxs(:,1)); %neg cell
                            plot(trajplot{1}{1}(:,5)./max(trajplot{1}{1}(:,5)), '.-','Color', 'k','Linewidth',trajline);
                            plot(trajplot{ru}{1}(:,5)./max(trajplot{ru}{1}(:,5)),'.-','Color',[.9 .9 .4], 'Linewidth',trajline);
                            tit = {sprintf('%s PFC(%d %d %d) CA1(%d %d %d)  RipCor(%0.5f) SpCor(%0.5f)', animals{pairdata(1,i).CA1sigidxs(ru-1,1)}, pairdata(1,i).PFCidx(1,2), pairdata(1,i).PFCidx(1,3), pairdata(1,i).PFCidx(1,4),...
                                pairdata(1,i).CA1sigidxs(ru-1,2), pairdata(1,i).CA1sigidxs(ru-1,3), pairdata(1,i).CA1sigidxs(ru-1,4), sigrset(ru-1,1),CorrCoefstructSET(ru-1,1))}; %label w pfc indx, ca1 idx, spcorr
                        else %comb plot
                            titcom = [];
                            plot((trajplot{1}{1}(:,5)./max(trajplot{1}{1}(:,5))), '.-','Color', 'k','Linewidth',trajline); %normalized pfc cell
                            if exist('combnormPOS');
                                plot(combnormPOS,'.-','Color',[.19 .39 .70], 'Linewidth',trajline);
                                titcom = [titcom sprintf('POSC(%0.5f)  ', CorrCoefcombPOS(CorrCoefcombPOSCOUNT,1))];
                            end
                            if exist('combnormNEG');
                                plot(combnormNEG,'.-','Color',[.93 .79 0], 'Linewidth',trajline);
                                titcom = [titcom sprintf('NEGC(%0.5f)  ', CorrCoefcombNEG(CorrCoefcombNEGCOUNT,1))];
                            end
                            if exist('combnormPOSNEGarith');
                                plot(combnormPOSNEGarith,'.-','Color',[0 1 0], 'Linewidth',trajline);
                                titcom = [titcom sprintf('POSNEGC(%0.5f)  ', CorrCoefcombPOSNEG(CorrCoefcombPOSNEGCOUNT,1))];
                            end
                            title(titcom, 'Fontsize', titlesize,'FontName',fonttype)
                        end
                        set(gca,'xtick',[])
                        set(gca,'xticklabel',[])
                        set(gca,'ytick',[])
                        set(gca,'yticklabel',[])
                        if ru <= length(trajplot); %if its not the combined stuff
                            title(tit, 'Fontsize', titlesize,'FontName',fonttype)
                        end
                        yLimits = get(gca,'YLim');  %# Get the range of the y axis
                        trj1 = length(trajplotEpAvg{1}{1}(:,5));
                        trj2 = trj1 + length(trajplotEpAvg{1}{2}(:,5));
                        trj3 = trj2 + length(trajplotEpAvg{1}{3}(:,5));
                        trj4 = trj3 + length(trajplotEpAvg{1}{4}(:,5));
                        line('XData',[trj1 trj1], 'YData', yLimits, 'LineStyle', '-', 'LineWidth', trajline, 'Color',[.8 .8 .8]);
                        line('XData',[trj2 trj2], 'YData', yLimits, 'LineStyle', '-', 'LineWidth', trajline, 'Color',[.8 .8 .8]);
                        line('XData',[trj3 trj3], 'YData', yLimits, 'LineStyle', '-', 'LineWidth', trajline, 'Color',[.8 .8 .8]);
                        line('XData',[trj4 trj4], 'YData', yLimits, 'LineStyle', '-', 'LineWidth', trajline, 'Color',[.8 .8 .8]);
                        %                         set(gca,'FontSize',2,'fontWeight','bold')
                        %collect the spCorr vals along side the noise corr
                        %vals
                        if ru <= length(trajplot); %if its not the combined stuff
                            AllCorrRip_Sp = [AllCorrRip_Sp; sigrset(ru-1,1) CorrCoefstructSET(ru-1,1)];
                        end
                    end
                    
                    figfile = [figdir,savefigfilename,'/','Trajs',animals{pairdata(1,i).PFCidx(1,1)}, num2str(i)];
                    if savefigs==1
                        print('-dpdf', figfile, '-r300');
                    end
                    if ~cyclemaps == 0
                        keyboard;
                    end
                    close
                end %if there are any ca1 cells in this
            end %for plot maps
            %dun touch these
            
            
            
            
            
            %-----------------------------------------------------------------------------------
            %this won't work any more in the current state as I've changed trajplot to be concatenated trajectories and nanmeaned over epochs..
            %although I can recover this by using the old trajplot (trajplotEpAvg) and altering this code to work epoch by epoch.. update.. now
            %im using the concatenated version of the trajectories..%holding off on the maps for now. update... trying the mpapsagain with combined
            %epochs.. days w diff 2 tracks will be messed up but this is just for visual for now
            
            if plotPFCCA1maps ==1;
                if savefigs == 1;
                    mkdir(figdir,[savefigfilename '_maps'])
                end
                trajindexshortlist = [PFCShortList(i); CA1ShortList; CA1ShortListb]; %The indices of the cells in set in order PFC posCA1 negCA1
                %                 load('/mnt/data25/sjadhav/HPExpt/ProcessedDataDR/AllAn_PFCCA1_ripplepos_DR_rippleposindex_vel4tet2filt');
                load(loadrippos); %
                %                 load('/mnt/data25/sjadhav/HPExpt/ProcessedDataDR/AllAn_PFCCA1_ripplepos_DR_vel4tet0_rippleposindex_vel4tet1filt');
                for oui = 1:length(indexALLlist);
                    indexALLShortList(oui,:) = str2num(sprintf('%d',indexALLlist(oui,:)));
                end
                for ru = 1:length(trajplotEpAvg); %all the cells in this set;; USING data averages across epochs****
                    subplot(3,length(trajplotEpAvg),ru); hold on; %make fig that fits all cells in set
                    for o = 1:length(trajplotEpAvg{ru}); %plot each trajectory
                        plot(trajplotEpAvg{ru}{o}(:,5),[clr{o} '.-'],'Linewidth',2);
                    end
                    if ru ==1; %add legend to pfc plot
                        legfix = legend('OL','IL','OR','IR');
                        legend(legfix, 'boxoff');
                        hText = findobj(legfix, 'type', 'text');
                        set(hText(4),'color',  [0 0 1],'FontSize', 7);
                        set(hText(3),'color',  [1 0 0]);
                        set(hText(2),'color',  [0 1 0]);
                        set(hText(1),'color',  [1 0 1]);
                        linesInPlot = findobj(legfix, 'type', 'line');
                        set(linesInPlot(1:8),'XData',[0 0]);
                        
                    elseif ru-1 <= poscelsinset; %length(pairdata(1,i).CA1posidxs(:,1)); %pos cell %elseif ru-1 <= length(pairdata(1,i).CA1posidxs(:,1)); %blue the pos CA1 indxs
                        typclr = [0 .7 .93];
                        set(gca,'Color',typclr);
                        set(gcf, 'InvertHardCopy', 'off');
                        
                    elseif ru-1 <= poscelsinset + negcelsinset; %length(pairdata(1,i).CA1sigidxs(:,1)); %neg cell %else %yellow the neg CA1 indxs
                        typclr = [1 1 .67];
                        set(gca,'Color',typclr);
                        set(gcf, 'InvertHardCopy', 'off');
                        
                    end
                    set(gca,'xtick',[])
                    set(gca,'xticklabel',[])
                    set(gca,'ytick',[])
                    set(gca,'yticklabel',[])
                    %firing map
                    subplot(3, length(mapplotEpAvg), ru+ length(trajplot)); hold on;
                    imagesc(mapplotEpAvg{ru});
                    axis([0, length(mapplotEpAvg{ru}(1,:)), 0, length(mapplotEpAvg{ru}(:,1))]);
                    if ru ==1;
                        tit = {sprintf('%s PFC %d %d %d', animals{pairdata(1,i).PFCidx(1,1)},pairdata(1,i).PFCidx(1,2), pairdata(1,i).PFCidx(1,3), pairdata(1,i).PFCidx(1,4))};  %label w index
                    else
                        tit = {sprintf('%s ca1 %d %d %d', animals{pairdata(1,i).CA1sigidxs(ru-1,1)}, pairdata(1,i).CA1sigidxs(ru-1,2), pairdata(1,i).CA1sigidxs(ru-1,3), pairdata(1,i).CA1sigidxs(ru-1,4)); sprintf('Sp(%0.5f)', CorrCoefstructSET(ru-1,1))}; %label w index
                    end
                    set(gca,'xtick',[])
                    set(gca,'xticklabel',[])
                    set(gca,'ytick',[])
                    set(gca,'yticklabel',[])
                    title(tit,  'Fontsize', titlesize,'FontName',fonttype)
                    %plot rip pos
                    currc1=[]; curr1allrips = [];
                    if ~isempty(find(trajindexshortlist(ru) == indexALLShortList));
                        currc1 = ripplepos{find(trajindexshortlist(ru) == indexALLShortList)};
                        currc1allrips = allripplemod(find(trajindexshortlist(ru) == indexALLShortList)).rippos;
                    end
                    if ~isempty(currc1allrips);
                        subplot(3, length(mapplotEpAvg), ru+ length(trajplot)*2); hold on;
                        plot(currc1allrips(:,1), currc1allrips(:,2),'.','LineWidth',0.5,'MarkerEdgeColor', [.8 .8 .8]);
                        if ~isempty(currc1);
                            plot(currc1(:,1), currc1(:,2),'x','LineWidth',2,'MarkerEdgeColor', 'k');
                        end
                        axis([min(currc1allrips(:,1))-5,max(currc1allrips(:,1))+5, min(currc1allrips(:,2))-5, max(currc1allrips(:,2))+5]);
                    end
                    set(gca,'xtick',[])
                    set(gca,'xticklabel',[])
                    set(gca,'ytick',[])
                    set(gca,'yticklabel',[])
                end
                figfile = [figdir,[savefigfilename '_maps'],'/','Maps_rippos',animals{pairdata(1,i).PFCidx(1,1)}, num2str(i)];
                if savefigs==1
                    print('-djpeg', figfile,'-r300');
                end
                if ~cyclemaps == 0
                    keyboard;
                end
                close
            end
            
         
            
            
%             __________________________________________________________________________________________________________________________________________________
            
        else
            pairsskippedPFC3Hz = pairsskippedPFC3Hz + length(pairdata(i).CA1sigidxs(:,1));
%             pfcskipped3Hz = pfcskipped3Hz +1;
            'pfc cell sets that didnt have peak rate >3Hz'
            PFCShortList(i) %print the skipped cell
        end %if pfc cell didn't have peakthresh rate
        
    end %for each pfc cell
end
% checking that the numbers make sense... the ca1skipped3Hz + pairsskippedPFC3Hz + #pos pairs used +#neg pairs used should = 84
%         pfcskipped3Hz
        ca1skipped3Hz
%         ca1posgood
%         ca1neggood
        pairsskippedPFC3Hz

%save the mod and unmod structures to run RipCorrVSSpatialCorr.m
%--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
if mod1unmod0 ~=2;
    clear modRs unmodRs
    if mod1unmod0 == 1;
        modRs = [sigbetasscatter, spatcorrvalsscatter];
        save(sprintf('%smodRs_signonsig',savedir), 'modRs');
    elseif mod1unmod0 == 0;
        unmodRs = [sigbetasscatter, spatcorrvalsscatter];
        save(sprintf('%sunmodRs_signonsig',savedir), 'unmodRs');
    end
end


%GATHER R VALUES FOR EACH ANIMAL(POS, NEG, SHUFF)
%--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%gather means of trajs  and run stats


    
    % posR = []; negR = []; allRShuff = [];
    %     for yu = 1:length(CorrCoefstructPOS); %collect mean corrcoef for each pos pair into list
    %         posR = [posR; (CorrCoefstructPOS{yu,1}.r(:))];
    %     end
    %
    %     for qu = 1:length(CorrCoefstructNEG); %collect mean corrcoef for each neg pair into list
    %         negR = [negR; (CorrCoefstructNEG{qu,1}.r(:))];
    %     end
    %
    %     for wu = 1:length(CorrCoefstructShuff); %collect mean corrcoef for each pair into list
    %         allRShuff = [allRShuff; (CorrCoefstructShuff{wu,1}.r(:))];
    %     end
    
    %     %for scatter betas vs spatial corr vals
    %     for juu = 1:length(CorrCoefstruct)
    %         sigbetas(juu,1) = CorrCoefstruct{juu}.sigbetas;
    %         spatcorrvals(juu,1) = CorrCoefstruct{juu}.r;
    %     end
    
    %POS
    hpa = 0; hpb = 0; hpc = 0; ndl = 0; posRHPa = []; posRHPb = []; posRHPc = []; posRNdl = [];
    for fu = 1:length(CorrCoefstructPOS); %POS
        if CorrCoefstructPOS{fu,1}.pfcindex(1,1) == 1;
            hpa = hpa+1;
            CorrCoefstructPOSHPa{hpa,1} = CorrCoefstructPOS{fu,1};
            %             posRHPa(hpa, 1) = mean(CorrCoefstructPOS{fu,1}.r(:));
            posRHPa = [posRHPa; (CorrCoefstructPOS{fu,1}.r(:))]; %vert stack
        elseif CorrCoefstructPOS{fu,1}.pfcindex(1,1) == 2;
            hpb = hpb+1;
            CorrCoefstructPOSHPb{hpb,1} = CorrCoefstructPOS{fu,1};
            %             posRHPb(hpb, 1) = mean(CorrCoefstructPOS{fu,1}.r(:));
            posRHPb = [posRHPb; (CorrCoefstructPOS{fu,1}.r(:))]; %vert stack
        elseif CorrCoefstructPOS{fu,1}.pfcindex(1,1) == 3;
            hpc = hpc+1;
            CorrCoefstructPOSHPc{hpc,1} = CorrCoefstructPOS{fu,1};
            %             posRHPc(hpc, 1) = mean(CorrCoefstructPOS{fu,1}.r(:));
            posRHPc = [posRHPc; (CorrCoefstructPOS{fu,1}.r(:))]; %vert stack
        elseif CorrCoefstructPOS{fu,1}.pfcindex(1,1) == 4;
            ndl = ndl+1;
            CorrCoefstructPOSNdl{ndl,1} = CorrCoefstructPOS{fu,1};
            %             posRNdl(ndl, 1) = mean(CorrCoefstructPOS{fu,1}.r(:));
            posRNdl = [posRNdl; (CorrCoefstructPOS{fu,1}.r(:))]; %vert stack
        end
    end
    
    %NEG
    hpa = 0; hpb = 0; hpc = 0; ndl = 0; fu = 0; negRHPa = []; negRHPb = []; negRHPc = []; negRNdl = [];
    for fu = 1:length(CorrCoefstructNEG); %NEG
        if CorrCoefstructNEG{fu,1}.pfcindex(1,1) == 1;
            hpa = hpa+1;
            CorrCoefstructNEGHPa{hpa,1} = CorrCoefstructNEG{fu,1};
            negRHPa = [negRHPa; (CorrCoefstructNEG{fu,1}.r(:))]; %vert stack
        elseif CorrCoefstructNEG{fu,1}.pfcindex(1,1) == 2;
            hpb = hpb+1;
            CorrCoefstructNEGHPb{hpb,1} = CorrCoefstructNEG{fu,1};
            negRHPb = [negRHPb; (CorrCoefstructNEG{fu,1}.r(:))]; %vert stack
        elseif CorrCoefstructNEG{fu,1}.pfcindex(1,1) == 3;
            hpc = hpc+1;
            CorrCoefstructNEGHPc{hpc,1} = CorrCoefstructNEG{fu,1};
            negRHPc = [negRHPc; (CorrCoefstructNEG{fu,1}.r(:))]; %vert stack
        elseif CorrCoefstructNEG{fu,1}.pfcindex(1,1) == 4;
            ndl = ndl+1;
            CorrCoefstructNEGNdl{ndl,1} = CorrCoefstructNEG{fu,1};
            negRNdl = [negRNdl; (CorrCoefstructNEG{fu,1}.r(:))]; %vert stack
        end
    end
    
    %SHUFF
    wu = 0; allRShuffHPa = []; allRShuffHPb = []; allRShuffHPc = []; allRShuffNdl = [];
    for wu = 1:length(CorrCoefstructShuffHPa); %collect mean corrcoef for each pair into list
        allRShuffHPa = [allRShuffHPa; (CorrCoefstructShuffHPa{wu,1}.r(:))]; %vert stack
    end
    wu = 0;
    for wu = 1:length(CorrCoefstructShuffHPb); %collect mean corrcoef for each pair into list
        allRShuffHPb = [allRShuffHPb; (CorrCoefstructShuffHPb{wu,1}.r(:))]; %vert stack
    end
    wu = 0;
    for wu = 1:length(CorrCoefstructShuffHPc); %collect mean corrcoef for each pair into list
        allRShuffHPc = [allRShuffHPc; (CorrCoefstructShuffHPc{wu,1}.r(:))]; %vert stack
    end
    wu = 0;
    for wu = 1:length(CorrCoefstructShuffNdl); %collect mean corrcoef for each pair into list
        allRShuffNdl = [allRShuffNdl; (CorrCoefstructShuffNdl{wu,1}.r(:))]; %vert stack
    end
    

    %GET RID OF NANs AND GET P VALUES USING RANK SUM AND KS TEST
    %--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    %     posR = posR(~isnan(posR));
    %     negR = negR(~isnan(negR));
    posRHPa = posRHPa(~isnan(posRHPa));
    negRHPa = negRHPa(~isnan(negRHPa));
    posRHPb = posRHPb(~isnan(posRHPb));
    negRHPb = negRHPb(~isnan(negRHPb));
    posRHPc = posRHPc(~isnan(posRHPc));
    negRHPc = negRHPc(~isnan(negRHPc));
    posRNdl = posRNdl(~isnan(posRNdl));
    negRNdl = negRNdl(~isnan(negRNdl));
    %     allRShuff = allRShuff(~isnan(allRShuff));
    allRShuffHPa = allRShuffHPa(~isnan(allRShuffHPa));
    allRShuffHPb = allRShuffHPb(~isnan(allRShuffHPb));
    allRShuffHPc = allRShuffHPc(~isnan(allRShuffHPc));
    allRShuffNdl = allRShuffNdl(~isnan(allRShuffNdl));
    
    %VERTICALLY CONCATENATE R VALUES FOR EACH AN TO CREATE GROUPED DATA
    posRAllAn = [posRHPa; posRHPb; posRHPc; posRNdl];
    negRAllAn = [negRHPa(negRHPa<.67); negRHPb; negRHPc; negRNdl];
    ShuffRAllAn = [allRShuffHPa; allRShuffHPb; allRShuffHPc; allRShuffNdl];
    
    
        %see note at top... horrible hack to get uncorr struct. 
    if plotCorrCoef == 1;
        if saveuncorr;
            clear ShuffRAllAn;
            ShuffRAllAn = posRAllAn;
            save(sprintf('%sShuffRAllAn',savedir), 'ShuffRAllAn');
            return
        else
            clear ShuffRAllAn;
            load(sprintf('%sShuffRAllAn',savedir));
        end
   
    
    
    
    
    %STATS
    [hposnegALL pposnegALL] = kstest2(posRAllAn,negRAllAn); %GROUPED DATA
    rsposnegALL = ranksum(posRAllAn,negRAllAn);
    [hposALL pposALL] = kstest2(posRAllAn,ShuffRAllAn);
    [hnegALL pnegALL] = kstest2(negRAllAn,ShuffRAllAn);
    rsposALL = ranksum(posRAllAn,ShuffRAllAn);
    rsnegALL = ranksum(negRAllAn,ShuffRAllAn);
    %     [hpos ppos] = kstest2(posR,allRShuff);
    %     [hneg pneg] = kstest2(negR,allRShuff);
    %     rspos = ranksum(posR,allRShuff);
    %     rsneg = ranksum(negR,allRShuff);
    
    %RS for individual animals.. won't work if less than 2 vals per category
    %     [hposHPa pposHPa] = kstest2(posRHPa,allRShuffHPa); %HPa
    %     [hposnegHPa pposnegHPa] = kstest2(posRHPa,negRHPa);
    %     [hnegHPa pnegHPa] = kstest2(negRHPa,allRShuffHPa);
    %     rsposHPa = ranksum(posRHPa,allRShuffHPa);
    %     rsnegHPa = ranksum(negRHPa,allRShuffHPa);
    %     rsposnegHPa = ranksum(posRHPa,negRHPa);
    %     [hposHPb pposHPb] = kstest2(posRHPb,allRShuffHPb); %HPb
    %
    %         [hposnegHPb pposnegHPb] = kstest2(posRHPb,negRHPb);
    %         [hnegHPb pnegHPb] = kstest2(negRHPb,allRShuffHPb);
    %         rsposHPb = ranksum(posRHPb,allRShuffHPb);
    %         rsposnegHPb = ranksum(posRHPb,negRHPb);
    %         rsnegHPb = ranksum(negRHPb,allRShuffHPb);
    %
    %     [hposHPc pposHPc] = kstest2(posRHPc,allRShuffHPc); %HPc
    %     [hnegHPc pnegHPc] = kstest2(negRHPc,allRShuffHPc);
    %     [hposnegHPc pposnegHPc] = kstest2(posRHPc,negRHPc);
    %     rsposHPc = ranksum(posRHPc,allRShuffHPc);
    %     rsposnegHPc = ranksum(posRHPc,negRHPc);
    %     rsnegHPc = ranksum(negRHPc,allRShuffHPc);
    %     [hposNdl pposNdl] = kstest2(posRNdl,allRShuffNdl); %Nadal
    %     [hposnegNdl pposnegNdl] = kstest2(posRNdl,negRNdl);
    %     [hnegNdl pnegNdl] = kstest2(negRNdl,allRShuffNdl);
    %     rsposNdl = ranksum(posRNdl,allRShuffNdl);
    %     rsposnegNdl = ranksum(posRNdl,negRNdl);
    %     rsnegNdl = ranksum(negRNdl,allRShuffNdl);
    
    
    %update 3/28/14.. Using a Krusak Wallis test (a non-parametric version of a 1way anova for more than 2 groups).. then correcting for multiple comparisons.
    
    % posstr = zeros(1, length(posRAllAn));
    % negstr = ones(1,length(negRAllAn)) .*2
    
    [pKW table statsKW] = kruskalwallis([posRAllAn' ShuffRAllAn' negRAllAn'], [ones(1, length(posRAllAn)) ones(1,length(ShuffRAllAn)).*2 ones(1,length(negRAllAn)).*3])
    [c, m, h, gnames] = multcompare(statsKW, 'estimate', 'kruskalwallis', 'ctype', 'hsd', 'display', 'on', 'alpha', 0.06) % change the alpha values to determine significance range.
    %update 4/11/14.. now using only simultaneously recorded neurons.. current p vals. neg-shuff >.05 (not sig); pos-neg < .0001 pos-shuff <.0001
    %update 4/11/14... just use the ranksum values and 0.05/3 as the mult compare sig thresh..
    %old results w sampling: neg-shuff <.05; pos-neg < .0001 pos-shuff <.00001
    
    
    %update 4/3/14... what % of times is the mean of a subset with matched # samples from the shuffled set met or exceeded the magnititude of the mean of the test distribution
    if runsampling == 1;
        for yy = 1:1000000;
            set4pos(yy) = nanmean(datasample(ShuffRAllAn, length(posRAllAn), 'Replace', false));
            set4neg(yy) = nanmean(datasample(ShuffRAllAn, length(negRAllAn), 'Replace', false));
            set4posneg(yy) = nanmean(datasample(posRAllAn, length(negRAllAn), 'Replace', false));
        end
        pval4pos = length(set4pos(set4pos >= nanmean(posRAllAn)))/length(set4pos);
        pval4neg = length(set4neg(set4neg <= nanmean(negRAllAn)))/length(set4neg);
        pval4posneg = length(set4posneg(set4posneg <= nanmean(negRAllAn)))/length(set4neg);
    end
    
    %--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    %PLOT BAR FIGS
    %______________________________________________________________________________________
    %plot all anims grouped
    
    sigthresh = 0.05/3
    'rsposnegALL' 
    rsposnegALL
    'rsnegALL'
    rsnegALL
    'rsposALL'
     rsposALL
    
    figure
    hold on
    %     subplot(2,1,1)
    
    bar([.5 1 1.5], [mean(posRAllAn) mean(ShuffRAllAn) mean(negRAllAn)], .9, 'EdgeColor', 'none')
    %     bar(1, mean(posRAllAn) , 'FaceColor', [0 .6 .83], 'EdgeColor', 'none')
    % bar(1, mean(posRAllAn) , 'FaceColor', [.0 .50 .75], 'EdgeColor', 'none')
    
    %     bar(2, mean(ShuffRAllAn), 'FaceColor', [.5 .5 .5], 'EdgeColor', 'none')
    %     bar(3, mean(negRAllAn), 'FaceColor', [1 1 .67], 'EdgeColor', 'none')
    %     bar(3, mean(negRAllAn), 'FaceColor', [.93 .80 .38], 'EdgeColor', 'none')
    %     errorbar2([1 2 3], [mean(posRAllAn) mean(ShuffRAllAn) mean(negRAllAn) ],  [std(posRAllAn) std(ShuffRAllAn) std(negRAllAn)] , 0.3, 'k')
    errorbar2([.5 1 1.5], [mean(posRAllAn) mean(ShuffRAllAn) mean(negRAllAn) ],  [stderr(posRAllAn) stderr(ShuffRAllAn) stderr(negRAllAn)] , 0.3, 'k')
    set(gcf,'PaperPositionMode','auto');
    set(gcf, 'Position', [100 100 500 500])
    
    %now adding combined cells data
    %     CorrCoefcombPOS = CorrCoefcombPOS(~isnan(CorrCoefcombPOS));
    %     CorrCoefcombNEG = CorrCoefcombNEG(~isnan(CorrCoefcombNEG));
    %     CorrCoefcombPOSNEG = CorrCoefcombPOSNEG(~isnan(CorrCoefcombPOSNEG));
    %     bar(4, mean(CorrCoefcombPOS), 'FaceColor', [0 .6 .83], 'EdgeColor', 'none')
    %     bar(5, mean(CorrCoefcombNEG), 'FaceColor', [.9 .9 .57], 'EdgeColor', 'none')
    %     bar(6, mean(CorrCoefcombPOSNEG), 'FaceColor', [0 1 0], 'EdgeColor', 'none')
    %     errorbar2([4 5 6], [mean(CorrCoefcombPOS) mean(CorrCoefcombNEG) mean(CorrCoefcombPOSNEG) ],  [std(CorrCoefcombPOS) std(CorrCoefcombNEG) std(CorrCoefcombPOSNEG)] , 0.3, 'k')
    %     rsposC = ranksum(CorrCoefcombPOS,posRAllAn);
    %     rsnegC = ranksum(CorrCoefcombNEG,negRAllAn);
    %     rsposnegC = ranksum(CorrCoefcombPOSNEG,posRAllAn);
    
    %data points
    %     plot([1], [posRAllAn ] , 'ko', 'MarkerSize',5,'LineWidth',.5)
    %     plot([3], [negRAllAn ] , 'ko', 'MarkerSize',5,'LineWidth',.5)
    %     plot([4], [CorrCoefcombPOS ] , 'ko', 'MarkerSize',5,'LineWidth',.5)
    %     plot([5], [CorrCoefcombNEG ] , 'ko', 'MarkerSize',5,'LineWidth',.5)
    %     plot([6], [CorrCoefcombPOSNEG ] , 'ko', 'MarkerSize',5,'LineWidth',.5)
    
    %     sigstar({[1,2]},[rsposALL,rsnegALL],1)
    
    %     sigstar({[1,2],[2,3], [1,3], [1,4], [3,5],[1,6]},[rsposALL,rsnegALL,rsposnegALL, rsposC, rsnegC, rsposnegC ])
    % sigstar({[1,2],[2,3], [1,3]},[rsposALL,rsnegALL,rsposnegALL,])
    
    %     errorbar2([1 2 3], [mean(posRAllAn) mean(ShuffRAllAn) mean(negRAllAn) ],  [stderr(posRAllAn) stderr(ShuffRAllAn) stderr(negRAllAn)] , 0.3, 'k')
    %         xlim([0.3 2.7])
    set(gca, 'Fontsize', titlesize,'FontName',fonttype)
    %     set(gca, 'xtick', [1 2 3 4 5 6], 'xticklabel', {'POS', 'AllShuff', 'NEG', 'POSCOM', 'NEGCOM', 'POSNEGCOM'})
    set(gca, 'xtick', [.5 1 1.5], 'xticklabel', {'Positive', 'Non-significant', 'Negative'})
    ylabel('Spatial Correlation')
    %     title({['Rip Corr PFC-CA1 All Animals']; sprintf('ranksum pos-shuffP = %d  neg-shuffP = %d pos-negP = %d',rsposALL, rsnegALL, rsposnegALL); sprintf('KStest pos-shuffP = %d neg-shuffP = %d pos-negP = %d ',pposALL,pnegALL, pposnegALL)}, 'FontSize', 10)
    %     title({['Rip Corr PFC-CA1 All Animals']; sprintf('RS pos:shuff(%0.5f) neg:shuff(%0.5f) pos:neg(%0.5f)',rsposALL, rsnegALL, rsposnegALL); sprintf('pos:posC(%0.5f) negC:neg(%0.5f) posnegC:pos(%0.5f) ',rsposC,rsnegC, rsposnegC)}, 'FontSize', 15)
    %         title({['PFC-CA1 Pairs']}, 'FontSize', titlesize)
    set(gca,'XLim',[.2 1.8]);
    %     set(gca,'XLim',[.25 3.75]);
    set(gca,'YLim',[-.2 .2]);
    
    figfile = [figdir,savefigfilename,'/',sprintf('%s_AllAn', savefigfilename)];
    if savefigs==1
        print('-dpdf', figfile, '-r300');
    end
    if ~cyclemaps == 0
        keyboard;
    end
    close
    
    %violin plot
    figure; hold on; distributionPlot([{posRAllAn} {ShuffRAllAn} {negRAllAn}], 'color', [{[.0 .50 .75]} {[.3 .3 .3]} {[.93 .80 .38]}], 'showMM', 4)
    sigstar({[1,2],[2,3], [1,3]},[rsposALL,rsnegALL,rsposnegALL,])
    set(gca, 'Fontsize', titlesize,'FontName',fonttype)
    set(gca, 'xtick', [1 2 3], 'xticklabel', {'rip-POS', 'Shuffled', 'rip-NEG'})
    ylabel('Corr Coef')
    set(gca,'XLim',[.25 3.75]);
    set(gca,'YLim',[-.6 .75]);
    title({['PFC-CA1 Pairs']}, 'FontSize', titlesize)
    
    figfile = [figdir,savefigfilename,'/',sprintf('%s_AllAnViolin', savefigfilename)];
    if savefigs==1
        print('-dpdf', figfile, '-r300');
    end
    if ~cyclemaps == 0
        keyboard;
    end
    close
    
    %         %violin plot with combined
    %     figure; hold on; distributionPlot([{posRAllAn} {ShuffRAllAn} {negRAllAn} {CorrCoefcombPOS} {CorrCoefcombNEG} {CorrCoefcombPOSNEG}], 'color', [{[.0 .50 .75]} {[.3 .3 .3]} {[.93 .80 .38]} {[.0 .50 .75]} {[.3 .3 .3]} {[0 1 0]}], 'showMM', 4)
    %     set(gca, 'Fontsize', titlesize,'FontName',fonttype)
    %     ylabel('Corr Coef')
    %             title({['PFC-CA1 Pairs']}, 'FontSize', titlesize)
    %     figfile = [figdir,savefigfilename,'/',sprintf('%s_AllAnViolinComb', savefigfilename)];
    %     if savefigs==1
    %         print('-dpdf', figfile, '-r300');
    %     end
    %     if ~cyclemaps == 0
    %         keyboard;
    %     end
    %     close
    
    
    % % scatter old
    %     figure
    %     hold on
    %   scatter(sigbetasscatter, spatcorrvalsscatter)%,'MarkerEdgeColor','k','MarkerFaceColor',[.5 .5 .5],'LineWidth',1)
    %   %regression
    %       [b,bint,r,rint,stats] = regress(sigbetasscatter, spatcorrvalsscatter);
    %       [r2 p2]= corrcoef(sigbetasscatter(spatcorrvalsscatter>0), spatcorrvalsscatter(spatcorrvalsscatter>0));
    %     xpts = min(sigbetasscatter):0.01:max(sigbetasscatter);
    %     bfit = b(1)*xpts;
    %     plot(xpts,bfit,'k-','LineWidth',1);  % Theta vs Rip
    %       set(gca,'XLim',[min(sigbetasscatter) max(sigbetasscatter)]);
    %   set(gca,'YLim',[min(spatcorrvalsscatter) max(spatcorrvalsscatter)]);
    %   ylabel('SpatCorrCoef')
    %  xlabel('RipCorrR')
    %   title({['RipCorrR x SpatCorrCoefs']; sprintf('R^2(%0.5f) p(%0.5f)', r2(2,1), p2(2,1))});
    %
    %
    %     figfile = [figdir,savefigfilename,'/',sprintf('%s_AllAnScatter', savefigfilename)];
    %     if savefigs==1
    %         print('-dpdf', figfile, '-r300');
    %     end
    %     if ~cyclemaps == 0
    %         keyboard;
    %     end
    %     close
    
    return
    %each animal... NOT USING THIS ANYMORE
%     figure(2)
%     subplot(2,2,1);
%     hold on
%     bar(1, mean(posRHPa), 'FaceColor', [0 .7 .93], 'EdgeColor', 'none')
%     bar(2, mean(allRShuffHPa),'k', 'EdgeColor', 'none')
%     bar(3, mean(negRHPa),'FaceColor', [1 1 .67], 'EdgeColor', 'none')
%     %     errorbar2([1 2 3], [mean(posRHPa) mean(allRShuffHPa) mean(negRHPa) ],  [std(posRHPa) std(allRShuffHPa) std(negRHPa)] , 0.3, 'k')
%     %     errorbar2([1 2 3], [mean(posRHPa) mean(allRShuffHPa) mean(negRHPa) ],  [stderr(posRHPa) stderr(allRShuffHPa) stderr(negRHPa)] , 0.3, 'k')
%     %         xlim([0.3 2.7])
%     plot([1], [posRHPa ] , 'ko', 'MarkerSize',5,'LineWidth',.5)
%     plot([3], [negRHPa ] , 'ko', 'MarkerSize',5,'LineWidth',.5)
%     
%     set(gca, 'fontsize', 10)
%     set(gca, 'xtick', [1 2 3], 'xticklabel', {'POS', 'HPaShuff', 'NEG'})
%     ylabel('CorrCoef')
%     %     title({['Rip Corr PFC-CA1 HPA']; sprintf('ranksum posP = %0.5f  negP = %0.5f posnegP = %0.5f',rsposHPa, rsnegHPa, rsposnegHPa); sprintf('KStest posP = %0.5f negP = %0.5f negP = %0.5f',pposHPa,pnegHPa, pposnegHPa )}, 'FontSize', 10)
%     %     title({['Rip Corr PFC-CA1 HPA']; sprintf('RS posP(%0.5f) negP(%0.5f) posnegP(%0.5f)',rsposHPa, rsnegHPa, rsposnegHPa)},'Fontsize', titlesize,'FontName',fonttype)
%     title({['Rip Corr PFC-CA1 HPA']},'Fontsize', titlesize,'FontName',fonttype)
%     
%     subplot(2,2,2);
%     hold on
%     bar(1, mean(posRHPb), 'FaceColor', [0 .7 .93], 'EdgeColor', 'none')
%     bar(2, mean(allRShuffHPb),'k', 'EdgeColor', 'none')
%     bar(3, mean(negRHPb),'FaceColor', [1 1 .67], 'EdgeColor', 'none')
%     %     errorbar2([1 2 3], [mean(posRHPb) mean(allRShuffHPb) mean(negRHPb) ],  [std(posRHPb) std(allRShuffHPb) std(negRHPb)] , 0.3, 'k')
%     %     errorbar2([1 2 3], [mean(posRHPb) mean(allRShuffHPb) mean(negRHPb) ],  [stderr(posRHPb) stderr(allRShuffHPb) stderr(negRHPb)] , 0.3, 'k')
%     %         xlim([0.3 2.7])
%     plot([1], [posRHPb] , 'ko', 'MarkerSize',5,'LineWidth',.5)
%     plot([3], [negRHPb] , 'ko', 'MarkerSize',5,'LineWidth',.5)
%     set(gca, 'fontsize', 10)
%     set(gca, 'xtick', [1 2 3], 'xticklabel', {'POS', 'HPbShuff', 'NEG'})
%     ylabel('CorrCoef')
%     %     title({['Rip Corr PFC-CA1 HPb']; sprintf('ranksum posP = %0.5f  negP = %0.5f posnegP = %0.5f',rsposHPb, rsnegHPb, rsposnegHPb); sprintf('KStest posP = %0.5f negP = %0.5f posnegP = %0.5f ',pposHPb,pnegHPb, pposnegHPb)},'FontSize', 10)
%     %     title({['Rip Corr PFC-CA1 HPb']; sprintf('RS posP(%0.5f) negP(%0.5f) posnegP(%0.5f)',rsposHPb, rsnegHPb, rsposnegHPb)},'Fontsize', titlesize,'FontName',fonttype)
%     title({['Rip Corr PFC-CA1 HPb']},'Fontsize', titlesize,'FontName',fonttype)
%     
%     subplot(2,2,3);
%     hold on
%     bar(1, mean(posRHPc), 'FaceColor', [0 .7 .93], 'EdgeColor', 'none')
%     bar(2, mean(allRShuffHPc),'k', 'EdgeColor', 'none')
%     bar(3, mean(negRHPc),'FaceColor', [1 1 .67], 'EdgeColor', 'none')
%     %     errorbar2([1 2 3], [mean(posRHPc) mean(allRShuffHPc) mean(negRHPc) ],  [std(posRHPc) std(allRShuffHPc) std(negRHPc)] , 0.3, 'k')
%     %     errorbar2([1 2 3], [mean(posRHPc) mean(allRShuffHPc) mean(negRHPc) ],  [stderr(posRHPc) stderr(allRShuffHPc) stderr(negRHPc)] , 0.3, 'k')
%     %         xlim([0.3 2.7])
%     plot([1], [posRHPc] , 'ko', 'MarkerSize',5,'LineWidth',.5)
%     plot([3], [negRHPc] , 'ko', 'MarkerSize',5,'LineWidth',.5)
%     set(gca, 'fontsize', 10)
%     set(gca, 'xtick', [1 2 3], 'xticklabel', {'POS', 'HPcShuff', 'NEG'})
%     ylabel('CorrCoef')
%     %     title({['Rip Corr PFC-CA1 HPc']; sprintf('ranksum posP = %0.5f  negP = %0.5f posnegP = %0.5f',rsposHPc, rsnegHPc,rsposnegHPc); sprintf('KStest posP = %0.5f negP = %0.5f posnegP = %0.5f',pposHPc,pnegHPc, pposnegHPc)}, 'FontSize', 10)
%     %     title({['Rip Corr PFC-CA1 HPc']; sprintf('RS posP(%0.5f) negP(%0.5f) posnegP(%0.5f)',rsposHPc, rsnegHPc,rsposnegHPc)}, 'Fontsize', titlesize,'FontName',fonttype)
%     title({['Rip Corr PFC-CA1 HPc']}, 'Fontsize', titlesize,'FontName',fonttype)
%     
%     subplot(2,2,4);
%     hold on
%     bar(1, mean(posRNdl), 'FaceColor', [0 .7 .93], 'EdgeColor', 'none')
%     bar(2, mean(allRShuffNdl),'k', 'EdgeColor', 'none')
%     bar(3, mean(negRNdl),'FaceColor', [1 1 .67], 'EdgeColor', 'none')
%     %     errorbar2([1 2 3], [mean(posRNdl) mean(allRShuffNdl) mean(negRNdl) ],  [std(posRNdl) std(allRShuffNdl) std(negRNdl)] , 0.3, 'k') %pos and neg only have 1
%     %     errorbar2([1 2 3], [mean(posRNdl) mean(allRShuffNdl) mean(negRNdl) ],  [stderr(posRNdl) stderr(allRShuffNdl) stderr(negRNdl)] , 0.3, 'k')
%     %         xlim([0.3 2.7])
%     plot([1], [posRNdl] , 'ko', 'MarkerSize',5,'LineWidth',.5)
%     plot([3], [negRNdl] , 'ko', 'MarkerSize',5,'LineWidth',.5)
%     
%     set(gca, 'fontsize', 10)
%     set(gca, 'xtick', [1 2 3], 'xticklabel', {'POS', 'NdlShuff', 'NEG'})
%     ylabel('CorrCoef')
%     %     title({['Rip Corr PFC-CA1 Ndl']; sprintf('ranksum posP = %0.5f  negP = %0.5f posnegP = %0.5f',rsposNdl, rsnegNdl, rsposnegNdl); sprintf('KStest posP = %0.5f negP = %0.5f posnegP = %0.5f',pposNdl,pnegNdl, pposnegNdl)}, 'FontSize', 10)
%     %     title({['Rip Corr PFC-CA1 Ndl']; sprintf('RS posP(%0.5f)  negP(%0.5f) posnegP(%0.5f)',rsposNdl, rsnegNdl, rsposnegNdl)}, 'Fontsize', titlesize,'FontName',fonttype)
%     title({['Rip Corr PFC-CA1 Ndl']}, 'Fontsize', titlesize,'FontName',fonttype)
%     
%     figfile = [figdir,savefigfilename,'/',sprintf('%s_EachAn', savefigfilename)];
%     if savefigs==1
%         print('-djpeg', figfile, '-r300');
%     end
%     if ~cyclemaps == 0
%         keyboard;
%     end
%     close
end




