
warning('off','all');
% clear; close all;

runfiltflds = 0;

runanyscripts = 1;% won't run if runfiltflds
runscriptmaps = runanyscripts; 
runPECv2 = runanyscripts;
runscriptsparsity = runanyscripts;
runripposscript = 1;
runscriptISI = 0; 
runPEC = 0; %don't use this

plotanything = 1; %must be on for any plotting
savefigs = 0;
pausefigs = 1; %pause at each fig, keyboard press to continue

plotCorrCoef = 0; %population pairwise cc
plotPFCripmodMapsTrajsISImeanrate = 1; %plot the maps and all the other stuff.. update, right now it just plots maps and trajs and overlap
plotPECall = 0; %population PE stats, right now this includes scatter plots that loren asked for

plotPFCCA1Trajs = 0; %not really using this anymore
plotPFCCA1MapsTrajs =0; %not really using this anymore
plotPFCripmodMapsTrajs = 0; %not really using this anymore
plotPFCripmodMeanRates = 0; %not really using this anymore
plotPFCripmodISI = 0; %%not really using this anymore



savedir = '/data19/sjadhav/HPExpt/ProcessedDataDR/';
% savefilename = sprintf('Flds_%s',date); %has to match saved data name. specify GLM or Corr pairs and filters used: velocity filter <=.....nrip >= (#tetrodes ripples detected); std > of ripple detection power
savefilename = 'Flds_jan26';
% savefilename = 'Flds_feb9';
% savefigfilename = sprintf('SpatRipcorr_%s', date); %fig name
savefigfilename = sprintf('RipExcInhMaps_%s', date); %fig name
savefile = [savedir savefilename]; %area = 'PFC'; %clrunmod = 'r'; clrmod = 'b'; % PFC
figdir = '/data19/sjadhav/HPExpt/Figures_DR/';  %figdir = '/data19/sjadhav/HPExpt/Figures_DR/';
% loadrippos = '/data19/sjadhav/HPExpt/ProcessedDataDR/AllAn_PFCCA1_ripplepos_DR_vel5tet1'; %load the rip positions from DFSsj_HPexpt_getripalignspikingGRAllPosition_DR.m
loadpairindices = '/data19/sjadhav/HPExpt/HP_ProcessedData/corrindsForSpatial'; %load the correct indices!! specificy GLM or Corr Pairs in savefilename above..
minV = 3;

% If runscript, run Datafilter and save data
if runfiltflds == 1
    tic;
    
    %Animal selection
    %-----------------------------------------------------
    animals = {'HPa' 'HPb' 'HPc' 'Ndl' 'Rtl' 'Brg'};
    %             animals = {'HPa' 'HPb' 'HPc'};
    %         animals = {'HPc'};
    %     animals = {'HPa'};
    %         animals = {'nadal'};
    
    %Filter creation
    %-----------------------------------------------------
    % Epoch filter
    % -------------
    
    dayfilter = ''; %leave blank to take all days from HP animals and Ndl
    runepochfilter = 'isequal($type, ''run'') && ~isequal($environment, ''lin'')';
    
    % %Cell filter
    % %-----------
    %     placecellfilter = '(strcmp($area, ''PFC'') || (strcmp($area, ''CA1'') || && ($numspikes > 100))';  % not mod/unmod
    %     placecellfilter = '(strcmp($area, ''PFC'') && ($numspikes > 100) || strcmp($area, ''CA1'') && ($numspikes > 100))';
    %     placecellfilter = '( strcmp($tag, ''CA1Pyr'') || strcmp($tag, ''iCA1Pyr'') || strcmp($tag, ''PFC'')) && ($numspikes > 100)';
    placecellfilter = '( strcmp($area, ''CA1'') || strcmp($area, ''iCA1'') || strcmp($area, ''PFC'')) && ($numspikes > 100)';
    
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
%     psf = setfilterfunction(spatf, 'DFAsj_filtercalclinfields_tf',{'spikes', 'linpos'}, 'binsize', 2);
spks = setfilterfunction(spatf, 'DFAdr_spikes',{'spikes', 'linpos', 'pos'}, 'binsize', 1, 'std', 2); %spikes for ISI
fields = setfilterfunction(spatf, 'DFA_loadfields_DR', {'linfields', 'mapfields', 'cellinfo'}); %map and lin fields
spatinf = setfilterfunction(spatf, 'getspatialinfo_DR2', {'linpos', 'spikes'}, 6,  minV, 'peakrate', 3,'appendindex', 1, 'incltraj', [1:4]); %I changed it to output all trajs' summed per cell, not combined across eps yet
% modf = setfilterfunction(spatf,'DFAsj_getripalignspiking_position',{'spikes', 'ripples', 'tetinfo', 'pos'}); %ripple position

% Run analysis-----------------------
%     flds = runfilter(psf);  % Place Field Stability.. trajectories
spks = runfilter(spks);  %spatf = setfilterfunction(spatf, 'getspatialinfo_DR', {'linpos', 'spikes'},6,  minV, 'peakrate', 3,'appendindex', 1, 'incltraj', traj); %I changed it to output all trajs' summed per cell, not combined across eps yet spike data in order to get ISI
flds = runfilter(fields);  % Place Field Map
spatinf = runfilter(spatinf); %spatial information
% ripP = runfilter(modf); %ripple position

%     end
disp('Finished running filter script');
%--------------------- Finished Filter Function Run -------------------
%clear some variables so that when the mat gets loaded in, it doesn't override the intended variable vals.
clear  runfiltflds runscriptsparsity runanyscripts runripposscript runscript runPECv2 runPEC runscriptmaps plotPFCripmodISI plotPFCripmodMapsTrajsISImeanrate cyclefigs savefigs  plotPFCCA1MapsTrajs plotPFCCA1Trajs plotCorrCoef loadrippos  loadpairindices savefilename savedir  savefigfilename  savedata  figdir savefigs pausefigs plotanything plotPFCripmodMapsTrajs plotPFCripmodMeanRates
save(savefile);
toc;
return
% else
%     load(savefile);
end % end runfilter

%
% if ~exist('savedata')
%     return
% end

% -------------------------  Filter Format Done -------------------------

%prepping data structs
%________________________________________________________________________
tic
% if runscript == 1;
load(savefile);
if runanyscripts ==1;
    if runscriptmaps == 1;
        load(loadpairindices); %GR/SJ provided ep combined rip corr struct
        matpairs = cell2mat(arrayfun(@(x) [x.index x.allp_epcomb x.allr_epcomb], corrindsForSpatial, 'UniformOutput', false)'); %collect all pfc, ca1 pair indices, p vals, and r vals. %N[cells] x 8[an day ca1tet ca1cell pfctet pfccell pval rval]
        
        %collect ALL of the flds data into workable format
        cnt = 0;
        for anims = 1:length(animals);
            for cinds = 1:length(flds(anims).output{1});
                cnt = cnt+1;
                fldsdata(cnt,1) = flds(anims).output{1}(cinds);
                fldsdata(cnt,1).index = [anims fldsdata(cnt,1).index]; %store the animal # in with the other indices. [an day ep tet cell]
                matfldsinds(cnt,:) = [anims flds(anims).output{1}(cinds).index([1 3 4])]; %anim day tet cell
            end
        end
        
        %mean/combine the flds data across epochs and find the peak rate for each cell..
        [uniqmatallcells] = unique(matfldsinds, 'rows', 'stable');
        for i = 1:length(uniqmatallcells);
            fldsmatch{i} = find(~any(bsxfun(@minus,matfldsinds,uniqmatallcells(i,:)),2)); % return the row indices of the flds matrix with rows that match the index vector from each row of the unique cells ind matrix. ' minus' because matching vectors will return a zero vec.
            %i could have also done the line above using ismember
            epcombcell = cell(1,4);
            skipcell = 0;
            for j = 1:length(fldsmatch{i})
                try
                    epcombcell = [epcombcell; cellfun(@(x) x(:,4:5)', fldsdata(fldsmatch{i}(j)).trajdata, 'UniformOutput', false)]; %this is grabbing the [unsmoothed and smoothed] trajs
                    epcombcellarea = fldsdata(fldsmatch{i}(j)).area; %this is grabbing the area
                catch
                    epcombcell{1} = NaN;  %if there isn't data for a particular traj or map data for this cell
                    skipcell = 1;
                end
            end
            if ~skipcell; %if there are traj and map data, truncate it and store it
                epcombcell = epcombcell(2:end,:); %clip the empty first row that was used to initialize
                minlentraj = mat2cell(repmat(min(cell2mat(cellfun(@(x) length(x), epcombcell,'UniformOutput',false)),[],1),length(epcombcell(:,1)),1), length(epcombcell(:,1)),4); %find the min length of each traj across epochs and create a row of these vals for each epoch.. turn into cell
                trunceps = cellfun(@(x,y) x(2,1:y), epcombcell, minlentraj, 'UniformOutput', false); %truncated epochs according to the shortest trajectory across epochs.
                truncepsNOTSMOOTH = cellfun(@(x,y) x(1,1:y), epcombcell, minlentraj, 'UniformOutput', false); %truncated epochs according to the shortest trajectory across epochs. UNSMOOTHED FOR PE
                truncepcomb = cell(1,4);
                truncepcombNOTSMOOTH = cell(1,4);
                for k = 1:4;
                    truncepcomb{1,k} = vertcat(trunceps{:,k});
                    truncepcombNOTSMOOTH{1,k} = vertcat(truncepsNOTSMOOTH{:,k});
                end
                fldsepcomb{i,1} =  cellfun(@(x) nanmean(x, 1), truncepcomb, 'UniformOutput', false); %ep combined fldsdata
                fldsepcomb{i,2} = uniqmatallcells(i, :); %horz concatenate pair indices for debugging
                fldspeak(i,1) = max(cell2mat(cellfun(@(x) max(max(max(x, 1))), truncepcomb, 'UniformOutput', false))); %finds the peak rate, BEFORE AVERAGING ACROSS EPOCHS, across trajs and across eps for each unique cell.. use for threshing check
                fldsepcombNOTSMOOTH{i,1} =  cellfun(@(x) nanmean(x, 1), truncepcombNOTSMOOTH, 'UniformOutput', false); %ep combined fldsdata not smoothed
                fldsepcombNOTSMOOTH{i,2} = uniqmatallcells(i, :); %horz concatenate pair indices for debugging
                if lower(epcombcellarea(1)) == 'p'; %if pfc
                    fldsepcombNOTSMOOTH{i,3} = 2;
                else
                    fldsepcombNOTSMOOTH{i,3} = 1; %else ca1
                end
            else
                fldsepcomb{i,1}{1}(1) = NaN;
                fldsepcomb{i,2} = uniqmatallcells(i, :);
                fldspeak(i,1) = NaN;
            end

            
            

            
            
            %%
            
            %         %% calculate the path equivalent coefficient (PEC) for each cell.. based on measure in Frank et al 2000.. this currently won't work for only 3 trajs.. need to fix..
            %         %also need to include some filter on including only those trajs with a 'field'.. meh i don't think i need to do this as im taking the max anyway
            %         %i DO need to remove the center arm part for shared paths.. done i think
            %         if runPEC && ~skipcell;
            %             PEtrajs = fldsepcombNOTSMOOTH{i,1}; % don't overwrite the flds data
            %             minlentraj = mat2cell(repmat(min(cell2mat(cellfun(@(x) length(x), PEtrajs,'UniformOutput',false))),1,length(PEtrajs)),1,length(PEtrajs)); %shortest trajectory
            %             trunctrjs = cellfun(@(x,y) x(1:y), PEtrajs, minlentraj, 'UniformOutput', false); %truncated trajs according to the shortest trajectory i'm clipping from the end because i haven't flipped the inbound trajs yet, so these should all be in terms of distance from center well.
            %             if length(fldsepcomb{i,1}) == 4; %i'm not sure which traj is which if there's less than 4 right now.. but i think this data is somewhere around here.. to do
            %                 PEtrajs{2} = fliplr(PEtrajs{2}); PEtrajs{4} = fliplr(PEtrajs{4});
            %             end
            %             permtraj = mat2cell(nchoosek(1:length(PEtrajs),2),length(nchoosek(1:length(PEtrajs),2)),1); % permutations of trajs
            %             gaussfilter = fspecial('gaussian',[1 6],1); %create 6 point gaussian filter with std 1
            %             trunctrjsSMOOTH = cellfun(@(x) conv(x,gaussfilter,'same'), trunctrjs, 'UniformOutput', false); % smooth
            %             %this needs to be a for loop because i need to throw out the center arm for 1-3 and 2-4 pairs
            %             for h = 1:length(permtraj);
            %                 if permtraj{h} == [1 3]; %if it's a combo with common center arm direction of movement, i.e. outleft and outright is 1,3  %clip the first half of 1,3
            %                     outleft = trunctrjsSMOOTH{permtraj{h}(1)};
            %                     outright = trunctrjsSMOOTH{permtraj{h}(2)};
            %                     outleft = outleft(floor(length(outleft)/2):end);
            %                     outright = outright(floor(length(outright)/2):end);
            %                     PEcorr{h} = corrcoef([outleft' outright'], 'rows', 'pairwise');
            %                 elseif permtraj{h} == [2 4]; %if it's a combo with common center arm direction of movement, i.e. infromleft and infromright is 2,4
            %                     %clip the last half of 2,4...
            %                     inleft = trunctrjsSMOOTH{permtraj{h}(1)};
            %                     inright = trunctrjsSMOOTH{permtraj{h}(2)};
            %                     inleft = inleft(1:floor(length(inleft)/2));
            %                     inright = inright(1:floor(length(inright)/2));
            %                     PEcorr{h} = corrcoef([inleft' inright'], 'rows', 'pairwise');
            %                 else %just use the whole traj
            %                     PEcorr{h} = corrcoef([trunctrjsSMOOTH{permtraj{h}(1)}' trunctrjsSMOOTH{permtraj{h}(2)}'], 'rows', 'pairwise');
            %                 end
            %             end
            %             %             PEcorr = cellfun(@(x) corrcoef([trunctrjsSMOOTH{x(1)}' trunctrjsSMOOTH{x(2)}'], 'rows', 'pairwise'), permtraj,'UniformOutput',false);
            %             PEcorrmax = max(cell2mat(cellfun(@(x) x(1,2), PEcorr,'UniformOutput',false)));
            %             PEtrajsShuff = [trunctrjs cellfun(@(x) [x(floor(length(x)/2)+1:end) fliplr(x(1:floor(length(x)/2)))], trunctrjs,'UniformOutput',false)]; %A....BC....D becomes C....DB....A.. Also i'm putting the real trajs as the first 4 cells
            %             permtrajsSHUFF = combvec(1:length(trunctrjs), length(trunctrjs)+1:length(trunctrjs)*2); % create a list of permutations of two vectors.. 1:4 vs 5:8 typically
            %             permtrajsSHUFFs = mat2cell(permtrajsSHUFF,1,length(permtrajsSHUFF)); %cellify
            %             trunctrjsSHUFFSMOOTH = cellfun(@(x) conv(x,gaussfilter,'same'), PEtrajsShuff, 'UniformOutput', false); % smooth
            %             PEcorrSHUFF = cellfun(@(x) corrcoef([trunctrjsSHUFFSMOOTH{x(1)}' trunctrjsSHUFFSMOOTH{x(2)}'], 'rows', 'pairwise'), permtrajsSHUFFs,'UniformOutput',false);
            %             PEcorrSHUFFmax = max(cell2mat(cellfun(@(x) x(1,2), PEcorrSHUFF,'UniformOutput',false)));
            %             PEcoef = PEcorrmax; % - PEcorrSHUFFmax;
            %             fldsepcomb{i,3} = PEcoef; %stick the PEC next to the traj data and indices
            %             PEcoefall(i,:) = [fldsepcomb{i,2} PEcoef];
            %         end
            
            

            
            %% calc Path Equivalent overlap coeficient.. this is based on the measure from battaglia 2004 and Singer et al 2010
            if runPECv2 && ~skipcell;
                PEtrajs = fldsepcombNOTSMOOTH{i,1}; % don't overwrite the flds data
                minlentraj = mat2cell(repmat(min(cell2mat(cellfun(@(x) length(x), PEtrajs,'UniformOutput',false))),1,length(PEtrajs)),1,length(PEtrajs)); %shortest trajectory
                trunctrjs = cellfun(@(x,y) x(1:y), PEtrajs, minlentraj, 'UniformOutput', false); %truncated trajs according to the shortest trajectory i'm clipping from the end because i haven't flipped the inbound trajs yet, so these should all be in terms of distance from center well.
                if length(trunctrjs) == 4; %i'm not sure which traj is which if there's less than 4 right now.. but i think this data is somewhere around here.. to do
                    trunctrjs{2} = fliplr(trunctrjs{2}); trunctrjs{4} = fliplr(trunctrjs{4});
                else
                    keyboard
                end
                gaussfilter = fspecial('gaussian',[1 6],1); %create 6 point gaussian filter with std 1
                trunctrjsSMOOTH = cellfun(@(x) conv(x,gaussfilter,'same'), trunctrjs, 'UniformOutput', false); % smooth
                %what the next 4 lines are doing is normalizing the smoothoccnormFR by the total smoothoccnormFR for all trajs for which the total smoothoccnormFR is not zero
                trunctrjsexclude = cell2mat(cellfun(@(x) nansum(x) ~= 0, trunctrjsSMOOTH, 'UniformOutput', false)); % logical indices of all non zero trajectories.. need this to exclude empty trajs from normalization which turns everything to nan because it divides by 0
                
                trunctrjsSMOOTHnormtemp = cellfun(@(x) x./nansum(x), trunctrjsSMOOTH(trunctrjsexclude), 'UniformOutput', false); % normalize the firing rates
                trunctrjsSMOOTHnorm = trunctrjsSMOOTH;
                trunctrjsSMOOTHnorm(trunctrjsexclude) = trunctrjsSMOOTHnormtemp; %only the total smoothoccnormFR that are not zero will be normalized
                
                permtraj = mat2cell([1 4; 2 3],2,1); %get overlap measure for the trajs in the same direction.. e.g. outleft(1) vs infromright(4)
                
                shufftrajs = cellfun(@(x) [x(floor(length(x)/2)+1:end) fliplr(x(1:floor(length(x)/2)))], trunctrjs,'UniformOutput',false); %A....BC....D becomes C....DB....A..
                %%trying a couple different ways to do a shuffle
                %             shiftmat = mat2cell([1:100],1,100);
                %             shufftrajs3 = cellfun(@(x) circshift(x,[0 y]), trunctrjs{1}, shiftmat, 'UniformOutput',false); %
                %           shufftrajs2 = cellfun(@(x) repmat(nanmean(x),1,length(x)), trunctrjs,'UniformOutput',false);
                shufftrajssmooth = cellfun(@(x) conv(x,gaussfilter,'same'), shufftrajs, 'UniformOutput', false); % smooth
                shufftrajssmoothnormtemp = cellfun(@(x) x./nansum(x), shufftrajssmooth(trunctrjsexclude), 'UniformOutput', false); % normalize the firing rates
                shufftrajssmoothnorm = shufftrajssmooth;
                shufftrajssmoothnorm(trunctrjsexclude) = shufftrajssmoothnormtemp;
                %PE2 means the style of the second generation path equivalence measure using the min overlap
                shuffPE2corrnorm = cell2mat(cellfun(@(x) (2*nansum(nanmin([trunctrjsSMOOTHnorm{x(1)};shufftrajssmoothnorm{x(2)}])))/(nansum(nansum([trunctrjsSMOOTHnorm{x(1)};shufftrajssmoothnorm{x(2)}]))), permtraj, 'UniformOutput', false)); % (2 * min of traj vecs)/sum
                shuffPE2corrnormmin = cell2mat(cellfun(@(x) nanmin([trunctrjsSMOOTHnorm{x(1)};shufftrajssmoothnorm{x(2)}]), permtraj, 'UniformOutput', false)); %shuffled traj 1 vs real traj 4. shuffled traj 2 vs real traj 3.
                shuffPE2corr = cell2mat(cellfun(@(x) (2*nansum(nanmin([trunctrjsSMOOTH{x(1)};shufftrajssmooth{x(2)}])))/(nansum(nansum([trunctrjsSMOOTH{x(1)};shufftrajssmooth{x(2)}]))), permtraj, 'UniformOutput', false)); % (2 * min of traj vecs)/sum
                shuffPE2corrmin = cell2mat(cellfun(@(x) nanmin([trunctrjsSMOOTH{x(1)};shufftrajssmooth{x(2)}]), permtraj, 'UniformOutput', false)); %shuffled traj 1 vs real traj 4. shuffled traj 2 vs real traj 3.
                PE2corr = cell2mat(cellfun(@(x) (2*nansum(min([trunctrjsSMOOTH{x(1)};trunctrjsSMOOTH{x(2)}])))/(nansum(nansum([trunctrjsSMOOTH{x(1)};trunctrjsSMOOTH{x(2)}]))), permtraj, 'UniformOutput', false)); % (2 * min of traj vecs)/sum
                PE2corrnorm = cell2mat(cellfun(@(x) (2*nansum(nanmin([trunctrjsSMOOTHnorm{x(1)};trunctrjsSMOOTHnorm{x(2)}])))/(nansum(nansum([trunctrjsSMOOTHnorm{x(1)};trunctrjsSMOOTHnorm{x(2)}]))), permtraj, 'UniformOutput', false)); % (2 * min of traj vecs)/sum
                PE2corrnormmin = cell2mat(cellfun(@(x) nanmin([trunctrjsSMOOTHnorm{x(1)};trunctrjsSMOOTHnorm{x(2)}]), permtraj, 'UniformOutput', false)); % (2 * min of traj vecs)/sum
                PE2corrnormMinusShuff = arrayfun(@(x,y) x-y, PE2corrnorm, shuffPE2corrnorm);
                PE2corrMinusShuff = arrayfun(@(x,y) x-y, PE2corr, shuffPE2corr);
                %PE1 means the style of the first generation path equivalence measure using the corrcoef
                PE1corrnormtemp = cellfun(@(x) corrcoef([trunctrjsSMOOTHnorm{x(1)}' trunctrjsSMOOTHnorm{x(2)}'], 'rows', 'pairwise'), permtraj,'UniformOutput',false);
                PE1corrnorm = cellfun(@(x) x(1,2), PE1corrnormtemp, 'UniformOutput', false);
                PE1corrnorm = cell2mat(PE1corrnorm);
                
                PE1corrtemp = cellfun(@(x) corrcoef([trunctrjsSMOOTH{x(1)}' trunctrjsSMOOTH{x(2)}'], 'rows', 'pairwise'), permtraj,'UniformOutput',false);
                PE1corr = cellfun(@(x) x(1,2), PE1corrtemp, 'UniformOutput', false);
                PE1corr = cell2mat(PE1corr);
                %structure containing ind as well as a bunch of metrics and data attempted... i think we've settled on using the normalized second gen PE (without shuff), called PE2corrnorm
                PEcoefallv2{i,1} = fldsepcomb{i,2}; PEcoefallv2{i,2} = PE2corr; PEcoefallv2{i,3} = PE2corrnorm;
                PEcoefallv2{i,4} = fldsepcombNOTSMOOTH{i,3}; %ca1 = 1, pfc = 2.
                PEcoefallv2{i,5} = PE2corrnormmin; PEcoefallv2{i,6} = {trunctrjsSMOOTHnorm}; PEcoefallv2{i,7} = trunctrjsSMOOTH;
                PEcoefallv2{i,8} = shuffPE2corrnorm; PEcoefallv2{i,9} = shuffPE2corrnormmin; PEcoefallv2{i,10} = PE2corrnormMinusShuff; PEcoefallv2{i,11} = PE1corrnorm;
                PEcoefallv2{i,12} = PE2corrMinusShuff; PEcoefallv2{i,13} = PE1corr;
                PEcoefallv2INDS(i,[1:4]) = fldsepcomb{i,2};
                %                 PEcoefallv2INDS(i,5) = fldsepcombNOTSMOOTH{i,3};
            else
                PEcoefallv2{i,1} = fldsepcomb{i,2}; PEcoefallv2{i,2} = nan; PEcoefallv2{i,3} = nan;PEcoefallv2{i,4} = nan;
                PEcoefallv2INDS(i,:) = uniqmatallcells(i, :);
            end
                       %% runscriptsparsity
            if runscriptsparsity && ~skipcell;
                areatrajs = cell2mat(fldsepcomb{i,1}); %concat all trajs
                arealength = length(areatrajs(~isnan(areatrajs))); %length of occupied bins
                areathresh = nanmax(areatrajs)*.25; %1/4 of peak rate
                sparsityabsval = nansum(areatrajs > areathresh);
                sparsityfrac = nansum(areatrajs > areathresh)/arealength;
                sparsityALLdata{i,1} = fldsepcomb{i,2}; 
                sparsityALLdata{i,2} = [sparsityabsval sparsityfrac]; 
                sparsityALLdata{i,3} = fldsepcombNOTSMOOTH{i,3};
            end
            
        end %
    end
    %%
    if runscriptISI
        cnt = 0;
        for anims = 1:length(animals);
            for cinds = 1:length(spks(anims).output{1});
                cnt = cnt+1;
                spksdata(cnt,1).goodspikes = spks(anims).output{1}(cinds).goodspikes;
                spksdata(cnt,1).index = [anims spks(anims).output{1}(cinds).index]; %store the animal # in with the other indices. [an day ep tet cell]
                matspksinds(cnt,:) = [anims flds(anims).output{1}(cinds).index([1 2 3 4])]; %anim day ep tet cell
            end
        end
    end
    % end %if runscript
    
    %%
    if runripposscript
        %     load /data19/sjadhav/HPExpt/ProcessedDataDR/_modf_rippos_Jan26 %made with /home/droumis/MATLAB/DFSsj_HPexpt_getripalignspikingGRAllPosition_DR.m
        load /data19/sjadhav/HPExpt/ProcessedDataDR/AllAn_PFCCA1_ripplepos_DR_Jan26modf_rippos_Jan26.mat %modf
        cnt = 0;
        for anims = 1:length(animals);
            for cinds = 1:length(modf(anims).output{1});
                cnt = cnt+1;
                modfdata(cnt,1) = modf(anims).output{1}(cinds);
                modfdata(cnt,1).index = [anims modfdata(cnt,1).index]; %store the animal # in with the other indices. [an day ep tet cell]
                matmodfinds(cnt,:) = modfdata(cnt,1).index;
            end
        end
    end
    %% spatial pair corr
                %for each pair, find their flds data, run a spatial corr, and stick the result next to the rip corr r vals
            badlength = [];
            for i = 1:length(matpairs);
                [ ~, pfcmatch] = ismember(matpairs(i,[1 2 5 6]), uniqmatallcells, 'rows');
                [ ~, ca1match] = ismember(matpairs(i,[1 2 3 4]), uniqmatallcells, 'rows');
                pfctrajdata = fldsepcomb{pfcmatch, 1};
                ca1trajdata = fldsepcomb{ca1match, 1};
                pfcpeak(i,1) = fldspeak(pfcmatch);
                ca1peak(i,1) = fldspeak(ca1match);
                minlentrajpair = mat2cell(repmat(min(cell2mat(cellfun(@(x) length(x), [pfctrajdata; ca1trajdata],'UniformOutput',false)),[],1),2,1), 2,4); %find the min length of each traj across epochs
                truncepspair{i} = cellfun(@(x,y) x(1:y), [pfctrajdata; ca1trajdata], minlentrajpair, 'UniformOutput', false); %truncated epochs according to the shortest trajectory across epochs.\
                truncepspair{i}{1,2} = fliplr(truncepspair{i}{1,2}); truncepspair{i}{1,4} = fliplr(truncepspair{i}{1,4}); truncepspair{i}{2,2} = fliplr(truncepspair{i}{2,2}); truncepspair{i}{2,4} = fliplr(truncepspair{i}{2,4}); %flip the inbound trajectories so they start from outer arm... for plotting
                allpairdata{i,1} = cell2mat(truncepspair{i})'; %first cell of alldata are vert concatenated vectors of all traj firing of the two cells... first col is the pfc cell.
                corrc = corrcoef(allpairdata{i}, 'rows', 'pairwise'); %'pairwise' means that it doesn't use a row if nan in either row..
                matpairs(i, 9) = corrc(1,2);
                matpairs(i, 10) = fldspeak(pfcmatch);
                matpairs(i, 11) = fldspeak(ca1match);
                allpairdata{i,2} = {matpairs(i,:)}; %second cell of alldata is all the matpairs info for that pair
                allpairdata{i,3} = mat2cell(cumsum(cell2mat(minlentrajpair(1,:))),1,4); %third cell of alldata are the ending points of each traj.. used for plotting trajplots
                
                pfcmapmatches = find(ismember(matfldsinds, matpairs(i,[1 2 5 6]), 'rows')); %this is different from the use of ismember above.. the first argument is the entire mat of indices and the second is the one to be matched.. this gives me ALL the matching row vecs..which i need now bc im working with multiple epochs with the maps
                pfcmapdata = arrayfun(@(x) fldsdata(x).mapdata.smoothedspikerate, pfcmapmatches, 'UniformOutput', false); %get the map data from all epochs of this cell
                
                ca1mapmatches = find(ismember(matfldsinds, matpairs(i,[1 2 3 4]), 'rows')); %'find' just spits out the indices of the nonzero elements, i.e. the row matches
                ca1mapdata = arrayfun(@(x) fldsdata(x).mapdata.smoothedspikerate, ca1mapmatches, 'UniformOutput', false); %get the map data from all epochs of this cell
                
                %          alldata{i,4} = [{fldsdata(pfcmapmatches(1)).mapdata.smoothedspikerate}; {fldsdata(ca1mapmatches(1)).mapdata.smoothedspikerate}]; %pfc, ca1 map data for the first run epoch of the pair.. this really should be an average across epochs
                allpairdata{i,4} = [{pfcmapdata}; {ca1mapdata}]; %pfc, ca1 map data for ALL the epochs
            end
    
end %runanyscript
toc

%%
if plotanything;
    if savefigs
        mkdir(figdir,savefigfilename)
    end
    if plotCorrCoef;
        figure; hold on;
        subplot(1,2,1)
        scatter(matpairs(:,8), matpairs(:,9),'.k'); hold on;
        [b,bint,r,rint,stats] = regress(matpairs(:,9), [ones(size(matpairs(:,8)),1) matpairs(:,8)]);
        plot(min(matpairs(:,8)):.1:max(matpairs(:,8)), b(1)+b(2)*[min(matpairs(:,8)):.1:max(matpairs(:,8))])
        xlabel('rip corr r'); ylabel('spat corr r'); title(sprintf('r^2(%0.3f) p(%0.3f)', stats(1), stats(4))); set(gcf,'PaperPositionMode','auto'); set(gcf, 'Position', [200 200 500 400])
        
        subplot(1,2,2)
        possigrippairs = matpairs(matpairs(:,7)<0.05 & matpairs(:,8) > 0,9);
        negsigrippairs = matpairs(matpairs(:,7)<0.05 & matpairs(:,8) < 0,9);
        nonsigrippairs = matpairs(matpairs(:,7)>0.05, 9);
        bar([mean(negsigrippairs) mean(nonsigrippairs) mean(possigrippairs)]); hold on;
        errorbar2([1 2 3], [mean(negsigrippairs) mean(nonsigrippairs) mean(possigrippairs)],  [stderr(negsigrippairs) stderr(nonsigrippairs) stderr(possigrippairs)] , 0.3, 'k')
        
        if savefigs==1
            figfile = [figdir,savefigfilename,'/',sprintf('ScatterBars_%s', savefigfilename)];
            print('-dpdf', figfile, '-r300');
        end
        
        %stats
        [pKW table statsKW] = kruskalwallis([negsigrippairs' nonsigrippairs' possigrippairs'], [ones(1, length(negsigrippairs)) ones(1,length(nonsigrippairs)).*2 ones(1,length(possigrippairs)).*3]);
        [c, m, h, gnames] = multcompare(statsKW, 'estimate', 'kruskalwallis', 'ctype', 'hsd', 'display', 'on', 'alpha', 0.01); % change the alpha values to determine significance range.
        
        if pausefigs; keyboard; end %pause
    end
    %% plotPFCCA1Trajs not really using this anymore
    
    if plotPFCCA1Trajs;
        [B IX] = sort(matpairs(:,9)); %sort pairs by ascending spat corr
        alldatasort = allpairdata(IX,:); %sort pairs by ascending spat corr
        for i = 1:16; %length(alltrajdata);
            subplot(4,4,i); hold on;
            pfcdata = alldatasort{i}(:,1);
            ca1data = alldatasort{i}(:,2);
            plot(pfcdata./max(pfcdata), 'k'); hold on; %normalize each cell's firing rate from 0-1.
            plot(ca1data./max(ca1data), 'r');
            title({sprintf('Spat Corr(%0.3f) Rip Corr(%0.3f)',alldatasort{i,2}{1}(1,9), alldatasort{i,2}{1}(1,8)); sprintf('Indx pfc(%d %d %d %d) ca1(%d %d %d %d)', alldatasort{i,2}{1}([1 2 5 6]), alldatasort{i,2}{1}([1 2 3 4]))});
            set(gca,'xtick',[], 'xticklabel', [], 'ytick', [], 'yticklabel', []);
            if alldatasort{i,2}{1}(1,8) > 0 %if the rip corr is > 0, turn background green
                set(gca,'color',[ .3 .9 .7]);
            end
            cellfun(@(x) line('XData',[x x], 'YData', [0 1], 'LineStyle', '-', 'LineWidth', 2, 'Color',[.8 .8 .8]), alldatasort{i,3}, 'UniformOutput', false); % add traj edge lines to plot
            if pausefigs; keyboard; end %pause
            
            if savefigs
                figfile = [figdir,savefigfilename,'/',sprintf('pair%d', i)];
                print('-dpng', figfile, '-r300');
            end
        end
    end
    
    %% plotPFCCA1Maps with traj plot at bottom..
    if plotPFCCA1MapsTrajs
        for i = 1:length(allpairdata);
            close all;
            pfcmapdata = allpairdata{i,4}{1}{1}; %this map is currently only the first epoch.. but the rest of the data is in here
            ca1mapdata = allpairdata{i,4}{2}{1}; %this map is currently only the first epoch.. but the rest of the data is in here
            pfctrajdata = allpairdata{i,1}(:,1)./max(allpairdata{i,1}(:,1)); %normalize each cell's firing rate from 0-1.
            ca1trajdata = allpairdata{i,1}(:,2)./max(allpairdata{i,1}(:,2)); %normalize each cell's firing rate from 0-1.
            tmpfig = figure; hold on;
            subplot(3,4,[1 2 5 6]);
            imagesc(pfcmapdata); set(gca,'xtick',[], 'xticklabel', [], 'ytick', [], 'yticklabel', []);
            subplot(3,4,[3 4 7 8]);
            imagesc(ca1mapdata); set(gca,'xtick',[], 'xticklabel', [], 'ytick', [], 'yticklabel', []);
            subplot(3,4,[9:12]);
            plot(pfctrajdata,'k'); hold on;
            plot(ca1trajdata,'r'); set(gca,'xtick',[], 'xticklabel', [], 'ytick', [], 'yticklabel', []);
            title('Black: PFC                           Red: CA1')
            cellfun(@(x) line('XData',[x x], 'YData', [0 1], 'LineStyle', '-', 'LineWidth', 2, 'Color',[.6 .6 .6]), allpairdata{i,3}, 'UniformOutput', false); % add traj edge lines to plot
            supertitle({sprintf('Spat Corr(%0.3f) Rip Corr(%0.3f)',allpairdata{i,2}{1}(1,9), allpairdata{i,2}{1}(1,8)); sprintf('Indx pfc(%d %d %d %d) ca1(%d %d %d %d)', allpairdata{i,2}{1}([1 2 5 6]), allpairdata{i,2}{1}([1 2 3 4]))});
            
            if pausefigs; keyboard; end %pause
            if savefigs
                figfile = [figdir,savefigfilename,'/',sprintf('An%dD%dpair%d', allpairdata{i,2}{1}([1]), allpairdata{i,2}{1}([2]), i)];
                print('-dpng', figfile, '-r300');
            end
            close(tmpfig)
        end
    end
    
    %% Plot ripple excited vs ripple inhibited pfc cells not really using this
    if plotPFCripmodMapsTrajs
        %         load swrmodinds.mat
        load /home/droumis/MATLAB/swrmodinds_Jan27th.mat
        %        GIDEON FIXED THIS
        %         %patch to fix gideon's reordering animals. puts Borg at the end so that it's hpa hpb, hpc, nadal, rosenthal, borg. don't run this after he fixes the order
        %         PFCindsExc(PFCindsExc(:,1) == 4, 1) = 7;
        %         PFCindsExc(PFCindsExc(:,1) == 5, 1) = 4;
        %         PFCindsExc(PFCindsExc(:,1) == 6, 1) = 5;
        %         PFCindsExc(PFCindsExc(:,1) == 7, 1) = 6;
        %
        %         PFCindsInh(PFCindsInh(:,1) == 4, 1) = 7;
        %         PFCindsInh(PFCindsInh(:,1) == 5, 1) = 4;
        %         PFCindsInh(PFCindsInh(:,1) == 6, 1) = 5;
        %         PFCindsInh(PFCindsInh(:,1) == 7, 1) = 6;
        
        PFCindsmod = [PFCindsExc;PFCindsInh];
        
        for i = 1:length(PFCindsmod);
            close all;
            pfcmapmatches = find(ismember(matfldsinds, PFCindsmod(i,:), 'rows')); %find all the epochs for this pfc cell
            pfcmapdata = arrayfun(@(x) fldsdata(x).mapdata.smoothedspikerate, pfcmapmatches, 'UniformOutput', false); %get the map data from all epochs of this cell
            tmpfig = figure; hold on;
            for j = 1:length(pfcmapdata);
                subplot(1,length(pfcmapdata),j);
                imagesc(pfcmapdata{j}); set(gca,'xtick',[], 'xticklabel', [], 'ytick', [], 'yticklabel', []);
                %                 title(pfcmeanrate(i){j});
            end
            if i <= length(PFCindsExc);
                colr = [1 0 .6];
                supertitle({sprintf('Rip Exc (%d %d %d %d)', PFCindsmod(i,:))}, colr); % I added a color var to supertitle as the second input
            else colr = [0 0 1];
                supertitle({sprintf('Rip Inh (%d %d %d %d)', PFCindsmod(i,:))}, colr); % I added a color var to supertitle as the second input
            end %magenta title if exc, blue if inh
            
            if pausefigs; keyboard; end %pause
            if savefigs
                %                 figfile = [figdir,savefigfilename,'/',sprintf('An%dD%dT%dC%d', PFCindsmod(i,1), PFCindsmod(i, 2), PFCindsmod(i, 3), PFCindsmod(i, 4))];
                figfile = [figdir,savefigfilename,'/',sprintf('pfcripmod%d', i)];
                print('-dpng', figfile, '-r300');
            end
            close(tmpfig)
        end
    end
    
    %% plotPFCripmodMeanRates
    
    if plotPFCripmodMeanRates
        close all;
        %         load swrmodinds.mat
        load /home/droumis/MATLAB/swrmodinds_Jan27th.mat
        %        GIDEON FIXED THIS
        %         %patch to fix gideon's reordering animals. puts Borg at the end so that it's hpa hpb, hpc, nadal, rosenthal, borg. don't run this after he fixes the order
        %         PFCindsExc(PFCindsExc(:,1) == 4, 1) = 7;
        %         PFCindsExc(PFCindsExc(:,1) == 5, 1) = 4;
        %         PFCindsExc(PFCindsExc(:,1) == 6, 1) = 5;
        %         PFCindsExc(PFCindsExc(:,1) == 7, 1) = 6;
        %
        %         PFCindsInh(PFCindsInh(:,1) == 4, 1) = 7;
        %         PFCindsInh(PFCindsInh(:,1) == 5, 1) = 4;
        %         PFCindsInh(PFCindsInh(:,1) == 6, 1) = 5;
        %         PFCindsInh(PFCindsInh(:,1) == 7, 1) = 6;
        
        for i = 1:length(PFCindsExc);
            pfcexcmatches = find(ismember(matfldsinds, PFCindsExc(i,:), 'rows')); %find all the epochs for this pfc cell
            pfcEXCmeanrate(i) = mean(cell2mat(arrayfun(@(x) mean(fldsdata(x).mapdata.spikerate(find(fldsdata(x).mapdata.occupancy))), pfcexcmatches, 'UniformOutput', false))); % mean the spikerates for all bins where occupancy is not zero.. then mean across all epochs
        end
        for i = 1:length(PFCindsInh);
            pfcinhmatches = find(ismember(matfldsinds, PFCindsInh(i,:), 'rows')); %find all the epochs for this pfc cell
            pfcINHmeanrate(i) = mean(cell2mat(arrayfun(@(x) mean(fldsdata(x).mapdata.spikerate(find(fldsdata(x).mapdata.occupancy))), pfcinhmatches, 'UniformOutput', false))); % mean the spikerates for all bins where occupancy is not zero.. then mean across all epochs
        end
        
        bar([mean(pfcINHmeanrate) mean(pfcEXCmeanrate)]); hold on;
        errorbar2([1 2], [mean(pfcINHmeanrate) mean(pfcEXCmeanrate)],  [stderr(pfcINHmeanrate) stderr(pfcEXCmeanrate)] , 0.3, 'k')
        
        if savefigs
            figfile = [figdir,savefigfilename,'/',sprintf('pfcripmodRATES%d', i)];
            print('-dpng', figfile, '-r300');
        end
        
        [h,p,ci,stats] = ttest2(pfcINHmeanrate, pfcEXCmeanrate);
        p
        %         [pKW table statsKW] = kruskalwallis([pfcINHmeanrate pfcEXCmeanrate], [ones(1, length(pfcINHmeanrate)) ones(1,length(pfcEXCmeanrate)).*2]);
        %         [c, m, h, gnames] = multcompare(statsKW, 'estimate', 'kruskalwallis', 'ctype', 'hsd', 'display', 'on', 'alpha', 0.01); % change the alpha values to determine significance range.
        %
        if pausefigs; keyboard; end %pause
    end
    
    
    
    
    %% plotPFCripmodISI  not rreally using this
    plotPFCripmodISI = 0;
    if plotPFCripmodISI
        %         %collect ALL of the spks data into workable format
        %         cnt = 0;
        %         for anims = 1:length(animals);
        %             for cinds = 1:length(spks(anims).output{1});
        %                 cnt = cnt+1;
        %                 spksdata(cnt,1).goodspikes = spks(anims).output{1}(cinds).goodspikes;
        %                 spksdata(cnt,1).index = [anims spks(anims).output{1}(cinds).index]; %store the animal # in with the other indices. [an day ep tet cell]
        %                 matspksinds(cnt,:) = [anims flds(anims).output{1}(cinds).index([1 2 3 4])]; %anim day ep tet cell
        %             end
        %         end
        close all;
        %         load swrmodinds.mat
        load /home/droumis/MATLAB/swrmodinds_Jan27th.mat
        %        GIDEON FIXED THIS
        %         %patch to fix gideon's reordering animals. puts Borg at the end so that it's hpa hpb, hpc, nadal, rosenthal, borg. don't run this after he fixes the order
        %         PFCindsExc(PFCindsExc(:,1) == 4, 1) = 7;
        %         PFCindsExc(PFCindsExc(:,1) == 5, 1) = 4;
        %         PFCindsExc(PFCindsExc(:,1) == 6, 1) = 5;
        %         PFCindsExc(PFCindsExc(:,1) == 7, 1) = 6;
        %
        %         PFCindsInh(PFCindsInh(:,1) == 4, 1) = 7;
        %         PFCindsInh(PFCindsInh(:,1) == 5, 1) = 4;
        %         PFCindsInh(PFCindsInh(:,1) == 6, 1) = 5;
        %         PFCindsInh(PFCindsInh(:,1) == 7, 1) = 6;
        
        PFCindsmod = [PFCindsExc;PFCindsInh];
        
        for i = 1:length(PFCindsmod);
            close all;
            pfcISImatches = find(ismember(matspksinds(:,[1 2 4 5]), PFCindsmod(i,:), 'rows')); %find all the epochs for this pfc cell
            pfcISIdata = arrayfun(@(x) diff(spksdata(x).goodspikes(:,1)), pfcISImatches, 'UniformOutput', false); % get isi for each cell's epochs
            tmpfig = figure; hold on;
            for j = 1:length(pfcISIdata);
                subplot(1,length(pfcISIdata),j);
                bins = linspace(min(pfcISIdata{j}),max(pfcISIdata{j}),length(pfcISIdata{j}));
                n = hist(pfcISIdata{j},length(pfcISIdata{j}));
                lower95 = length(n(cumsum(n)<(length(pfcISIdata{j})*0.90))); %only use the lower 95% of the isi bins.. i.e. discard the long tail that makes it hard to visualize
                bar(bins(1:lower95), n(1:lower95),'EdgeColor','none');
                %                 axis([0 bins(lowemean(pfcINHPE)r95) 0 max(n)]);
                xlabel('ISI(s)'); ylabel('count');
            end
            if i <= length(PFCindsExc);
                colr = [1 0 .6];
                supertitle({sprintf('Rip Exc (%d %d %d %d)', PFCindsmod(i,:))}, colr); % I added a color var to supertitle as the second input
            else colr = [0 0 1];
                supertitle({sprintf('Rip Inh (%d %d %d %d)', PFCindsmod(i,:))}, colr); % I added a color var to supertitle as the second input
            end %magenta title if exc, blue if inh
            
            if pausefigs; keyboard; end %pause
            if savefigs
                %                 figfile = [figdir,savefigfilename,'/',sprintf('An%dD%dT%dC%d', PFCindsmod(i,1), PFCindsmod(i, 2), PFCindsmod(i, 3), PFCindsmod(i, 4))];
                figfile = [figdir,savefigfilename,'/',sprintf('pfcripmodISI%d', i)];
                print('-dpng', figfile, '-r300');
            end
            close(tmpfig)
        end
    end
    
    %% plotPECall ... the dayswise stuff is currently incorrect because nadal starts at 8 and also rosenthal goes from 1-12
    
    if plotPECall
        %collect ALL of the spks data into workable format
        %         load swrmodinds.mat
        
        load /home/droumis/MATLAB/swrmodinds_Jan27th.mat
        skippedPE = []; clear pfcEXCPE pfcINHPE
        usePECval = 3; %3 is norm PE2, 11 is norm PE1
        
        for i = 1:length(PFCindsExc);
            try%skip cells without PEC
                %              mean(pfcINHPE)   pfcexcmatch = find(ismember(PEcoefall(:,1:4), PFCindsExc(i,:), 'rows'));
                pfcexcmatch = find(ismember(PEcoefallv2INDS(:,1:4), PFCindsExc(i,:), 'rows'));
                %                 pfcEXCPE(i,:) = PEcoefallv2{pfcexcmatch,3}';
                pfcEXCPE(i,:) = [PEcoefallv2{pfcexcmatch,usePECval}' PEcoefallv2{pfcexcmatch,1}];
                if PEcoefallv2{pfcexcmatch,1}(1,[1:4]) ~= sparsityALLdata{pfcexcmatch,1};
                    keyboard %pause if the sparsity data and the PE data don't have matching indices
                end
                pfcEXCSpar(i,:) = [sparsityALLdata{pfcexcmatch,2} nanmean(pfcEXCPE(i,[1:2]))];
            catch
                skippedPE = [skippedPE; PFCindsExc(i,:)];
            end
        end
        
        
        for i = 1:length(PFCindsInh);
            try %skip cells without PEC
                %                 pfcinhmatch = find(ismember(PEcoefall(:,1:4), PFCindsInh(i,:), 'rows'));
                %                 pfcINHPE(i) = PEcoefall(pfcinhmatch,5);
                pfcinhmatch = find(ismember(PEcoefallv2INDS(:,1:4), PFCindsInh(i,:), 'rows'));
                %                 pfcINHPE(i,:) = PEcoefallv2{pfcinhmatch,3}';
                pfcINHPE(i,:) = [PEcoefallv2{pfcinhmatch,usePECval}' PEcoefallv2{pfcinhmatch,1}];
                if PEcoefallv2{pfcinhmatch,1}(1,[1:4]) ~= sparsityALLdata{pfcinhmatch,1};
                    keyboard %pause if the sparsity data and the PE data don't have matching indices
                end
                pfcINHSpar(i,:) = [sparsityALLdata{pfcinhmatch,2} nanmean(pfcINHPE(i,[1:2]))];
            catch
                skippedPE = [skippedPE; PFCindsInh(i,:)];
            end
        end
        
        %collect all the PE values [val1 val2 animal day area]
        for i = 1:length(PEcoefallv2)
            try
%                 if i == 984; keyboard; end
                PEall(i,[1:2]) = PEcoefallv2{i,usePECval}'; %PE values
                PEall(i,[3:6]) = PEcoefallv2{i,1}(1,[1:4]); %PEval1 PEval2 animal day tet cell
                PEall(i,7) = PEcoefallv2{i,4}(1); %pfc = 2 ca = 1; %PEval1 PEval2 animal day tet cell PFC(2)orCA1(1)
                if PEcoefallv2{i,1}(1,[1:4]) ~= sparsityALLdata{i,1};
                    keyboard %pause if the sparsity data and the PE data don't have matching indices
                end
                SparAll(i,:) = [sparsityALLdata{i,1} sparsityALLdata{i,3} sparsityALLdata{i,2} nanmean(PEall(i,[1:2]))]; %sparsity values and average of PE
            catch
            end
        end
        % sort for the days before reshaping
        for i = 1:8;
            pfcPEdays{i} = PEall(PEall(:,4) == i & PEall(:,7) == 2,[1:2]);%& PEall(:,3) < 4
            pfcPEdays{i} = reshape(pfcPEdays{i},2*length(pfcPEdays{i}(:,1)),1);
            ca1PEdays{i} = PEall(PEall(:,4) == i & PEall(:,7) == 1,[1:2]); %& PEall(:,3) < 4
            ca1PEdays{i} = reshape(ca1PEdays{i},2*length(ca1PEdays{i}(:,1)),1);
            pfcEXCPEdays{i} = pfcEXCPE(pfcEXCPE(:,4) == i & pfcEXCPE(:,3) < 4,[1:2]);
            pfcEXCPEdays{i} = reshape(pfcEXCPEdays{i},2*length(pfcEXCPEdays{i}(:,1)),1);
            pfcINHPEdays{i} = pfcINHPE(pfcINHPE(:,4) == i & pfcINHPE(:,3) < 4,[1:2]);
            pfcINHPEdays{i} = reshape(pfcINHPEdays{i},2*length(pfcINHPEdays{i}(:,1)),1);
        end
        
        %reshape ripmodpfc
        pfcEXCPE = reshape(pfcEXCPE(:,[1:2]), 2*length(pfcEXCPE(:,1)),1);
        pfcINHPE = reshape(pfcINHPE(:,[1:2]),2*length(pfcINHPE(:,1)),1);
        
        %separate ca1 and pfc(excluding the ripmod pfcs)
        ca1PEall = PEall(PEall(:,5) == 1, [1:2]);
        ca1PEall = reshape(ca1PEall, 2*length(ca1PEall),1);
        pfcPEall = PEall(PEall(:,5) == 2, [1:6]);
        ripmodPFCinds = [PFCindsExc;PFCindsInh];
        pfcPEallexcluded = pfcPEall(find(~ismember(pfcPEall(:,[3:6]), ripmodPFCinds,'rows')),[1:2]); %exclude the pfc cells that are in the rip mod list
        pfcPEallexcluded = reshape(pfcPEallexcluded, 2*length(pfcPEallexcluded),1);
        
        ca1Sparall = SparAll(SparAll(:,5) == 1, [6:8]); %abs spar val, % spar val , mean PE val
        pfcSparall = SparAll(SparAll(:,5) == 2, [1:8]);
        pfcSparallexcl = pfcSparall(find(~ismember(pfcSparall(:,[1:4]), ripmodPFCinds,'rows')),[6:8]); %exclude the pfc cells that are in the rip mod list
        pfcSparallexclfrac = pfcSparallexcl(:,2); %fraction of area
        pfcSparallexclabs = pfcSparallexcl(:,1); %abs area
        
        %create histograms and normalize
        edges = [0:.1:1];
        pfcINHPEhist = histc(pfcINHPE,edges); pfcINHPEhist = pfcINHPEhist./nansum(pfcINHPEhist);
        pfcEXCPEhist = histc(pfcEXCPE,edges); pfcEXCPEhist = pfcEXCPEhist./nansum(pfcEXCPEhist);
        pfcPEallexcludedhist = histc(pfcPEallexcluded,edges); pfcPEallexcludedhist = pfcPEallexcludedhist./nansum(pfcPEallexcludedhist);
        ca1PEallhist = histc(ca1PEall,edges); ca1PEallhist = ca1PEallhist./nansum(ca1PEallhist);
        
        %PLOT
        figure; subplot(2,2,1);
        bar(1,nanmean(pfcINHPE),'facecolor','b'); hold on;
        bar(2,nanmean(pfcEXCPE),'facecolor','r');
        bar(3,nanmean(pfcPEallexcluded),'facecolor','g');
        bar(4,nanmean(ca1PEall),'facecolor','k');
        title('Inh Exc PFC CA1'); set(gca,'xtick',[], 'xticklabel', [])
        errorbar2([1 2 3 4], [nanmean(pfcINHPE) nanmean(pfcEXCPE) nanmean(pfcPEallexcluded) nanmean(ca1PEall)],  [stderr(pfcINHPE) stderr(pfcEXCPE) stderr(pfcPEallexcluded(~isnan(pfcPEallexcluded))) stderr(ca1PEall(~isnan(ca1PEall)))] , 0.3, 'k')
        
        subplot(2,2,2);
        %         histPE = bar(edges,[pfcINHPEhist pfcEXCPEhist pfcPEallexcludedhist ca1PEallhist]);
        %         set(histPE(1),'FaceColor','b'); set(histPE(2),'FaceColor','r'); set(histPE(3),'FaceColor','g'); set(histPE(4),'FaceColor','k');
        plot(edges,pfcINHPEhist, 'b',edges,pfcEXCPEhist,'r', edges, pfcPEallexcludedhist,'g', edges, ca1PEallhist,'k');
        axis([-.1 1 0 max(max([pfcINHPEhist pfcEXCPEhist pfcPEallexcludedhist ca1PEallhist]))]);
        xlabel('norm Overlap'); ylabel('% pairs');

        
        subplot(2,2,3)
        bar(1,nanmean(pfcINHSpar(:,2)),'facecolor','b'); hold on;
        bar(2,nanmean(pfcEXCSpar(:,2)),'facecolor','r');
        bar(3,nanmean(pfcSparallexcl(:,2)),'facecolor','g');
        bar(4,nanmean(ca1Sparall(:,2)),'facecolor','k');
        title('% Area >.25peak'); set(gca,'xtick',[], 'xticklabel', [])
        errorbar2([1 2 3 4], [nanmean(pfcINHSpar(:,2)) nanmean(pfcEXCSpar(:,2)) nanmean(pfcSparallexcl(:,2)) nanmean(ca1Sparall(:,2))],  [stderr(pfcINHSpar(:,2)) stderr(pfcEXCSpar(:,2)) stderr(pfcSparallexcl(~isnan(pfcSparallexcl(:,2)),2)) stderr(ca1Sparall(:,2))] , 0.3, 'k')

                subplot(2,2,4)
        bar(1,nanmean(pfcINHSpar(:,1)),'facecolor','b'); hold on;
        bar(2,nanmean(pfcEXCSpar(:,1)),'facecolor','r');
        bar(3,nanmean(pfcSparallexcl(:,1)),'facecolor','g');
        bar(4,nanmean(ca1Sparall(:,1)),'facecolor','k');
        title('abs Area >.25peak'); set(gca,'xtick',[], 'xticklabel', [])
        errorbar2([1 2 3 4], [nanmean(pfcINHSpar(:,1)) nanmean(pfcEXCSpar(:,1)) nanmean(pfcSparallexcl(:,1)) nanmean(ca1Sparall(:,1))],  [stderr(pfcINHSpar(:,1)) stderr(pfcEXCSpar(:,1)) stderr(pfcSparallexcl(~isnan(pfcSparallexcl(:,1)),1)) stderr(ca1Sparall(:,1))] , 0.3, 'k')

        
        % scatter plots of PE2 vs Sparsity
        figure
        subplot(2,2,1)
        scatter(pfcINHSpar(:,2),pfcINHSpar(:,3), '.b'); hold on;
         [b,bint,r,rint,stats] = regress(pfcINHSpar(:,3), [ones(size(pfcINHSpar(:,2)),1) pfcINHSpar(:,2)]);
         plot(min(pfcINHSpar(:,3)):.01:max(pfcINHSpar(:,3)), b(1)+b(2)*[min(pfcINHSpar(:,3)):.01:max(pfcINHSpar(:,3))],'b')
         axis([0 1 0 1]);
        xlabel('%spacecover'); ylabel('PEoverlap'); title(sprintf('r^2(%0.3f) p(%0.3f)', stats(1), stats(4)));
        
        subplot(2,2,2)
        scatter(pfcEXCSpar(:,2),pfcEXCSpar(:,3), '.r'); hold on;
                 [b,bint,r,rint,stats] = regress(pfcEXCSpar(:,3), [ones(size(pfcEXCSpar(:,2)),1) pfcEXCSpar(:,2)]);
         plot(min(pfcEXCSpar(:,3)):.01:max(pfcEXCSpar(:,3)), b(1)+b(2)*[min(pfcEXCSpar(:,3)):.01:max(pfcEXCSpar(:,3))],'r')
         axis([0 1 0 1]);
        xlabel('%spacecover'); ylabel('PEoverlap'); title(sprintf('r^2(%0.3f) p(%0.3f)', stats(1), stats(4)));
        
        subplot(2,2,3)
        scatter(pfcSparallexcl(:,2),pfcSparallexcl(:,3),'.g'); hold on;
                 [b,bint,r,rint,stats] = regress(pfcSparallexcl(:,3), [ones(size(pfcSparallexcl(:,2)),1) pfcSparallexcl(:,2)]);
         plot(min(pfcSparallexcl(:,3)):.01:max(pfcSparallexcl(:,3)), b(1)+b(2)*[min(pfcSparallexcl(:,3)):.01:max(pfcSparallexcl(:,3))],'g')
         axis([0 1 0 1]);
        xlabel('%spacecover'); ylabel('PEoverlap'); title(sprintf('r^2(%0.3f) p(%0.3f)', stats(1), stats(4)));
        
        subplot(2,2,4)
        scatter(ca1Sparall(:,2), ca1Sparall(:,3), '.k'); hold on;
                 [b,bint,r,rint,stats] = regress(ca1Sparall(:,3), [ones(size(ca1Sparall(:,2)),1) ca1Sparall(:,2)]);
         plot(min(ca1Sparall(:,3)):.01:max(ca1Sparall(:,3)), b(1)+b(2)*[min(ca1Sparall(:,3)):.01:max(ca1Sparall(:,3))],'k')
         axis([0 1 0 1]);
        xlabel('%spacecover'); ylabel('PEoverlap'); title(sprintf('r^2(%0.3f) p(%0.3f)', stats(1), stats(4)));
        
% % -------- plotting across days
%         subplot(3,2,3)
%         bar([nanmean(pfcINHPEdays{1}) nanmean(pfcINHPEdays{2}) nanmean(pfcINHPEdays{3}) nanmean(pfcINHPEdays{4}) nanmean(pfcINHPEdays{5}) nanmean(pfcINHPEdays{6}) nanmean(pfcINHPEdays{7}) nanmean(pfcINHPEdays{8})],'facecolor', 'b'); hold on;
%         errorbar2([1 2 3 4 5 6 7 8], [nanmean(pfcINHPEdays{1}) nanmean(pfcINHPEdays{2}) nanmean(pfcINHPEdays{3}) nanmean(pfcINHPEdays{4}) nanmean(pfcINHPEdays{5}) nanmean(pfcINHPEdays{6}) nanmean(pfcINHPEdays{7}) nanmean(pfcINHPEdays{8})],...
%             [stderr(pfcINHPEdays{1}) stderr(pfcINHPEdays{2}) stderr(pfcINHPEdays{3}) stderr(pfcINHPEdays{4}) stderr(pfcINHPEdays{5}) stderr(pfcINHPEdays{6}) stderr(pfcINHPEdays{7}) stderr(pfcINHPEdays{8})],0.3,'k');
%         ylabel('norm Overlap'); xlabel('day');
%         
%         subplot(3,2,4)
%         bar([nanmean(pfcEXCPEdays{1}) nanmean(pfcEXCPEdays{2}) nanmean(pfcEXCPEdays{3}) nanmean(pfcEXCPEdays{4}) nanmean(pfcEXCPEdays{5}) nanmean(pfcEXCPEdays{6}) nanmean(pfcEXCPEdays{7}) nanmean(pfcEXCPEdays{8})],'facecolor', 'r'); hold on;
%         errorbar2([1 2 3 4 5 6 7 8], [nanmean(pfcEXCPEdays{1}) nanmean(pfcEXCPEdays{2}) nanmean(pfcEXCPEdays{3}) nanmean(pfcEXCPEdays{4}) nanmean(pfcEXCPEdays{5}) nanmean(pfcEXCPEdays{6}) nanmean(pfcEXCPEdays{7}) nanmean(pfcEXCPEdays{8})],...
%             [stderr(pfcEXCPEdays{1}) stderr(pfcEXCPEdays{2}) stderr(pfcEXCPEdays{3}) stderr(pfcEXCPEdays{4}) stderr(pfcEXCPEdays{5}) stderr(pfcEXCPEdays{6}) stderr(pfcEXCPEdays{7}) stderr(pfcEXCPEdays{8})],0.3,'k');
%         
%         
%         subplot(3,2,5);
%         bar([nanmean(pfcPEdays{1}) nanmean(pfcPEdays{2}) nanmean(pfcPEdays{3}) nanmean(pfcPEdays{4}) nanmean(pfcPEdays{5}) nanmean(pfcPEdays{6}) nanmean(pfcPEdays{7}) nanmean(pfcPEdays{8})],'facecolor', 'g'); hold on;
%         errorbar2([1 2 3 4 5 6 7 8], [nanmean(pfcPEdays{1}) nanmean(pfcPEdays{2}) nanmean(pfcPEdays{3}) nanmean(pfcPEdays{4}) nanmean(pfcPEdays{5}) nanmean(pfcPEdays{6}) nanmean(pfcPEdays{7}) nanmean(pfcPEdays{8})],...
%             [stderr(pfcPEdays{1}) stderr(pfcPEdays{2}) stderr(pfcPEdays{3}) stderr(pfcPEdays{4}) stderr(pfcPEdays{5}) stderr(pfcPEdays{6}) stderr(pfcPEdays{7}) stderr(pfcPEdays{8})],0.3,'k');
%         
%         subplot(3,2,6);
%         bar([nanmean(ca1PEdays{1}) nanmean(ca1PEdays{2}) nanmean(ca1PEdays{3}) nanmean(ca1PEdays{4}) nanmean(ca1PEdays{5}) nanmean(ca1PEdays{6}) nanmean(ca1PEdays{7}) nanmean(ca1PEdays{8})],'facecolor','k'); hold on;
%         errorbar2([1 2 3 4 5 6 7 8], [nanmean(ca1PEdays{1}) nanmean(ca1PEdays{2}) nanmean(ca1PEdays{3}) nanmean(ca1PEdays{4}) nanmean(ca1PEdays{5}) nanmean(ca1PEdays{6}) nanmean(ca1PEdays{7}) nanmean(ca1PEdays{8})],...
%             [stderr(ca1PEdays{1}) stderr(ca1PEdays{2}) stderr(ca1PEdays{3}) stderr(ca1PEdays{4}) stderr(ca1PEdays{5}) stderr(ca1PEdays{6}) stderr(ca1PEdays{7}) stderr(ca1PEdays{8})],0.3,'b');
% % ---------------        
        
        %         hist(pfcINHPE)
        %         figure
        %         hist(pfcEXCPE)
        
        %%also need to plot the PE values for CA1 cells, as well as all the rest of the pFC cells to compare..
        
        if savefigs
            figfile = [figdir,savefigfilename,'/',sprintf('pfcripmodPEC_norm')];
            print('-dpng', figfile, '-r300');
        end
        [pKW table statsKW] = kruskalwallis([pfcINHPE' pfcEXCPE' pfcPEallexcluded' ca1PEall'], [ones(1, length(pfcINHPE)) ones(1,length(pfcEXCPE)).*2 ones(1,length(pfcPEallexcluded)).*3 ones(1,length(ca1PEall)).*4]);
        [c, m, h, gnames] = multcompare(statsKW, 'estimate', 'kruskalwallis', 'ctype', 'hsd', 'display', 'on', 'alpha', 0.01); % change the alpha values to determine significance range.
        if pausefigs; keyboard; end %pause
        [pKW table statsKW] = kruskalwallis([pfcINHSpar(:,2)' pfcEXCSpar(:,2)' pfcSparallexcl(:,2)' ca1Sparall(:,2)'], [ones(1, length(pfcINHSpar(:,2))) ones(1,length(pfcEXCSpar(:,2))).*2 ones(1,length(pfcSparallexcl(:,2))).*3 ones(1,length(ca1Sparall(:,2))).*4]);
        [c, m, h, gnames] = multcompare(statsKW, 'estimate', 'kruskalwallis', 'ctype', 'hsd', 'display', 'on', 'alpha', 0.05); % change the alpha values to determine significance range.
        if pausefigs; keyboard; end %pause
        [pKW table statsKW] = kruskalwallis([pfcINHSpar(:,1)' pfcEXCSpar(:,1)' pfcSparallexcl(:,1)' ca1Sparall(:,1)'], [ones(1, length(pfcINHSpar(:,1))) ones(1,length(pfcEXCSpar(:,1))).*2 ones(1,length(pfcSparallexcl(:,1))).*3 ones(1,length(ca1Sparall(:,1))).*4]);
        [c, m, h, gnames] = multcompare(statsKW, 'estimate', 'kruskalwallis', 'ctype', 'hsd', 'display', 'on', 'alpha', 0.05); % change the alpha values to determine significance range.
        
        
        %         [h,p,ci,stats] = ttesset(gca,'xtick',[], 'xticklabel', [],t2(pfcINHPE, pfcEXCPE);
        %         p
        
    end
    
    %% plotPFCripmodMapsTrajsISImeanrate
    %add the
    
    if plotPFCripmodMapsTrajsISImeanrate
        load /data19/sjadhav/HPExpt/ProcessedDataDR/_modf_rippos_Jan26 %,made with /home/droumis/MATLAB/DFSsj_HPexpt_getripalignspikingGRAllPosition_DR.m
        cnt = 0;
        for anims = 1:length(animals);
            for cinds = 1:length(flds(anims).output{1});
                cnt = cnt+1;
                fldsdata(cnt,1) = flds(anims).output{1}(cinds);
                fldsdata(cnt,1).index = [anims fldsdata(cnt,1).index]; %store the animal # in with the other indices. [an day ep tet cell]
                matfldsinds(cnt,:) = [anims flds(anims).output{1}(cinds).index([1 3 4])]; %anim day tet cell
            end
        end
        %         load swrmodinds.mat
        load /home/droumis/MATLAB/swrmodinds_Jan27th.mat
        
        PFCindsmod = [PFCindsExc;PFCindsInh];
        trajsnotflipped = [];
        not4trajs = [];
        for i = 1:length(PFCindsmod);
            close all;
            pfcmapmatches = find(ismember(matfldsinds, PFCindsmod(i,:), 'rows')); %find all the epochs for this pfc cell
            pfcmapdata = arrayfun(@(x) fldsdata(x).mapdata.smoothedspikerate, pfcmapmatches, 'UniformOutput', false); %get the map data from all epochs of this cell
            %             pfcISImatches = find(ismember(matspksinds(:,[1 2 4 5]), PFCindsmod(i,:), 'rows')); %find all the epochs for this pfc cell
            %             pfcISIdata = arrayfun(@(x) diff(spksdata(x).goodspikes(:,1)), pfcISImatches, 'UniformOutput', false); % get isi for each cell's epochs
            pfcmeanrates = cell2mat(arrayfun(@(x) mean(fldsdata(x).mapdata.spikerate(find(fldsdata(x).mapdata.occupancy))), pfcmapmatches, 'UniformOutput', false)); % mean the spikerates for all bins where occupancy is not zero..
            %             pfctrajdata = arrayfun(@(x) (fldsdata(x).trajdata), pfcmapmatches, 'UniformOutput', false); %get the traj data from all epochs of this cell
            pfcPEmatch = find(ismember(PEcoefallv2INDS(:,1:4), PFCindsmod(i,:), 'rows')); %find all the epochs for this pfc cell
            pfcRipPosMatches = find(ismember(matmodfinds(:,[1 2 4 5]), PFCindsmod(i,:), 'rows')); %find all the epochs for this pfc cell
            pfcRipPosdata = arrayfun(@(x) modfdata(x).ripPos, pfcRipPosMatches, 'UniformOutput', false); %get the rip pos data from all epochs of this cell
            pfcRipPosSpikedata = arrayfun(@(x) modfdata(x).trialResps, pfcRipPosMatches, 'UniformOutput', false); %get the rip pos data from all epochs of this cell
            
            %previously was using traj data directly from fldsdata
            %             for j = 1:length(pfctrajdata);
            %                 k = pfctrajdata{j};
            %                 try %error if missing a traj
            %                     pfctrajdata{j}{1,2} = flipud(pfctrajdata{j}{1,2});
            %                     pfctrajdata{j}{1,4} = flipud(pfctrajdata{j}{1,4});
            %                 catch
            %                     trajsnotflipped = [trajsnotflipped; PFCindsmod(i,:)];
            %                 end
            %
            %                 for l = 1:length(k);
            %                     trajlengths{j}(l) = length(k{l});
            %                 end
            %                 trajlengths{j} = cumsum(trajlengths{j});
            %             end
            try
                %             pfctrajdataconcat = cellfun(@(x) cell2mat(x'),pfctrajdata,'UniformOutput',false); %concatenate the trajectories
                currtrajs = {PEcoefallv2{pfcPEmatch,7}}; %need to put into additional cell wrap in order to concatenate next
                pfctrajdataconcat = cellfun(@(x) cell2mat(x),currtrajs,'UniformOutput',false); %concatenate the trajectories
                pfctrajdataconcat = cell2mat(pfctrajdataconcat);
                pfctrajlengths = cellfun(@(x) length(x),PEcoefallv2{pfcPEmatch,7},'UniformOutput',false); %get lengths of trajs
                pfctrajlengths = cell2mat(pfctrajlengths)';
                pfctrajlengths = cumsum(pfctrajlengths);
                
                tmpfig = figure; hold on;
                for j = 1:length(pfcmapdata);
                    subplot(4,length(pfcmapdata),j);
                    imagesc(pfcmapdata{j}); set(gca,'xtick',[], 'xticklabel', [], 'ytick', [], 'yticklabel', []);
                    epmeans(j) = pfcmeanrates(j);
                    if j == 1;
                        if i <= length(PFCindsExc);
                            colr = [1 0 0];
                            title(sprintf('EXC(%d %d %d %d) uFR(%0.2f) PEC:(%0.2f, %0.2f)',PFCindsmod(i,:), nanmean(epmeans), PEcoefallv2{pfcPEmatch,3}'),'Color', colr, 'FontSize', 9);
                        else colr = [0 0 1];
                            title(sprintf('INH(%d %d %d %d) uFR(%0.2f) PEC:(%0.2f, %0.2f)',PFCindsmod(i,:), nanmean(epmeans), PEcoefallv2{pfcPEmatch,3}'),'Color', colr, 'FontSize', 9);
                        end
                    end
                    
                end
                
                %magenta title if exc, blue if inh
                %                 %Ripple position
                %                 subplot(3,length(pfcmapdata),j+length(pfcmapdata));set(gca,'xtick',[], 'xticklabel', [], 'ytick', [], 'yticklabel', []);
                %                 scatter(pfcRipPosdata{j}(:,1), pfcRipPosdata{j}(:,2), 10, pfcRipPosSpikedata{j},'fill'); %pfcRipPosSpikedata acts as color LUT... scales based on # of spikes %'MarkerEdgeColor','b','MarkerFaceColor','k','LineWidth',1.5 );
                %                 axis([min(pfcRipPosdata{j}(:,1)) max(pfcRipPosdata{j}(:,1)) min(pfcRipPosdata{j}(:,2)) max(pfcRipPosdata{j}(:,2))]);
                %                 set(gca,'xtick',[], 'xticklabel', [], 'ytick', [], 'yticklabel', []); % title('#spikes/ripple');
                %                 colorbar('Location', 'eastoutside');
                
                pfctrajNORM = cellfun(@(x) cell2mat(x'),PEcoefallv2{pfcPEmatch,6},'UniformOutput',false); %concatenate the normalized trajectories
                pfctrajNORM = cell2mat(pfctrajNORM);
                subplot(4,length(pfcmapdata),1+length(pfcmapdata):length(pfcmapdata)*2);
                plot(pfctrajNORM(1,:),'b'); hold on;
                plot(pfctrajNORM(4,:),'r');
                area(PEcoefallv2{pfcPEmatch,5}(1,:),'FaceColor',[.5 .5 .5]); alpha(0.5);
                
                subplot(4,length(pfcmapdata),1+length(pfcmapdata)*2:length(pfcmapdata)*3);
                plot(pfctrajNORM(2,:),'g'); hold on;
                plot(pfctrajNORM(3,:),'c');
                area(PEcoefallv2{pfcPEmatch,5}(2,:),'FaceColor',[.5 .5 .5]); alpha(0.5);
                
                
                %                 pfctrajconcatNORM = cellfun(@(x) cell2mat(x'),PEcoefallv2{pfcPEmatch,6},'UniformOutput',false); %concatenate the normalized trajectories
                %                 pfctrajconcatNORM = cell2mat(pfctrajconcatNORM);
                %                 mintrajsconcat = [PEcoefallv2{pfcPEmatch,5}(1,:) PEcoefallv2{pfcPEmatch,5}(2,:) PEcoefallv2{pfcPEmatch,5}(2,:) PEcoefallv2{pfcPEmatch,5}(1,:)];
                %
                %                 plot(pfctrajconcatNORM,'k'); alpha(0.5); hold on;
                %                 plot(mintrajsconcat,'c'); alpha(0.5);
                
                subplot(4,length(pfcmapdata),1+length(pfcmapdata)*3:length(pfcmapdata)*4);
                plot(pfctrajdataconcat,'k'); hold on;
                plot(1:pfctrajlengths(1),pfctrajdataconcat(1:pfctrajlengths(1)),'b'); hold on;
                plot(pfctrajlengths(1):pfctrajlengths(2),pfctrajdataconcat(pfctrajlengths(1):pfctrajlengths(2)),'g');
                plot(pfctrajlengths(2):pfctrajlengths(3),pfctrajdataconcat(pfctrajlengths(2):pfctrajlengths(3)),'c');
                plot(pfctrajlengths(3):pfctrajlengths(4),pfctrajdataconcat(pfctrajlengths(3):pfctrajlengths(4)),'r');
                axis([0 length(pfctrajdataconcat) 0 max(pfctrajdataconcat)]);
                %                 arrayfun(@(x) line('XData',[x x], 'YData', [0 max(pfctrajdataconcat{j}(:,5))], 'LineStyle', '-', 'LineWidth', 2, 'Color',[.6 .6 .6]), trajlengths{j}, 'UniformOutput', false); % add traj edge lines to plot
                arrayfun(@(x) line('XData',[x x], 'YData', [0 max(pfctrajdataconcat)], 'LineStyle', '-', 'LineWidth', 2, 'Color',[.6 .6 .6]), pfctrajlengths, 'UniformOutput', false); % add traj edge lines to plot
                set(gca,'xtick',[], 'xticklabel', []); %title('C-L  L-C  C-R  R-C');
                
                
                %
                % %ISI data
                %  subplot(4,length(pfcISIdata),j+length(pfcISIdata)*3);
                %                 bins = linspace(min(pfcISIdata{j}),max(pfcISIdata{j}),length(pfcISIdata{j}));
                %                 n = hist(pfcISIdata{j},length(pfcISIdata{j}));
                %                 lower95 = length(n(cumsum(n)<(length(pfcISIdata{j})*0.90))); %only use the lower 95% of the isi bins.. i.e. discard the long tail that makes it hard to visualize
                %                 bar(bins(xlabel('ISI(s)'); ylabel('count');1:lower95), n(1:lower95),'EdgeColor','none');
                %                 axis([0 bins(lower95) 0 max(n)]);
                %                 xlabel('ISI(s)'); ylabel('count'); title(sprintf('mean %0.1f Hz',pfcmeanrates(j)));
                
                %         end
                
                
                
                %             if i <= length(PFCindsExc);
                %                 colr = [1 0 0];
                %                 supertitle({sprintf('SWR Excited (%d %d %d %d)', PFplot(1:pfctrajlengths(1),pfctrajdataconcat(1:pfctrajlengths(1));Cindsmod(i,:))}, colr); % I added a color var to supertitle as the second input
                %             else colr = [0 0 1];
                %                 supertitle({sprintf('SWR Inhibited (%d %d %d %d)', PFCindsmod(i,:))}, colr); % I added a color var to supertitle as the second input
                %             end %magenta title if exc, blue if inh
                
                if pausefigs; keyboard; end %pause
                if savefigs
                    %                 figfile = [figdir,savefigfilename,'/',sprintf('An%dD%dT%dC%d', PFCindsmod(i,1), PFCindsmod(i, 2), PFCindsmod(i, 3), PFCindsmod(i, 4))];
                    figfile = [figdir,savefigfilename,'/',sprintf('pfcripmodmapisi%d', i)];
                    print('-dpng', figfile, '-r300');
                end
                close(tmpfig)
            catch
                not4trajs = [not4trajs; PFCindsmod(i,:)];
            end
        end
    end
end








