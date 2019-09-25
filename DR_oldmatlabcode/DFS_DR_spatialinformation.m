

%to do..
%rerun linfields without using shantanu's premade structures.
%currently i have to run the filter everyday bc it's using the date function to load the name.. change this?
%include rip positions in 2d maps

%this now also plots pfc rip exc and inh maps


warning('off','all');
% clear; close all;

runfiltflds = 1;
runscriptspatinfo = 0;

plotspatinfo = 0;

plotanything = 0;
savefigs = 0;
pausefigs = 0; %pause at each fig, keyboard press to continue

savedir = '/data19/sjadhav/HPExpt/ProcessedDataDR/';
% savefilename = sprintf('Flds_%s',date); %has to match saved data name. specify GLM or Corr pairs and filters used: velocity filter <=.....nrip >= (#tetrodes ripples detected); std > of ripple detection power
% savefilename = 'Flds_jan26';
savefilename = 'spatinfo_feb9';
% savefigfilename = sprintf('SpatRipcorr_%s', date); %fig name
savefigfilename = sprintf('RipModspatinfo_%s', date); %fig name
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
% spks = setfilterfunction(spatf, 'DFAdr_spikes',{'spikes', 'linpos', 'pos'}, 'binsize', 1, 'std', 2); %spikes for ISI
% fields = setfilterfunction(spatf, 'DFA_loadfields_DR', {'linfields', 'mapfields', 'cellinfo'}); %map and lin fields
spatinf = setfilterfunction(spatf, 'getspatialinfo_DR2', {'linpos', 'spikes'}, 6,  minV, 'peakrate', 3,'appendindex', 1, 'incltraj', [1:4]); %I changed it to output all trajs' summed per cell, not combined across eps yet
% modf = setfilterfunction(spatf,'DFAsj_getripalignspiking_position',{'spikes', 'ripples', 'tetinfo', 'pos'}); %ripple position

% Run analysis-----------------------
%     flds = runfilter(psf);  % Place Field Stability.. trajectories
% spks = runfilter(spks);  %spatf = setfilterfunction(spatf, 'getspatialinfo_DR', {'linpos', 'spikes'},6,  minV, 'peakrate', 3,'appendindex', 1, 'incltraj', traj); %I changed it to output all trajs' summed per cell, not combined across eps yet spike data in order to get ISI
% flds = runfilter(fields);  % Place Field Map
spatinf = runfilter(spatinf); %spatial information
% ripP = runfilter(modf); %ripple position

%     end
disp('Finished running filter script');
%--------------------- Finished Filter Function Run -------------------
%clear some variables so that when the mat gets loaded in, it doesn't override the intended variable vals.
clear  runfiltflds runripposscript runscript runPECv2 runPEC runscriptmaps plotPFCripmodISI plotPFCripmodMapsTrajsISImeanrate cyclefigs savefigs  plotPFCCA1MapsTrajs plotPFCCA1Trajs plotCorrCoef loadrippos  loadpairindices savefilename savedir  savefigfilename  savedata  figdir savefigs pausefigs plotanything plotPFCripmodMapsTrajs plotPFCripmodMeanRates
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
