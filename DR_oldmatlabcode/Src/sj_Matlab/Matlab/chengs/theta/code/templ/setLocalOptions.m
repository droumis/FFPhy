clear all
global fmaux % file manager auxillary variable

% set search paths
addpath('/bach/theta/common')
addpath('/bach/theta/code')
addpath('/bach/AdaptFilter2/main')

fmaux.prefix='ter';

%fmaux.cellselect='/bach/Ter/data2/cellselect-Adaptx3-0t0-05d';
%fmaux.cellselect='/bach/Ter/data2/cellselect-all';
fmaux.cellselect='/bach/Ter/data2/cellselect-test';
%fmaux.cellselect='/bach/Ter/data2/cellselect-one';

%fmaux.tEEG= 14;
fmaux.tEEG= 2;

% directories
fmaux.datadir= '/bach/Ter/data/';
fmaux.EEGdir= '/bach/Ter/data/EEG/';
fmaux.data2dir= '/bach/Ter/data2/';
%fmaux.resultsdir= '/bach/Ter/results/'
%fmaux.adaptdir= '/bach/Ter/adapt/';

fmaux.loaded=[];

