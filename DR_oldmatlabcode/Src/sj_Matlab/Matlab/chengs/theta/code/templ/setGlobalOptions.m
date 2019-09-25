clear all
global fmaux % file manager auxillary variable

% set search paths
addpath('/bach/theta/common')
addpath('/bach/theta/code')
addpath('/bach/AdaptFilter2/main')

fmaux.prefix='ter';
fmaux.cellselect='/bach/Ter/data2/cellselect-all';
fmaux.tEEG= 2;

% directories
fmaux.datadir= '/bach/Ter/data/';
fmaux.EEGdir= '/bach/Ter/data/EEG/';
fmaux.data2dir= '/bach/Ter/data2/';
fmaux.loaded=[];
