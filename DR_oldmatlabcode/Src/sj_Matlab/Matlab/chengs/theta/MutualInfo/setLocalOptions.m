clear all
global fmaux % file manager auxillary variable

% set search paths
addpath('/bach/theta/common')
addpath('/bach/theta/code')
addpath('/bach/AdaptFilter/main')

fmaux.prefix='';

fmaux.selectid='100s';
fmaux.select=['./select-' fmaux.selectid];

% directories
fmaux.datadir= './';
fmaux.EEGdir= './';
fmaux.data2dir= './';

fmaux.loaded=[];

