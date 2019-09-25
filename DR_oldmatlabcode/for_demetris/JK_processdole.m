%% script to process dol data
clear all; close all; clc;


% day process
%sj_dayprocess_cont2(
currentdir = pwd;
cd /opt/data40/jeff/Dole
directoryname = '/opt/data40/jeff/Dole/Dol/';
fileprefix = 'dol';

days = 15;
tet_rip = [1:5 7:12]';
tet_theta = [1:5 7:12]';
for kday = days
    if kday<10
        daydir = ['Dole0' num2str(kday)];
    else
        daydir = ['Dole' num2str(kday)];
    end
    makedayparms(daydir,30,200,'pos',0)
    jk_dayprocess(daydir,'/opt/data40/jeff/Dole/Dol', 'dol',kday, 'cmperpix', 43/154)
    sj_dayprocess_cont2(daydir,'/opt/data40/jeff/Dole/Dol', 'dol',kday, 'cmperpix', 43/154)
    tetdaylist = [kday*ones(size(tet_rip)) tet_rip];
    jk_rippledayprocess(directoryname,fileprefix,kday,'daytetlist',tetdaylist)
    tetdaylist = [kday*ones(size(tet_theta)) tet_theta];
    jk_thetadayprocess(directoryname, fileprefix,kday,'daytetlist',tetdaylist)
end

%%
tetrodes = 1:5; % specify all tetrodes
% extract ripples
mindur = 0.015; % 15 ms minimum duration
nstd = 7; % 2 std
for d = days
    JK_extractripples(directoryname, fileprefix, d, tetrodes, mindur, nstd);
end
tetrodes = [1:5 7:12]; % specify all tetrodes
% extract theta
mindur = 1; % 1 s minimum duration
nstd = 1; % 1 std
for d = days
    JK_extracthightheta(directoryname, fileprefix, d, tetrodes, mindur, nstd);
end


cd(currentdir)