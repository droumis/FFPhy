%% script to process bob data
clear; close all;


% day process
%sj_dayprocess_cont2(
currentdir = pwd;
rawdatafile = '/data19/droumis/bob/';
cd (rawdatafile)
% directoryname = '/data19/droumis/bob/bob_proc';
animdirect = '/data19/droumis/bob_proc/'
fileprefix = 'bob';

% days = 25;
days = [13 15 16];
% tet_rip = [1 2 3 6 7 8 9 12]';
% tet_theta = [1 2 3 6 7 8 9 12]';
% cmperpix = [.16 .48 .16 .48 .16 .48 .16]  %alternating cam between sleep (.16cm/px) and run (.48cm/px).
thresh = 40;
maxallowedamp = 200;
for kday = days
    cd (rawdatafile)
    if kday<10
        daydir = ['bob0' num2str(kday)];
    else
        daydir = ['bob' num2str(kday)];
    end
%     makedayparms(daydir,40,200,'pos',0);
%     makedayparms(daydir,40,200)
    sj_makedayparms_pc2(daydir ,thresh, maxallowedamp);
%     sj_makedayparms_pc2(dayfolder ,thresh, maxallowedamp, varargin)
%     DR_kk_dayprocess(daydir,animdirect,'bob',kday,'cmperpix', cmperpix)
%     %run this after clustering
%     
%     
%     
%     jk_dayprocess(daydir,'/data19/droumis/bob/bobstim_proc', 'bob',kday, 'cmperpix', 43/154)
%     sj_dayprocess_cont2(daydir,'/data19/droumis/bob/bobstim_proc', 'bob',kday, 'cmperpix', 43/154)
%     tetdaylist = [kday*ones(size(tet_rip)) tet_rip];
%     jk_rippledayprocess(directoryname,fileprefix,kday,'daytetlist',tetdaylist)
%     tetdaylist = [kday*ones(size(tet_theta)) tet_theta];
%     jk_thetadayprocess(directoryname, fileprefix,kday,'daytetlist',tetdaylist)
end

%% holding off on ripples until after clustered
% tetrodes = [8 9]; % specify all tetrodes
% % extract ripples
% mindur = 0.015; % 15 ms minimum duration
% nstd = 5; % 2 std
% 
% for d = days
%     JK_extractripples(directoryname, fileprefix, d, tetrodes, mindur, nstd);
% end
% 
% 
% %extractripples(directoryname, fileprefix, day, tetrode, min_suprathresh_duration, nstd, varargin)
% %option: %            'excludenoise' - set the std threshold in vargargin +1 and
% %            minimum noise duration in varargin +2. if you want to
% %            specificy sample rate or start ind of eeg go to 
% %            JY_findhighamplitudeeeg function in this script and insert
% %            based on JY_findhighamplitudeeeg( eeg,samprate,eegstart,
% %            threshold, duration ). DR added 6/16/13
% 
% %% in progress
% % for d = days
% %     DR_extractripples_excludenoise(directoryname, fileprefix, d, tetrodes, mindur, nstd, 'excludenoise', .1, 0.2);
% % end
% % tetrodes = tet_rip; % specify all tetrodes
% 
% % 
% % % % extract theta
% % % mindur = 1; % 1 s minimum duration
% % % nstd = 1; % 1 std
% % % for d = days
% % %     JK_extracthightheta(directoryname, fileprefix, d, tetrodes, mindur, nstd);
% % % end
% % 
% 
% cd(currentdir)