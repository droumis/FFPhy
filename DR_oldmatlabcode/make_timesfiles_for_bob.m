% clear all;
% close all;
% clc;


cd /data19/droumis/bob


%% bob18

dirname = 'bob25';
names{1} = '1  All points 00:46:00-02:29:10';
names{2} = '2  run1  00:46:00-00:51:00';
names{3} = '3  sleep1  00:51:00-00:56:00';
names{4} = '4  run2 00:56:00-01:01:00';
names{5} = '5  sleep2  01:01:00-01:06:00';
names{6} = '6  run3  01:06:00-01:11:00';
names{7} = '7  sleep3  01:11:00-01:16:00';
names{8} = '8  run4  01:16:00-01:21:00';
names{9} = '9 sleep4 01:21:00-01:26:00';
names{10} = '10 run5  01:27:00-01:32:00';
names{11} = '11 sleep5  01:32:00-01:37:00';
names{12} = '12 run6 01:38:00-01:43:00';
names{13} = '13  sleep6 01:43:00-01:48:00';
names{14} = '14  run6  01:49:00-01:54:00';
names{15} = '15  sleep7 01:54:00-01:59:00';
names{16} = '16  run7  01:59:00-02:04:00';
names{17} = '17  sleep8  02:04:00-02:09:00';
names{18} = '18  run8  02:09:00-02:14:00';
names{19} = '19  sleep9  02:14:00-02:19:00';
names{20} = '20 run9 02:19:00-02:24:00';
names{21} = '21 sleep10  02:24:10-02:29:10';

for n=1:length(names)
    currnames=names{n};
    [T,R]=strtok(currnames,'-');
    ranges(n,:)=timetrans({T(end-7:end) R(2:end)}, 10000, 2);

end
% % Display 
save(fullfile(dirname, 'times.mat'),'names','ranges')
disp(['Saved Times File for ' dirname]);

%%
% 
% Dole08
% dirname = 'Dole08';
% names{1} = '1  All points 03:20:50-03:35:50';
% names{2} = '2  Sleep1  03:20:50-03:35:50';
% names{3} = '3  Run1  03:38:30-03:53:30';
% names{4} = '4  Sleep2  03:55:00-04:10:00';
% names{5} = '5  Run2  04:12:10-04:27:10';
% names{6} = '6  Sleep3  04:28:20-04:43:20';
% names{7} = '7  Run3  04:44:10-04:59:10';
% names{8} = '8  Sleep4  05:01:45-05:16:45';
% for n=1:length(names)
%     currnames=names{n};
%     [T,R]=strtok(currnames,'-');
%     ranges(n,:)=timetrans({T(end-7:end) R(2:end)}, 10000, 2);
% 
% end
% Display 
% save(fullfile(dirname, 'times.mat'),'names','ranges')
% disp(['Saved Times File for ' dirname]);


% %% Dole12
% clear all;
% dirname = 'Dole12';
% names{1} = '1  All points 00:01:50-00:56:08';
% names{2} = '2  Sleep1  00:41:08-00:56:08';
% names{3} = '3  Run1  00:56:24-01:11:24';
% % names{4} = '4  Sleep2  01:11:54-01:26:54';
% 
% for n=1:length(names)
%     currnames=names{n};
%     [T,R]=strtok(currnames,'-');
%     ranges(n,:)=timetrans({T(end-7:end) R(2:end)}, 10000, 2);
% 
% end
% % % Display 
% save(fullfile(dirname, 'times.mat'),'names','ranges')
% disp(['Saved Times File for ' dirname]);
% 
% 
% 

% %% Dole13
% clear all;
% dirname = 'Dole13';
% names{1} = '1  All points 00:01:50-00:16:50';
% names{2} = '2  Sleep1  00:01:50-00:16:50';
% 
% for n=1:length(names)
%     currnames=names{n};
%     [T,R]=strtok(currnames,'-');
%     ranges(n,:)=timetrans({T(end-7:end) R(2:end)}, 10000, 2);
% 
% end
% % % Display 
% save(fullfile(dirname, 'times.mat'),'names','ranges')
% disp(['Saved Times File for ' dirname]);
% 
% 
% %% Dole15
% clear all;
% dirname = 'Dole15';
% names{1} = '1  All points 00:27:04-00:42:04';
% names{2} = '2  Sleep1  00:27:04-00:42:04';
% 
% for n=1:length(names)
%     currnames=names{n};
%     [T,R]=strtok(currnames,'-');
%     ranges(n,:)=timetrans({T(end-7:end) R(2:end)}, 10000, 2);
% 
% end
% % % Display 
% save(fullfile(dirname, 'times.mat'),'names','ranges')
% disp(['Saved Times File for ' dirname]);

