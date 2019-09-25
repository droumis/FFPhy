clear all; clc; close all; 

daydir = '/data19/droumis/bob/bob16';
cd (daydir)
pwd

%%%%%%%%%%%%%%%%%%%
%Better Way - Generate ranges from names

% % bob13
% names{1} = '1>> All points 00:37:15-02:42:05';
% names{2} = '2  Sleep1  00:37:15-00:52:15';
% names{3} = '3  Run1  00:53:10-01:08:10';
% names{4} = '4  Sleep2  01:20:35-01:35:35';
% names{5} = '5  Run2  01:36:10-01:51:10';
% names{6} = '6  Sleep3  01:54:50-02:09:50';
% names{7} = '7  Run3  02:10:20-02:25:20';
% names{8} = '8  Sleep4  02:27:05-02:42:05';


%bob15
% names{1} = '1>> All points 0:17:15-2:09:00';
% names{2} = '2  Sleep1  00:17:15-00:32:15';
% names{3} = '3  Run1  00:33:50-00:48:50';
% names{4} = '4  Sleep2  00:49:20-01:04:20';
% names{5} = '5  Run2  01:06:15-01:21:15';
% names{6} = '6  Sleep3  01:22:00-01:37:00';
% names{7} = '7  Run3  01:38:00-01:53:00';
% names{8} = '8  Sleep4  01:54:00-02:09:00';

%bob16
names{1} = '1>> All points 00:20:35-02:37:10';
names{2} = '2  Sleep1  00:20:35-00:35:35';
names{3} = '3  Run1  00:59:20-01:14:20';
names{4} = '4  Sleep2  01:17:35-01:32:35';
names{5} = '5  Run2  01:33:15-01:48:15';
names{6} = '6  Sleep3  01:49:00-02:04:00';
names{7} = '7  Run3  02:04:45-02:19:45';
names{8} = '8  Sleep4  02:22:10-02:37:10';


for n=1:length(names);
    
    currnames=names{n};
    [T,R]=strtok(currnames,'-');
    ranges(n,:)=timetrans({T(end-7:end) R(2:end)}, 10000, 2);

end

% % Display 
names{1};
names{2};
names{3};

ranges

save times names ranges;
disp('Saved Times File');


%%%% Enter ranges and times manually

% ranges(1,:)=timetrans({'00:06:49' '01:30:37'}, 10000, 2);
% names{1} = '1>> All points';  
% 
% ranges(2,:)=timetrans({'00:06:49' '00:21:55'}, 10000, 2);
% names{2} = '2  Sleep1 00:06:49-00:21:55';
% 
% ranges(3,:)=timetrans({'00:25:25' '00:40:26'}, 10000, 2);
% names{3} = '3  Run1 00:25:25-00:40:26';
% 
% ranges(4,:)=timetrans({'00:41:35' '00:56:38'}, 10000, 2);
% names{4} = '4  Sleep2 00:41:35-00:56:38';
% 
% ranges(5,:)=timetrans({'00:58:59' '01:14:00'}, 10000, 2);
% names{5} = '5  Run2 00:58:59-01:14:00';
% 
% ranges(6,:)=timetrans({'01:15:08' '01:30:37'}, 10000, 2);
% names{6} = '6  Sleep3 01:15:08-01:30:37';
% 
% names;
% ranges;
% 
% save times;