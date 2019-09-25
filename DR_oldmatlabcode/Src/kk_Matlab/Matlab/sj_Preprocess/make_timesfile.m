



%%%%%%%%%%%%%%%%%%%
%Better Way - Generate ranges from names

clear all

names{1} = '1>> All points 00:48:15-02:00:42';
names{2} = '2  Sleep1 00:48:15-00:53:00';
names{3} = '3  Sleep2 00:54:42-00:58:29';
names{4} = '4  Run1 01:06:42-01:13:12';
names{5} = '5  Run2 01:14:23-01:20:32';
names{6} = '6  Sleep3 01:21:31-01:31:12';
names{7} = '7  Sleep4 01:34:18-01:35:42';
names{8} = '8  Run3 01:37:17-01:40:12';
names{9} = '9  Run4 01:41:32-01:52:04';
names{10}= '10  Sleep5 01:53:22-02:00:42';


%

for n=1:length(names)
%for n=2:6       
    currnames=names{n};
    [T,R]=strtok(currnames,'-');
    ranges(n,:)=timetrans({T(end-7:end) R(2:end)}, 10000, 2);
end

% Display 

ranges
names

save times names ranges
disp('Saved Times File');

clear all

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