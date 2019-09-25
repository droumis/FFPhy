



%%%%%%%%%%%%%%%%%%%
%Better Way - Generate ranges from names


names{1} = '1>> All points 00:05:52-02:08:45';
names{2} = '2  epoch1  00:05:52-00:25:52';
names{3} = '3  epoch2  00:27:05-00:47:05';
names{4} = '4  epoch3  00:48:10-01:08:10';
names{5} = '5  epoch4  01:09:10-01:30:25';
names{6} = '6  epoch5  01:31:20-01:51:35';
names{7} = '7  epoch6  01:53:45-01:51:35';
names{8} = '8  epoch7  01:42:48-01:52:48';


% names{1} = '1>> All points 00:05:52-02:08:45';
% names{2} = '2  Sleep1  00:05:52-00:25:52';
% names{3} = '3  Run1  00:27:05-00:47:05';
% names{4} = '4  Sleep2  00:48:10-01:08:10';
% names{5} = '5  Run2  01:09:10-01:30:25';
% names{6} = '6  Sleep3  01:31:20-01:51:35';
% names{7} = '7  Run3  01:53:45-01:51:35';
% names{8} = '8  Sleep4  01:42:48-01:52:48';



for n=1:length(names)
    
    currnames=names{n};
    [T,R]=strtok(currnames,'-');
    ranges(n,:)=timetrans({T(end-7:end) R(2:end)}, 10000, 2);

end

% % Display 
names{1}
names{2}
names{3}

ranges

save times names ranges
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