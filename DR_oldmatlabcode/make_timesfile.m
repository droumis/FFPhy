



%%%%%%%%%%%%%%%%%%%
%Better Way - Generate ranges from names


names{1} = '1>> All points 01:05:15-03:27:50';
names{2} = '2  Sleep1  01:05:15-01:36:10';
names{3} = '3  Run1  01:37:47-01:55:15';
names{4} = '4  Sleep2  01:56:00-02:18:09';
names{5} = '5  Run2  02:19:09-02:34:59';
names{6} = '6  Sleep3  02:36:18-03:01:30';
names{7} = '7  Run3  03:04:53-03:17:39';
names{8} = '8  Sleep4  03:18:20-03:27:50';



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