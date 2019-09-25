



%%%%%%%%%%%%%%%%%%%
%Better Way - Generate ranges from names

names{1} = '1>> All points 00:13:51-01:43:13';
names{2} = '2  Sleep1 00:13:51-00:29:12';
names{3} = '3  Run1  00:34:10-00:49:24';
names{4} = '4  Sleep2  00:51:07-01:06:21';
names{5} = '5  Run2 01:10:25-01:26:25';
names{6} = '6  Sleep3 01:28:01-01:43:13';

%

for n=1:length(names)
%for n=2:6       
    currnames=names{n};
    [T,R]=strtok(currnames,'-');
    ranges(n,:)=timetrans({T(end-7:end) R(2:end)}, 10000, 2);
end

% Display 
names{1}
names{2}
names{3}
names{4}
names{5}
names{6}
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