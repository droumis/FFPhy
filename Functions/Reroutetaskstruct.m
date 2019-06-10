function Reroutetaskstruct(directoryname,day, varargin)
%
%
%This program creates the 'task' structure for each day's data. This
%structure contains information about the exposure number for each epoch
%and the linearization coordinates for each run session. If a task file
%already exists for one of the days in index, then that file is
%loaded and written to to prevent unwanted deletion of data.
%
%      --INPUTS--
%
%    directoryname - example 'data99/user/animaldatafolder/', a folder
%                    containing processed matlab data for the animal
%
%
%     coordprogram - a string with the name of the function that will
%                    provide the track segment location info. This
%                    program should take three inputs - 1) directoryname
%                    2) pos.data for an epoch, and 3) the epoch index [day epoch]
%                    The output of the program is a cell
%                    array, where each cell describes one trajectory on
%                    the track (do not count foreward and backward
%                    movement as 2 trajectories).  Inside each cell is
%                    a three dimentional M-by-2-by-N matrix.  The
%                    matrix gives the x and y coodinates for the M
%                    segment endpoints, across N timeframes.  N should
%                    equal the length of pos.data.  The
%                    segment end points should start with the position
%                    that the user wants to be 0 on the linear scale and
%                    progress to the last segment.  All endpoints that
%                    are shared by two segments should only exist
%                    once per trajectory.  The function is called like
%                    this: coords = coordprogram(direcoryname,pos{day}{epoch}.data,[day epoch])
%    ---OPTIONS----
%
%        overwrite - 1 to overwrite old task file (default), 0 to keep the same file
%                    and just change the epochs defined in index
%   lowercasethree - variable prefix. Default is ''.
%    combineepochs - value of 1 combines all the epochs for the day when
%                    calling coordprogram. Default 0;
%
%


% get the current directory to go back to after the function is done
currentdir = pwd;

% define coordprogram specific reroute program
coordprogram='getcoord_reroute_arms'; % variant of the maze with arms from each corner
%coordprogram='getcoord_reroute'; % original reroute maze

% set default variables
lowercasethree = '';
overwrite = 1;
combineepochs = 0;

%set variable options
for option = 1:2:length(varargin)-1
    
    switch varargin{option}
        case 'lowercasethree'
            lowercasethree = varargin{option+1};
        case 'overwrite'
            overwrite = varargin{option+1};
        case 'combineepochs'
            combineepochs = varargin{option+1};
        otherwise
            error(['Option ''', varargin{option}, ''' not defined']);
    end
end



eval(['coordProgramHandle = @',coordprogram,';']);
task = [];
currday = 0;

% ---- File loading ----
% See if day number needs 0
dsz = '';
    if (day < 10)
        dsz = '0';
    end

% Specify data directory and load the file
%animalname = animaldir(1:end-1);
datadir = '/data14/jai/';
animalname = directoryname(1:end-1);
% Converts the day and epoch into text
dayt = num2str(day);

sfilename = strcat(datadir,directoryname,'/',animalname,'data',dsz,dayt,'.mat');

% Load the mat file

load(sfilename);

% define pos file

% go through each run epoch of the day

% load the times file, either automatically generated by NSpike process or
% manually generated with maketimesfile.m, which helps if there are breaks
% in the data

for k=1:size(data{1,day},2);
    
    if ~isempty(data{1,day}{1,k})
       
    if strcmp(data{1,day}{1,k}.Stats.epochtype,'Run')==1;
        
        posdata = data{1,day}{1,k}.Pos.correcteddata;
        tmplinearcoord = feval(coordProgramHandle,directoryname,posdata); % remember pixel dimension is hard coded
        
        S = size(posdata,1);
        task{day}{k}.type = 'reroute';
        %all time points are the same, so just grab the number of
        %points for the epoch
        for reshape = 1:length(tmplinearcoord)
            epochlinearcoord{reshape} = tmplinearcoord{reshape}(:,:,1:S);
        end
        task{day}{k}.linearcoord = epochlinearcoord;
        
        %if ~isempty(coordprogram)
        %        %run the user-defined coordinate program
        %        task{day}{k}.linearcoord = feval(coordProgramHandle,directoryname,posdata,[day k]);
        %end
    end
    
    k=k+1;
    else k=k+1;
    
    
end
%save the task file
dsz = '';
if (day < 10)
    dsz = '0';
end


outdir = strcat('/data14/jai/',directoryname,'/');
cd(outdir);
filename = sprintf('%stask%s%s.mat',animalname,dsz,num2str(day));
save(filename,'task');
cd(pwd);

end

