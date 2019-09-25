function createnovelobjecttaskstruct(directoryname,fileprefix,index, quadrants)
%
%createtaskstruct(directoryname,fileprefix,index, coordprogram,options)
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
%       fileprefix - animal specific prefix for each datafile
%    index - each row describes a run epoch - [day epoch]
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
%
%       quadrants - x y coordinates of the middle location. Quadrants is
%       used to define the quadrants in the following manner:
%
%                               x-------x-------x
%                               |       |       |
%                               |   1   |   2   |
%                               |       |       |
%                               x-------x-------x
%                               |       |       |
%                               |   4   |   3   |
%                               |       |       |
%                               x-------x-------x
%                           


%set default variables

%set variable options
% for option = 1:2:length(varargin)-1
%     
%     switch varargin{option}
%         case 'combineepochs'
%             combineepochs = varargin{option+1}; 
%         otherwise
%             error(['Option ''', varargin{option}, ''' not defined']);
%     end
% end

index = sortrows(index,1);

task = [];
currday = 0;

for i = 1:size(index,1)    
    day = index(i,1);
    epoch = index(i,2); 
    if (day ~= currday)
        if (currday ~= 0)
            %we have moved on to a new day, so save the previous day's data
            dsz = '';
            if (currday < 10)
                dsz = '0';
            end
            filename = [directoryname,fileprefix,'task',dsz,num2str(currday)];
            save(filename,'task');

            task = [];
        end
        currday = day;
        dsz = '';
        if (day < 10)
            dsz = '0';
        end
        %load the new day's data
        eval(['load ',directoryname,fileprefix,'pos', dsz, num2str(day), '.mat']);
    end
    
    %fill the structure with data
    task{day}{epoch}.type = 'run';
    
    %Figure out which quadrant the rat is in at each time
    quad = -1*ones(size(pos{day}{epoch}.data(:,2)));
    quad(pos{day}{epoch}.data(:,2) < quadrants(1) & pos{day}{epoch}.data(:,3) > quadrants(2)) = 1;
    quad(pos{day}{epoch}.data(:,2) > quadrants(1) & pos{day}{epoch}.data(:,3) > quadrants(2)) = 2;
    quad(pos{day}{epoch}.data(:,2) > quadrants(1) & pos{day}{epoch}.data(:,3) < quadrants(2)) = 3;
    quad(pos{day}{epoch}.data(:,2) < quadrants(1) & pos{day}{epoch}.data(:,3) < quadrants(2)) = 4;
    
    task{day}{epoch}.time = pos{day}{epoch}.data(:,1);
    task{day}{epoch}.quadrants = quad;
    task{day}{epoch}.linearcoords = quadrants;
    
    
end
%save the final day's task file
dsz = '';
if (currday < 10)
    dsz = '0';
end
filename = [directoryname,fileprefix,'task',dsz,num2str(currday)];
save(filename,'task');


