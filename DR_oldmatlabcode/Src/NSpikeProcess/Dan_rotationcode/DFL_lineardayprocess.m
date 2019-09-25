function DFL_lineardayprocess(directoryname,fileprefix,days, varargin)
%LINEARDAYPROCESS(directoryname,fileprefix,days, options)
%
%Runs linearizeposition for all run epochs in each day and saves the data in
%'linpos' in the directoryname folder.  See LINEARIZEPOSITION for the definitions
%of the options.
%
%directoryname - example 'data99/user/animaldatafolder/', a folder
%                containing processed matlab data for the animal
%
%fileprefix -    animal specific prefix for each datafile
%
%days -          a vector of experiment day numbers
%

lowercasethree = '';

%set variable options
for option = 1:2:length(varargin)-1
    switch varargin{option}
        case 'lowercasethree'
            lowercasethree = varargin{option+1};
            
    end
end


days = days(:)';

for day = days
    dsz = '';
    if (day < 10)
        dsz = '0';
    end
    eval(['load ',directoryname,fileprefix,'task',dsz,num2str(day), '.mat']);
    eval(['task = ',lowercasethree,'task;'])
    for i =  1:length(task{day})
        if ((~isempty(task{day}{i})) && (strcmp(task{day}{i}.type,'reroute')) )
            
            disp(['Day ',num2str(day), ', Epoch ',num2str(i)])
            
            index = [day i];
           
            % use modified version of lineardayprocess to assign each position to a track segment
            % 
            [linpos_temp{day}{i}.statematrix,linpos_temp{day}{i}.segmenttable, linpos_temp{day}{i}.trajwells,...
                linpos_temp{day}{i}.wellSegmentInfo, linpos_temp{day}{i}.segmentInfo, linpos_temp{day}{i}.trajmatrix,...
                linpos_temp{day}{i}.trialsegments, linpos_temp{day}{i}.trialsegments_raw] = DFL_linearizereroutepositionTemp(directoryname,...
                fileprefix, index, varargin);
            linpos_temp{day}{i}.statematrixfields = ['time startwell endwell segment segheaddir velocity lineardist'];
            linpos_temp{day}{i}.segmenttablefields = ['trajnum segnum segmentID'];
            
        end
    end
end
%group all trajectory trials together into one big matrix
trialsegments_all = {};
for day = days
    dsz = '';
    if (day < 10)
        dsz = '0';
    end
    eval(['load ',directoryname,fileprefix,'task',dsz,num2str(day), '.mat']);
    eval(['task = ',lowercasethree,'task;']);
    
    for i = 1:length(task{day})
        if ~isempty(task{day}{i}) && strcmp(task{day}{i}.type, 'reroute')
            trialsegments_all = [trialsegments_all linpos_temp{day}{i}.trialsegments];
        end
    end
end

%calculate trajectories
if ~isempty(trialsegments_all)
    [labeled_trial_segments segmenttable trajwells] = DFL_label_trial_segments(trialsegments_all);
else
    disp('Warning: No Trial Segments Found for any day');
end

%seperate trajectories back into each day, epoch
ind = 1;
for day = days
    
    linpos = [];
    dsz = '';
    if (day < 10)
        dsz = '0';
    end
    
    eval(['load ',directoryname,fileprefix,'task',dsz,num2str(day), '.mat']);
    eval(['task = ',lowercasethree,'task;']);
    
    for i = 1:length(task{day})
        if ~isempty(task{day}{i}) && strcmp(task{day}{i}.type, 'reroute') && ~isempty(trialsegments_all)
            trial_len = length(linpos_temp{day}{i}.trialsegments);
            if trial_len > 0
                linpos{day}{i} = linpos_temp{day}{i};
                linpos{day}{i}.trialsegments = labeled_trial_segments(:,ind:ind+trial_len-1);
                linpos{day}{i}.segmenttable = segmenttable;
                linpos{day}{i}.trajwells = trajwells;
                
                ind = ind + trial_len;
            end
        end
    end
    eval([lowercasethree,'linpos = linpos;']);
    filename=strcat(directoryname,fileprefix,'linpos',dsz,num2str(day),lowercasethree);
    save(filename,'linpos');
    
end
