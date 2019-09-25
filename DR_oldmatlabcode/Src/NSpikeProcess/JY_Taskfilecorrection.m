

function [Epochdata]=JY_Taskfilecorrection(animaldir,animalname,daylist,scalefactor)

% load existing task file and divide each cell value by a predefined value
% use when the wrong scaling factor was initially used to specify
% intersections on the track


for ii=daylist
    
    day=ii;
    
    % ---- File loading ----
    % See if day number needs 0
    dsz = '';
    if (day < 10)
        dsz = '0';
    end
    
    % Specify data directory and load the file
    %animalname = animaldir(1:end-1);
    datadir = '/data14/jai/';
    % Converts the day and epoch into text
    dayt = num2str(day);
    
    sfilename = strcat(datadir,animaldir,'/',animalname,'task',dsz,dayt,'.mat');
    
    % Load the mat file
    
    load(sfilename);
    %loop for each epoch
    for iii=1:size(task{1,day},2)
        
        if ~isempty(task{1,day}{1,iii})
            % check if an original version of the data exists
            
            if ~isfield(task{1,day}{1,iii},'linearcoord_original')
                
                task{1,day}{1,iii}.linearcoord_original=task{1,day}{1,iii}.linearcoord;
                task{1,day}{1,iii}.linearcoord=cellfun(@(x) x*scalefactor, task{1,day}{1,iii}.linearcoord, 'UniformOutput' ,false );
                
            else
                
                % only modify the original data
                task{1,day}{1,iii}.linearcoord=cellfun(@(x) x*scalefactor, task{1,day}{1,iii}.linearcoord_original, 'UniformOutput' ,false );
            end
        end
    end
    
    
    outdir = strcat(datadir,animaldir,'/');
    cd(outdir);
    filename = sprintf('%stask%s%s.mat',animalname,dsz,num2str(day));
    save(filename,'task');
end
end

