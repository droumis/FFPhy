function add_ripplephase_to_spikes(directoryname,fileprefix,days, varargin)

%
%Goes through each cell and adds the current ripple phase recorded on a 
%representative CA1 tetrode to the spikes structure
%
%directoryname - example '/data99/user/animaldatafolder/', a folder 
%                containing processed matlab data for the animal
%
%fileprefix -    animal specific prefix for each datafile (e.g. 'fre')
%
%days -          a vector of experiment day numbers 
%
%options -

%set variable options
for option = 1:2:length(varargin)-1
    switch varargin{option}
        case 'ref_tetrode'
            ref_tetrode = varargin{option+1};
    end
end


% check to see if the directory has a trailing '/'
if (directoryname(end) ~= '/')
    warning('directoryname should end with a ''/'', appending one and continuing');
    directoryname(end+1) = '/';
end

days = days(:)';
load '/usr/local/filtering/ripplefilter.mat'

for day = days
    %load up the spike file
    spikes = loaddatastruct(directoryname,fileprefix, 'spikes',day);
    
    %load the tetinfo and cellinfo file
    tetinfo = loaddatastruct(directoryname,fileprefix,'tetinfo');
    cellinfo = loaddatastruct(directoryname,fileprefix,'cellinfo');
    
    for epoch = 1:min(7,length(spikes{day}))
        if ~isempty(spikes{day}{epoch})
            
            %Go through each tetrode that has cells on it, 
            %assign both the global ripple phase to each spike
            cells = evaluatefilter(cellinfo{day}{epoch},'isequal($area,''CA1'') | isequal($area,''CA3'') & $numspikes > 0');
            tets = unique(cells(:,1));

            for t = tets'
                %Load the ripple file
                loadfile = sprintf('%s/EEGnonreference/%seeg%02d-%d-%02d.mat', directoryname, fileprefix,day,epoch,t);
                load(loadfile)
                ltime = geteegtimes(eeg{day}{epoch}{t});
                lphase = filtereeg2(eeg{day}{epoch}{t},ripplefilter);
                lphase = lphase.data(:,2);
                clear ripple loadfile

                %Go through each cell on the tetrode and assign local ripple phase
                for c = cells(cells(:,1)==t,2)'
                    if ~isempty(spikes{day}{epoch}{t}{c}.data)
                        s = spikes{day}{epoch}{t}{c}.data(:,1);
                        ind = lookup(s,ltime);
                        spikes{day}{epoch}{t}{c}.ripplephase = lphase(ind);
                    end
                end
            end
        end
    end
    
    %Save spikes file
    savefile = sprintf('%s/%sspikes%02d.mat', directoryname, fileprefix,day);
    save(savefile,'spikes');
        
end