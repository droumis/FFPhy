function toppeak = jy_spikeprocess(animdirect, daydirect, fileprefix, daynum, varargin)

%It will process the matclust files for the desired
%day, and put the files in the ANIMDIRECT directory.
%
%based on parts of DAYPROCESS(DAYDIRECT, ANIMDIRECT, VARPREFIX, FILEPREFIX, DAYNUM, OPTIONS)
%
%DAYDIRECT  -- folder with raw data files in /data14/jai/ 
%              no extra / needed eg. 'H2'
%ANIMDIRECT -- folder with the animals processed data in /data14/jai/ 
%              no extra / needed eg. 'H2_'
%FILEPREFIX -- name of animal eg. 'H2'
%DAYNUM     -- the day number for the experiment (starting with 1)
%OPTIONS    -- optional input in the form of: 'option', value, 'option', value
%          
%           SYSTEM -- 1 for old rig, 2 for nspike (default 2)
%           VARPREFIX -- the first three letters of the animals name, which is attached to all variable names.
%              I recommend putting '' for this, which will not include any name in the variables. (Default '')

currdir = pwd;

% set root directory
rootdir='/data14/jai/';

% set default values for the optional arguments
processeeg = 1;
system = 2;
varprefix = '';

% replace default values with any supplied ones
for option = 1:2:length(varargin)-1
    switch varargin{option}   
        case 'processeeg'
            processeeg = varargin{option+1};
        case 'system'
            system = varargin{option+1};
        case 'varprefix'
            varprefix = varargin{option+1};       
        otherwise
            error(['Option ''', varargin{option}, ''' not defined']);
    end
end

% specify raw and processed data directories
if (daynum < 10)
   daystring = ['0',num2str(daynum)];
else
   daystring = num2str(daynum);
end

daydirect=strcat(rootdir,daydirect,'/',fileprefix,daystring,'/');

animdirect=strcat(rootdir,animdirect,'/');

cd(animdirect)

eegdir = dir('EEG');
if (isempty(eegdir))
   %an EEG folder needs to be created
   !mkdir EEG
end

% get time file for day
tfilename = strcat(daydirect,'times.mat');
times = gettimes(tfilename);
runepochs = [];
eeg = [];

% get pos file
posfilename = strcat(animdirect,fileprefix,'data',daystring,'.mat');
load(posfilename);

% process clustered cells in each of the tetrode directories

cd(daydirect);         
spikes = [];  
   
tetfolders = dir;
for tet = 3:length(tetfolders)
  cd(currdir);
  cd(daydirect);   
  if tetfolders(tet).isdir
     dashfind = strfind(tetfolders(tet).name,'-');
     depth = str2num(tetfolders(tet).name(dashfind+1:end));
     tetnum = str2num(tetfolders(tet).name(1:dashfind-1));
     cd(tetfolders(tet).name);
     matclustfile = dir('matclust_*');
     
     if ~isempty(matclustfile)
        for epoch = 1:length(times)  
            if ~isempty(times(epoch).starttime)
                starttime = times(epoch).starttime;
                endtime = times(epoch).endtime;
                epochname = times(epoch).name;
                clear clustattrib clustdata
                load(matclustfile(1).name);
                starttime2 = timetrans({starttime},10000,2);
                endtime2 = timetrans({endtime},10000,2);
                if ~isempty(clustattrib.clusters)
                    for clustnum = 1:length(clustattrib.clusters)
                        if (~isempty(clustattrib.clusters{clustnum}))
                            %make sure that the cluster was defined for the current time epoch.  If not, make the cell empty.
                            if (is_cluster_defined_in_epoch(clustnum,epoch+1))
                                timepoints = clustdata.params(clustattrib.clusters{clustnum}.index,1);
                                if (size(clustdata.params,2) > 4)
                                    amps = clustdata.params(clustattrib.clusters{clustnum}.index,2:5);
                                else
                                    amps = nan(size(clustdata.params,1),1);
                                end
                                if (size(clustdata.params,2) > 5)
                                    spikewidth = mean(clustdata.params(clustattrib.clusters{clustnum}.index,6));
                                else
                                    spikewidth = nan(size(clustdata.params,1),1);
                                end
                                   
                                timepoints = timepoints(find((timepoints >= starttime2) & (timepoints <= endtime2)));
                                amps = amps(find((timepoints >= starttime2) & (timepoints <= endtime2)),:);
                                ampvar = var(amps);
                                [trash, maxvar] = max(ampvar);
                                amps = amps(:,maxvar);
                                timepoints = timepoints(:);
                                
                                
                                
                                if ~isempty(data{daynum}{epoch}.Pos.correcteddata)
                                    [spikepos, posindex] = lookuptime3(timepoints/10000, data{daynum}{epoch}.Pos.correcteddata(:,1),data{daynum}{epoch}.Pos.correcteddata(:,2:4)');
                                    findgoodpoints = find((spikepos(1,:)~=0)&(spikepos(2,:)~=0));
                                    spikepos = spikepos(:,findgoodpoints);
                                    posindex = posindex(findgoodpoints);
                                    timepoints = timepoints(findgoodpoints);
                                    amps = amps(findgoodpoints);
                                    if ~isempty(timepoints)

                                        spikes{daynum}{epoch}{tetnum}{clustnum}.data = timepoints/10000;
                                        spikes{daynum}{epoch}{tetnum}{clustnum}.data(:,2:4) = spikepos';
                                        spikes{daynum}{epoch}{tetnum}{clustnum}.data(:,6) = amps;
                                        spikes{daynum}{epoch}{tetnum}{clustnum}.data(:,7) = posindex;
                                        spikes{daynum}{epoch}{tetnum}{clustnum}.descript = 'spike data';
                                        spikes{daynum}{epoch}{tetnum}{clustnum}.fields = 'time x y dir not_used amplitude(highest variance channel) posindex';
                                        spikes{daynum}{epoch}{tetnum}{clustnum}.depth = depth;
                                        spikes{daynum}{epoch}{tetnum}{clustnum}.spikewidth = spikewidth;
                                        spikes{daynum}{epoch}{tetnum}{clustnum}.timerange = [starttime2 endtime2];
                                    else  %if the epoch was defined for the cluster, but no valid points were inside the boxes, make the data field empty
                                        spikes{daynum}{epoch}{tetnum}{clustnum}.data = [];
                                        spikes{daynum}{epoch}{tetnum}{clustnum}.descript = 'spike data';
                                        spikes{daynum}{epoch}{tetnum}{clustnum}.fields = 'time x y dir not_used amplitude(highest variance channel) posindex';
                                        spikes{daynum}{epoch}{tetnum}{clustnum}.depth = depth;
                                        spikes{daynum}{epoch}{tetnum}{clustnum}.spikewidth = spikewidth;
                                        spikes{daynum}{epoch}{tetnum}{clustnum}.timerange = [starttime2 endtime2];
                                    end
                                else %no valid pos info, so just save the spikes
                                    spikes{daynum}{epoch}{tetnum}{clustnum}.data = timepoints;
                                    spikes{daynum}{epoch}{tetnum}{clustnum}.descript = 'spike data';
                                    spikes{daynum}{epoch}{tetnum}{clustnum}.fields = 'time';
                                    spikes{daynum}{epoch}{tetnum}{clustnum}.depth = depth;
                                    spikes{daynum}{epoch}{tetnum}{clustnum}.spikewidth = spikewidth;
                                    spikes{daynum}{epoch}{tetnum}{clustnum}.timerange = [starttime2 endtime2];
                                end
                            else
                                spikes{daynum}{epoch}{tetnum}{clustnum} = [];
                            end
                        end
                    end
                end
            end
        end
     end
     cd ..
  end
end

% save files
    cd(animdirect);
    eval([varprefix,'spikes = spikes']);
    eval(['save ',fileprefix,'spikes',daystring,' ',varprefix,'spikes']);
   
    cd(currdir);
end
    
