function toppeak = sj_dayprocess_noeeg(daydirect, animdirect, fileprefix, daynum, varargin)


% Shantanu 07 Jan 2009
% Using dayprocess.m from /home/loren/Src/Matlab
% Adding multiunit and dio processin from Lorens and Calebs code in
% /home/loren/Src/Matlab/diomatlab

%%

%Run this program from the folder containing all of the animal's day folders.  It will process the matclust files for the desired
%day, and put the files in the ANIMDIRECT directory.
%
%DAYPROCESS(DAYDIRECT, ANIMDIRECT, VARPREFIX, FILEPREFIX, DAYNUM, OPTIONS)
%
%DAYDIRECT -- folder name where the day's data is stored
%ANIMDIRECT -- the path to where the animal's processed data will be stored -- example '/data99/student/Mil'
%FILEPREFIX -- also the first three letters of the animal's name (example 'con'), which will attach to
%              the beginning of the .mat files containing the variables.  I recommend using this.
%              If not, type ''.
%DAYNUM -- the day number for the experiment (starting with 1)
%OPTIONS -- optional input in the form of: 'option', value, 'option', value
%           optional arguments are:
%           DIODEPOS -- 0 uses back diode for pos, 1 uses front diode for
%               pos, and values in between use the proportionate distance
%               between the two diodes. (default 0.5, ie, average of diodes)
%           CMPERPIX -- size of each pixel in centimeters (must be specified)
%           POSFILT -- filter to use for pos smoothing when computing
%               velocity (default gaussian(30*0.5, 60))
%           PROCESSEEG -- 1 to process EEG, 0 to skip EEG processing (default 1)
%           SYSTEM -- 1 for old rig, 2 for nspike (default 2)
%           VARPREFIX -- the first three letters of the animals name, which is attached to all variable names.
%              I recommend putting '' for this, which will not include any name in the variables. (Default '')
%           REVERSEX -- mirror image pos info along x axis (default 0)
%           REVERSEY -- mirror image pos info along y axis (default 0)

%%

% set default values for the optional arguments
diodepos = 0.5;
cmperpix = NaN;
processeeg = 0;
system = 2;
posfilt = gaussian(30*0.5, 60);
varprefix = '';
reversex = 0;
reversey = 0;
processdio = 1;


% replace default values with any supplied ones
for option = 1:2:length(varargin)-1
    
    switch varargin{option}
        case 'diodepos'
            diodepos = varargin{option+1};
        case 'cmperpix'
            cmperpix = varargin{option+1};          
        case 'processeeg'
            processeeg = varargin{option+1};
        case 'system'
            system = varargin{option+1};
        case 'posfilt'
            posfilt = varargin{option+1};
        case 'varprefix'
            varprefix = varargin{option+1};    
        case 'reversex'
            reversex = varargin{option+1};
        case 'reversey'
            reversey = varargin{option+1};    
        
        otherwise
            error(['Option ''', varargin{option}, ''' not defined']);
    end
end

if (~isfinite(cmperpix))
    error('cmperpix must be specified');
end

[daypath,dayname] = fileparts(daydirect);

%%

currdir = pwd;
if (animdirect(end) == '/')
   animdirect = animdirect(1:end-1);
end
cd(animdirect)
eegdir = dir('EEG');
if (isempty(eegdir))
   %an EEG folder needs to be created
   !mkdir EEG
end
cd(currdir)

NUMTETRODES = 30;
currdir = pwd;
cd(daydirect);
times = gettimes('times');
runepochs = [];
eeg = [];

if (daynum < 10)
   daystring = ['0',num2str(daynum)];
else
   daystring = num2str(daynum);
end

%%

for epoch = 1:length(times)
    if ~isempty(times(epoch).starttime)
        starttime = times(epoch).starttime;
        endtime = times(epoch).endtime;
        epochname = times(epoch).name;
        % read in the position file
        posfiles = dir('*.p');


        % read in the raw diode information (front and back positions)
        % this looks through all the .p files and picks out the times that are
        % within the start and end times
        rawpos{daynum}{epoch} = readrawpos(posfiles, starttime,endtime);

        % interpolate over the raw positions to get location and direction
        % at each time
        if ~isempty(rawpos{daynum}{epoch}.data)
            pos{daynum}{epoch} = posinterp(rawpos{daynum}{epoch}, 'diode', diodepos,...
                'maxv', 300, 'maxdevpoints', 30, 'reversex', reversex, 'reversey', reversey);

            % multiply pos data by cmperpix to convert from pixel units to cm
            pos{daynum}{epoch}.cmperpixel = cmperpix;
            pos{daynum}{epoch}.data(:,2:3) = pos{daynum}{epoch}.data(:,2:3)*cmperpix;

            % now run functions to add additional behavior parameters to pos struct
            pos{daynum}{epoch} = addvelocity(pos{daynum}{epoch}, posfilt);
        else
            pos{daynum}{epoch}.data = [];
            pos{daynum}{epoch}.cmperpixel = cmperpix;
        end

        if (strncmp(epochname, 'run', 3))
            runepochs = [runepochs epoch];
            % assign the environment and task variables
        end
    end
end

%%

% get a list of all the eeg files in the directory
eegfiles = dir('*.eeg');
    
if (processeeg) %this part filters the eeg data 
    if ~isempty(eegfiles) % if there are eeg files to process, then proceed
        for eegfile = 1:length(eegfiles)
            % get the tetrode number and depth, which is contained in the
            % name of the eeg file
            eegfilename = eegfiles(eegfile).name;
            dashfind = strfind(eegfilename,'-');
            pointfind = strfind(eegfilename,'.');
            depth = str2num(eegfilename(dashfind(1)+1:pointfind-1));
            eegtet = str2num(eegfilename(1:dashfind-1));
            % prepare a string label for this eeg data
            if (eegtet < 10)
                eegstring = ['0',num2str(eegtet)];
            else
                eegstring = num2str(eegtet);
            end
            % now get the eeg data for each epoch
            for epoch = 1:length(times)
                if ~isempty(times(epoch).starttime)
                    ['epoch ', num2str(epoch),' tetrode ',num2str(eegtet)]
                    starttime = times(epoch).starttime;
                    endtime = times(epoch).endtime;
                    epochname = times(epoch).name;
                    if (system == 1)
                        eegstruct = readeeg(eegfilename, starttime, endtime, eegtet/NUMTETRODES);
                    elseif (system == 2)
                        eegstruct = nreadeeg(eegfilename, starttime, endtime);
                    end
                    eegstruct.depth = depth;
                    eval([varprefix,'eeg{',num2str(daynum),'}{',num2str(epoch),'}{',num2str(eegtet),'} = eegstruct;']);
                    tmpfilename = [animdirect,'/EEG/',fileprefix,'eeg',daystring,'-',num2str(epoch),'-',eegstring];
                    eval(['save ',tmpfilename, ' ',varprefix,'eeg']);
                    eval(['clear ',varprefix,'eeg']);
                end
            end
        end
    end
end

cd(currdir);
cd(daydirect);         
spikes = [];  
multi = [];
   

%% Multiunit Data - From Caleb and Loren's code

tetfolders = dir;
tetfolders(~[tetfolders.isdir]) = [];
% cd(daydirect);

tetfolders = {tetfolders.name};
tetfolderinfo = regexp(tetfolders,'(?<tetrode>\d\d)-(?<depth>\d\d\d)','names');

for tet = 1:length(tetfolders)
    if ~isempty(tetfolderinfo{tet})
        depth = str2num(tetfolderinfo{tet}.depth);
        tetnum = str2num(tetfolderinfo{tet}.tetrode);
        
        if exist(fullfile(tetfolders{tet},sprintf('%s-%02d_params.mat',dayname,tetnum)))
            rawdata = load(fullfile(tetfolders{tet},sprintf('%s-%02d_params.mat',dayname,tetnum)));
            if isempty(rawdata)
                error(sprintf('Raw param file not found for tetrode %d.',tetnum));
            end
            for epoch = 1:length(times)
                if ~isempty(times(epoch).starttime)
                    starttime = timetrans({times(epoch).starttime},10000,2);
                    endtime = timetrans({times(epoch).endtime},10000,2);
                    epochname = times(epoch).name;
                    multi{daynum}{epoch}{tetnum} = rawdata.filedata.params((find((rawdata.filedata.params(:,1) >= starttime) & (rawdata.filedata.params(:,1) <= endtime))),1);
                end
            end
        else
            multi{daynum}{epoch}{tetnum} = [];
        end
        
    end
end


%%


%%%%%% Spikes from ori dayprocess /home/loren/Src/Matlab


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
                                if ~isempty(pos{daynum}{epoch}.data)
                                    [spikepos, posindex] = lookuptime3(timepoints/10000, pos{daynum}{epoch}.data(:,1),pos{daynum}{epoch}.data(:,2:4)');
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




cd(animdirect);
eval([varprefix,'rawpos = rawpos;']);
eval(['save ',fileprefix,'rawpos',daystring,' ', varprefix,'rawpos']);
eval([varprefix,'pos = pos;']);
eval(['save ',fileprefix,'pos',daystring,' ',varprefix,'pos']);
eval([varprefix,'spikes = spikes;']);
eval(['save ',fileprefix,'spikes',daystring,' ',varprefix,'spikes']);

eval([varprefix,'multi = multi;']);
eval(['save ',fileprefix,'multi',daystring,' ',varprefix,'multi']);
%save(fullfile(animdirect,sprintf('%smulti%02d.mat',fileprefix,daynum)),[varprefix,'multi']);

cd(currdir);
%%

%if (processdio)
%    creatediostruct(daydirect,animdirect,fileprefix,daynum);
%end
%%

    
