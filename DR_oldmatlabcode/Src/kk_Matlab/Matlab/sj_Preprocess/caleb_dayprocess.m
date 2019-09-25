function toppeak = dayprocess(daydirect, animdirect,animprefix,fileprefix, daynum, varargin)

%Run this program from the folder containing all of the animal's day folders.  It will process the matclust files for the desired
%day, and put the files in the ANIMDIRECT directory.
%
%DAYPROCESS(DAYDIRECT, ANIMDIRECT, ANIMPREFIX, FILEPREFIX, DAYNUM)
%DAYPROCESS(DAYDIRECT, ANIMDIRECT, ANIMPREFIX, FILEPREFIX, DAYNUM, PROCESSEEG) -- set PROCESSEEG to 0 to skip eeg processing (default 1)
%DAYPROCESS(DAYDIRECT, ANIMDIRECT, ANIMPREFIX, FILEPREFIX, DAYNUM, PROCESSEEG, SYSTEM) -- SYSTEM is 1 for old rig, 2 for nspike (default 2)
%
%
%DAYDIRECT -- folder name where the day's data is stored
%ANIMDIRECT -- the path to where the animal's processed data will be stored -- example '/data99/student/Mil'
%ANIMPREFIX -- the first three letters of the animals name, which is attached to all variable names.
%              I recommend putting '' for this, which will not include any name in the variables.
%FILEPREFIX -- also the first three letters of the animal's name (example 'con'), which will attach to
%              the beginning of the .mat files containing the variables.  I recommend using this.
%              If not, type ''.
%DAYNUM -- the day number for the experiment (starting with 1)




if ~isempty(varargin)
   processeegfiles = varargin{1};
   if length(varargin) > 1
      systemused = varargin{2};
   else
      systemused = 2;
   end

else   
   processeegfiles = 1; %set to 1 if eeg files need to be read and filtered to theta files. It should be zero if this has already been done.
   systemused = 2;
end
NUMTETRODES = 30;
load thetafilter;
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

for epoch = 1:length(times)
   starttime = times(epoch).starttime;
   endtime = times(epoch).endtime;
   epochname = times(epoch).name;
   % read in the position file 
   posfile = dir('*.p');
   pnum = 1;
   if (length(posfile) > 1)
         % check for a .p file of the form 'animal##_n.p' where n is the
         % current epoch
         epochstr = sprintf('_%d', epoch);
         for j = 1:length(posfile)
            if (~isempty(findstr(posfile(j).name, epochstr)))
               pnum = j;
               break;
            elseif (isempty(findstr(posfile(j).name, '_')))
               % this is the base file, so it is the default
               pnum = j;
            end
         end
         %warning(sprintf('multiple .p files. Using %s', ...
         %      posfile(pnum).name))
   end
   % read in the raw diode information (front and back positions)
   
   
   %load([animdirect,'/',fileprefix,'rawpos',daystring]);
   %rawpos{daynum}{epoch} = readpos(posfile(pnum).name, starttime, endtime);
   rawpos{daynum}{epoch} = readrawpos(posfile(pnum).name, starttime,endtime);
   
   % interpolate over the raw positions to get location and direction
   % at each time
   pos{daynum}{epoch} = posinterp(rawpos{daynum}{epoch}, ...
               'pinterp', 1, 'dinterp', 1, 'diode', 0.5, 'maxv', 300, 'reversex', 1);
   if (strncmp(epochname, 'run', 3))
         runepochs = [runepochs epoch];
         % assign the environment and task variables
   end
 
 end

eegfiles = dir('*.eeg');

if (processeegfiles) %this part filters the eeg data and saves the theta files (which takes some time)
   if ~isempty(eegfiles)
   for eegfile = 1:length(eegfiles)
      eegfilename = eegfiles(eegfile).name;
      dashfind = strfind(eegfilename,'-');
      pointfind = strfind(eegfilename,'.');
      depth = str2num(eegfilename(dashfind(1)+1:pointfind-1));
      eegtet = str2num(eegfilename(1:dashfind-1));
      if (eegtet < 10)
         eegstring = ['0',num2str(eegtet)];
      else
         eegstring = num2str(eegtet);
      end  
      for epoch = 1:length(times)
         ['epoch ', num2str(epoch),' tetrode ',num2str(eegtet)]
         starttime = times(epoch).starttime;
         endtime = times(epoch).endtime;
         epochname = times(epoch).name;    
         if (systemused == 1)
            eegstruct = readeeg(eegfilename, starttime, endtime, eegtet/NUMTETRODES);
         elseif (systemused == 2)
            eegstruct = nreadeeg(eegfilename, starttime, endtime);
         end
         eegstruct.depth = depth;
         eval([animprefix,'eeg{',num2str(daynum),'}{',num2str(epoch),'}{',num2str(eegtet),'} = eegstruct;']);
         tmpfilename = [animdirect,'/EEG/',fileprefix,'eeg',daystring,'-',num2str(epoch),'-',eegstring];
         eval(['save ',tmpfilename, ' ',animprefix,'eeg']);
               
         eval([animprefix,'theta{',num2str(daynum),'}{',num2str(epoch),'}{',num2str(eegtet),'} = applyfilter(',animprefix,'eeg, [',num2str(daynum),' ',num2str(epoch),' ',num2str(eegtet),'], thetafilter, 200);']);
         eval([animprefix,'theta{',num2str(daynum),'}{',num2str(epoch),'}{',num2str(eegtet),'}.peaks = findeegpeaks(',animprefix,'theta, [',num2str(daynum),' ',num2str(epoch),' ',num2str(eegtet),']);']);
         tmpfilename = [animdirect,'/EEG/',fileprefix,'theta',daystring,'-',num2str(epoch),'-',eegstring];
         eval(['thetapeaks = ',animprefix,'theta{',num2str(daynum),'}{',num2str(epoch),'}{',num2str(eegtet),'}.data(',animprefix,'theta{',num2str(daynum),'}{',num2str(epoch),'}{',num2str(eegtet),'}.peaks);']);
         meanpeak(epoch,eegtet) = mean(thetapeaks);
         eval(['save ',tmpfilename, ' ',animprefix,'theta']);
         eval(['clear ',animprefix,'theta']);
         eval(['clear ',animprefix,'eeg']);
   
      end
   end
   end
   meanpeak = mean(meanpeak);
end
pushd(animdirect);
try
   %if the theta electrode was already defined, then use that
   load([fileprefix,'thetaelect',daystring]);
   toppeak = thetaelect;
catch
   %otherwise, find the tetrode with the largest average theta
   [trash,toppeak] = max(meanpeak);
   
end
popd;
if (toppeak < 10)
   eegstring = ['0',num2str(toppeak)];
else
   eegstring = num2str(toppeak);
end  
for epoch = 1:length(times)
    tmpfilename = [animdirect,'/EEG/',fileprefix,'theta',daystring,'-',num2str(epoch),'-',eegstring];

    eval(['load ',tmpfilename]);
    eval(['thetatet(',num2str(epoch),') = ',animprefix,'theta{',num2str(daynum),'}{',num2str(epoch),'}{',num2str(toppeak),'};']);
    
    eval(['clear ',animprefix,'theta']);
end

cd(currdir);
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
     load([animdirect,'/recordingareas']); 
     if ~isempty(matclustfile)
        for epoch = 1:length(times)  
           starttime = times(epoch).starttime;
           endtime = times(epoch).endtime;
           epochname = times(epoch).name;
           load(matclustfile(1).name);
             % generate the list of times corresponding to the peaks of the theta rhythm
           starttime2 = timetrans({starttime},10000,2);
           endtime2 = timetrans({endtime},10000,2);
           if ~isempty(clustattrib.clusters)
               for clustnum = 1:length(clustattrib.clusters)                  
                   if (~isempty(clustattrib.clusters{clustnum}))  
                       %make sure that the cluster was defined for the current time epoch.  If not, make the cell empty.
                       if (is_cluster_defined_in_epoch(clustnum,epoch+1))  
                           timepoints = clustdata.params(clustattrib.clusters{clustnum}.index,1);
                           amps = clustdata.params(clustattrib.clusters{clustnum}.index,2:5);
                           
                           timepoints = timepoints(find((timepoints >= starttime2) & (timepoints <= endtime2)));
                           amps = amps(find((timepoints >= starttime2) & (timepoints <= endtime2)),:);
                           ampvar = var(amps);
                           [trash, maxvar] = max(ampvar);
                           amps = amps(:,maxvar);
                           timepoints = timepoints(:);
                           [spikepos, posindex] = lookuptime3(timepoints/10000, pos{daynum}{epoch}.data(:,1),pos{daynum}{epoch}.data(:,2:4)');
                           findgoodpoints = find((spikepos(1,:)~=0)&(spikepos(2,:)~=0)&(spikepos(3,:)~=0));
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
                               spikes{daynum}{epoch}{tetnum}{clustnum}.fields = 'time x y dir thetaphase amplitude(highest variance channel) posindex';
                               spikes{daynum}{epoch}{tetnum}{clustnum}.depth = depth;
                           else  %if the epoch was defined for the cluster, but no valid points were inside the boxes, make the data field empty
                               spikes{daynum}{epoch}{tetnum}{clustnum}.data = [];
                               spikes{daynum}{epoch}{tetnum}{clustnum}.descript = 'spike data';
                               spikes{daynum}{epoch}{tetnum}{clustnum}.fields = 'time x y dir thetaphase amplitude(highest variance channel) posindex';
                               spikes{daynum}{epoch}{tetnum}{clustnum}.depth = depth;
                               spikes{daynum}{epoch}{tetnum}{clustnum}.location = recordingareas{daynum,tetnum};
                           end
                           %[epoch tetnum clustnum length(find(spikepos(1,:)==0))]
                           eegtimes = thetatet(epoch).peaks / thetatet(epoch).samprate + thetatet(epoch).starttime;
                           
                           if ~isempty(timepoints)
                               % make each peak correspond to an integer and interpolate spikes into the
                               % list of peaks. The resulting decimal fraction represents the distance
                               % between adjacent peaks and can be converted into a theta phase by
                               % multiplying by 2pi
                               % position locations
                               
                               phase = interp1(eegtimes, 1:length(thetatet(epoch).peaks), timepoints'/10000, 'linear');
                               phase = (phase - floor(phase)) * 2 * pi;
                               spikes{daynum}{epoch}{tetnum}{clustnum}.data(:,5) = phase;
                               spikes{daynum}{epoch}{tetnum}{clustnum}.location = recordingareas{daynum,tetnum};
                           end                   
                       else
                           spikes{daynum}{epoch}{tetnum}{clustnum} = [];
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
    
    thetaelect = toppeak;
    eval(['save ',fileprefix,'thetaelect',daystring,' thetaelect']);
    eval([animprefix,'rawpos = rawpos']);
    eval(['save ',fileprefix,'rawpos',daystring,' ', animprefix,'rawpos']);
    eval([animprefix,'pos = pos']);
    eval(['save ',fileprefix,'pos',daystring,' ',animprefix,'pos']);
    eval([animprefix,'spikes = spikes']);
    eval(['save ',fileprefix,'spikes',daystring,' ',animprefix,'spikes']);
   
    cd(currdir);
    
