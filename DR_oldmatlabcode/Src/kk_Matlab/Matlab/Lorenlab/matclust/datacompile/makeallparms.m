function makeallparms(thresh, maxallowedamp, system)
%
%makeallparms 
%makeallparms(thresh)
%makeallparms(thresh, maxallowedamp)
%makeallparms(thresh, maxallowedamp, system)
%run this program from the directory containing the folders for all the days
%it will read the .tt files and create a paramter .m file in the same subdirectory as the .tt file
%it will also create the m files containing the waveform info and the position info
%thresh - the threshold that at least one spike must excede in order to be included (default 0 microvolts)
%maxallowedamp - if the amplitude on any channel excedes this level (in microvolts), exclude the spike. Default 2500.
%system = 1 for the old rig with spike
%system = 2 for the new rigs with nspike (default)


if (nargin < 3)
    system = 2;  %threshhold: at least one spike must excede this amplitude (in micro volts)
end
if (nargin < 2)
    maxallowedamp = 2500;  %maximum allowed amplitude in micro volts.  This filters out noise and makes the graphing funtions faster 
end
if (nargin == 0)
    thresh = 0;
end



dayfolders = dir;
currdir = pwd;
for i = 3:length(dayfolders)
   if dayfolders(i).isdir
      disp(upper(dayfolders(i).name))
      cd(dayfolders(i).name);
      tetfolders = dir;
      for j = 3:length(tetfolders)
         if tetfolders(j).isdir
            cd(tetfolders(j).name);
            ttfile = dir('*.tt');
            if ~isempty(ttfile)
               disp(['      ',ttfile(1).name])
               dashfind = strfind(ttfile(1).name,'-');
               tetstring = ttfile(1).name(1:dashfind-1);
               [timestamps, waves] = readtt(ttfile(1).name);
               save([dayfolders(i).name,'-',tetstring],'timestamps', 'waves');
               if (~isempty(timestamps))
                  if (system == 1)
                     makeparams([dayfolders(i).name,'-',tetstring],str2num(tetstring));
                  elseif (system == 2)
                     nmakeparams([dayfolders(i).name,'-',tetstring],str2num(tetstring),waves,timestamps,thresh,maxallowedamp);
                  end             
               end
               clear tmp1;
            end
            cd ..
         end
      end
      cd(currdir);
   end
end   
