function sj_makedayparms_withpc(dayname,thresh, maxallowedamp, varargin)
%
% Shantanu - May 2012. Dont use this! Calls sj_nmakeparams_withpc.m
% PCs calculated without de-jitterring
% Insread, use sj_makedayparms_pc.m


%sj_makedayparms_withpc(dayname) 
%makedayparms(dayname,thresh)
%makedayparms(dayname,thresh, maxallowedamp)
%makedayparms(dayname,thresh, maxallowedamp, options)
%run this program from the directory containing the folders for all the days
%it will read the .tt files and create a paramter .m file in the same subdirectory as the .tt file
%it will also create the m files containing the waveform info and the position info
%thresh - the threshold that at least one spike must excede in order to be included (default 0 microvolts)
%maxallowedamp - if the amplitude on any channel excedes this level (in microvolts), exclude the spike. Default 2500.
% 
%Options
%'system' -- 1 for old rig, 2 for new rig (default)
%'pos' -- 0 for no position, 1 for position (default)

system = 3;   % default system = 3 for PC calculation and Klustakwik file generation
dogspikefile = 1; % Save for gsort (in addition to Klustakwik, which is default)
Samprange = [4:32];  %Important parameter: look at subset of spike sample for PC 

calcpos = 1;

% Manual Noise Removel 
donoise=1;
% To de-jitter waveforms before PC calculation for Klustakwik 
dejitter=0;

if (nargin > 3)
   for option = 1:2:length(varargin)-1
      switch varargin{option}
         case 'system'
            system = varargin{option+1}; 
         case 'pos'
            calcpos = varargin{option+1};
        
         otherwise
            error(['Option', varargin{option},'not defined']);
        
      end
   end
end


if (nargin < 3)
    maxallowedamp = 2500;  %maximum allowed amplitude in micro volts.  This filters out noise and makes the graphing funtions faster 
end
if (nargin < 2)
    thresh = 0;
end



dayfolders = dir;
currdir = pwd;
for i = 3:length(dayfolders)
   if dayfolders(i).isdir
      if strcmp(dayname,dayfolders(i).name)
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
                          if calcpos
                              nmakeparams([dayfolders(i).name,'-',tetstring],str2num(tetstring),waves,timestamps,thresh,maxallowedamp);
                          else
                              nmakeparamsNoPos([dayfolders(i).name,'-',tetstring],str2num(tetstring),waves,timestamps,thresh,maxallowedamp);
                          end
                      elseif (system == 3)  % for PC calc
                          if calcpos 
                              sj_nmakeparams_withpc([dayfolders(i).name,'-',tetstring],str2num(tetstring),waves,timestamps,thresh,maxallowedamp,dogspikefile, Samprange, donoise, dejitter);
                          else
                              sj_nmakeparamsNoPos_withpc([dayfolders(i).name,'-',tetstring],str2num(tetstring),waves,timestamps,thresh,maxallowedamp,dogspikefile, Samprange, donoise, dejitter);
                          end % end for PC calc
                      end
                  end
                  clear tmp1;
               end
               cd ..
            end
         end
      end
      cd(currdir);
   end
end   
