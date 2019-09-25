function sj_makedayparms_pc2(dayfolder ,thresh, maxallowedamp, varargin)
%
% Shantanu - ver 2. Changing argument calling to make it equivalent to makedayparms. 
% Give the day folder as argument instead of searching for it

% Shantanu - May 2012.
% From sj_makedayparms_dio.m - Getting rid of additional parameters and only keeping PCs

% eg. sj_makedayparms_pc('/data25/sjadhav/HPExpt/JW7/','JW7',4);


%run this program from the directory containing the folders for all the days
%it will read the .tt files and create a paramter .m file in the same subdirectory as the .tt file
%it will also create the m files containing the waveform info and the position info
%thresh - the threshold that at least one spike must excede in order to be included (default 0 microvolts)
%maxallowedamp - if the amplitude on any channel excedes this level (in microvolts), exclude the spike. Default 2500.
%
%Options
%'system' -- 1 for old rig, 2 for new rig (default), 3 for PC calculation
%'pos' -- 0 for no position, 1 for position (default)

system = 3;   % default system = 3 for Additonal Params and PC calculation and Klustakwik file generation
Samprange = 4:28;  %Parameter: look at subset of spike sample 
calcpos = 1;

dogspikefile = 0; % Save for gsort (in addition to Klustakwik, which is default)
% DIO removal
dodio = 0;
% Manual Noise Removel
donoise = 0;
% Do PC
dopc = 1;

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

if (nargin < 4)
    thresh = 0;
end

if (nargin < 5)
    maxallowedamp = 1000;  %maximum allowed amplitude in micro volts.  This filters out noise and makes the graphing funtions faster
end

currdir = pwd;
cd(dayfolder);
disp(upper(dayfolder));

if dodio==1
    % Load DIO file for day and get stimtimes
    % Window for stimulation artifact removal  - if dodio=1;
    pret=0; postt=5; % ms
    
    DIOfile = sprintf('%s/%sDIO%02d.mat', directdir, prefix, day);
    load(DIOfile);
    nepochs = size(DIO{day},2);
    if isempty(nepochs) || nepochs==0
        error('DIO file has no epochs!');
    end
    stim_time = [];
    for epoch=1:nepochs
        stim = DIO{day}{epoch}{15};
        if isempty(stim)
            stim = DIO{day}{epoch}{16};
        end
        stim_time = [stim_time; stim.pulsetimes(:,1)]; % This is in spike ts units - 1ms=10000 points (100 us resolution)
    end
end

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
            
            if dodio==1
                disp('      Removing Stimulation Artifacts');
                % Removal of stimulation artifacts by DIO file
                % Looping over stim_time is better that repmat as you could have ~1 million spikes
                artidxs=[];
                for s = 1:length(stim_time)
                    artidxs = [ artidxs; find( (timestamps>=(stim_time(s) - pret*10)) & (timestamps<=(stim_time(s) + postt*10)) )];
                end
                % If you want to plot and check the artifact waveforms
                %artwaves = waves(:,:,artidxs);
                %Nsamp = size(artwaves,1); Nch = size(artwaves,2); Nspk = size(artwaves,3);
                %artwaves = reshape(artwaves,Nsamp*Nch,Nspk);
                %figure; hold on; plot (artwaves,'r');
                timestamps(artidxs)=[]; waves(:,:,artidxs)=[];
            end
            
            save([dayfolder,'-',tetstring],'timestamps', 'waves');
            if (~isempty(timestamps))
                if (system == 1)
                    makeparams([dayfolder,'-',tetstring],str2num(tetstring));
                elseif (system == 2)
                    if calcpos
                        nmakeparams([dayfolder,'-',tetstring],str2num(tetstring),waves,timestamps,thresh,maxallowedamp);
                    else
                        nmakeparamsNoPos([dayfolder,'-',tetstring],str2num(tetstring),waves,timestamps,thresh,maxallowedamp);
                    end
                elseif (system == 3)  % for Addinal Params and PC calc
                    if calcpos
                        sj_nmakeparams_pc([dayfolder,'-',tetstring],str2num(tetstring),waves,timestamps,thresh,maxallowedamp, dopc, Samprange, donoise, dogspikefile);
                    else
                        sj_nmakeparamsNoPos_pc([dayfolder,'-',tetstring],str2num(tetstring),waves,timestamps,thresh,maxallowedamp, dopc, Samprange, donoise, dogspikefile);
                    end % end for PC calc
                end
            end
            clear tmp1;
        end
        cd ..
    end
end % for tetfolders






