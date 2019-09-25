function [Pos] = Posepochprocess(daydirect, animalname,day, epoch, varargin)

% Generates position data for the day
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

% 2011/11/18 - add recognition for 2 diodes

% set default values for the optional arguments
diodepos = 1; % position defined by the front diode to avoid non integer pixel values from taking a average of front and back diodes
cmperpix = NaN;
processeeg = 1;
system = 2;
posfilt = gaussian(30*0.5, 60);
varprefix = '';
reversex = 0;
reversey = 0;

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


cd(daydirect);
times = gettimes('times');
runepochs = [];


for epoch = epoch;
    if ~isempty(times(epoch).starttime)
        starttime = times(epoch).starttime;
        endtime = times(epoch).endtime;
        epochname = times(epoch).name;
        % read in the position file
        posfiles = dir('*.p');
        sprintf('Day %s Epoch %s',num2str(day), num2str(epoch))
        %posfile = sprintf('%s%s_%s.p',animalname,day,num2str(epoch));
        
        
        % read in the raw diode information (front and back positions)
        % this looks through all the .p files and picks out the times that are
        % within the start and end times
        rawpos = readrawpos(posfiles, starttime,endtime);
        
        if isempty(rawpos.data)
            error('no position data found in .p files');
        end
        
        %% Code to skip every other position (downsample) REMOVE WHEN DONE
        %         for aa = 1:50
        %         skipind = aa:51:size(rawpos.data,1);
        %             rawpos.data(skipind, 2:5) = repmat(0,length(skipind),4);
        %         end
        %%
        
        % perform Steve Kim's position reconstruction
        Rpos=rawpos.data(:,2:5);
        
%         % backup time stamp adjustment incase the wrong timestamp was used
%         % find out the offset between cputimestamp and postimestamp. The
%         % correct timestamp comes from postimestamp
%         dsz = '';
%         if (day < 10)
%             dsz = '0';
%         end
%         
%         filename=sprintf('%s%s%s.postimestamp', animalname,dsz,num2str(day));
%         fid=fopen(filename);
%         
%         headerstr = '';
%         while ~strncmp(headerstr, '%%ENDHEADER', 11)
%             headerstr = fgets(fid);
%             
%         end
%         
%         postime=fread(fid,'uint32');
%         
%         filename=sprintf('%s%s%s.cpupostimestamp', animalname,dsz,num2str(day));
%         fid=fopen(filename);
%         headerstr = '';
%         while ~strncmp(headerstr, '%%ENDHEADER', 11)
%             headerstr = fgets(fid);
%             
%         end
%         
%         cputime=fread(fid,'uint32');
%         
%         % match the times in rawpos with the corresponding time in postimes
%         
%         
%         
%         rawpos.data(:,1)=postime(lookup(rawpos.data(:,1),cputime/10000))/10000;
        %         % lens correction
        %         % for diode 1
        Rpos(Rpos(:,1)>320,1)=320;
        Rpos(Rpos(:,2)>240,2)=240;
        %         Rpos(Rpos(:,1)<1,1)=1;
        %         Rpos(Rpos(:,2)<1,2)=1;
        %
        %         % for diode 2
        Rpos(Rpos(:,3)>320,3)=320;
        Rpos(Rpos(:,4)>240,4)=240;
        %         Rpos(Rpos(:,3)<1,3)=1;
        %         Rpos(Rpos(:,4)<1,4)=1;
        %         min(Rpos)
        %         max(Rpos)
        
        pixind=pixelindex(Rpos(:,1:2),320,240);
        
        % find non out of bounds pixels
        inboundpix=find(pixind>0);
        coord=zeros(size(pixind,1),2);
        coord(inboundpix,:)=lenscorrtransf(pixind(inboundpix));
        Rpos(:,1:2)=coord;
        
        pixind=pixelindex(Rpos(:,3:4),320,240);
        inboundpix=find(pixind>0);
        coord=zeros(size(pixind,1),2);
        coord(inboundpix,:)=lenscorrtransf(pixind(inboundpix));
        Rpos(:,3:4)=coord;
        
        % if a 0 is in the position data for any x y points, we need to set
        % all the other x y points to 0 for the interpolation to work
        % correctly
        
        Rpos=[rawpos.data(:,1) Rpos(:,:)];
        
        % set out of bound
        
        %Rpos(((min(Rpos(:,2:5),[],2))==0),2:5)=0;
        %         zerorows=find(Rpos(:,2)<=0 | Rpos(:,3)<=0);
        %         Rpos(zerorows,2:5)=0; % values other than 0 are considered true positions and will not be correctly interpolated
        
        RRpos{1}.data=Rpos;
        
        % March 2013
        % shift pixels to avoid negative values 
        % need to avoid bad interpolation by SK_estimate_position
        % shift back later -- see below
        
        RRpos{1}.data(:,3)=RRpos{1}.data(:,3)+40;
        RRpos{1}.data(:,5)=RRpos{1}.data(:,5)+40;
        
        SKpos=[];
        %      if epoch~=5
        
        SKpos=SK_estimate_position(RRpos,'centimeters_per_pixel',cmperpix,'front_back_marker_weights',[1 0]);
        
        
        % shift pixels back to original position 
        SKpos{1}.data(:,3)=SKpos{1}.data(:,3)-40;
        
        if ~isempty(SKpos{1,1}.data(:,4));
        
        headdirection=SKpos{1,1}.data(:,4);
        
        else
            
            headdirection=zeros(size(SKpos{1,1}.data,1),1);
            
        end
        
        SKpos{1,1}.data(:,2:5)=SKpos{1,1}.data(:,2:5)/cmperpix;
        %       end
        % transform
        
        
        % need to use lens correction to correct for raw position
        
        % interpolate over the raw positions to get location and direction
        % at each time
        if ~isempty(rawpos.data)
            pos = posinterpcorrection(rawpos, 'diode', diodepos,...
                'maxv', 300, 'maxdevpoints', 30, 'reversex', reversex, 'reversey', reversey,'cmperpix',cmperpix);
            
            % fix out of bounds pixels since reconstruction program interpolates
            % beyond image boundary
            
            %             pos.data(pos.data(:,2)>320,1)=320;
            %             pos.data(pos.data(:,3)>240,2)=240;
            %             pos.data(pos.data(:,2)<1,1)=1;
            %             pos.data(pos.data(:,3)<1,2)=1;
            
            
            
            % multiply pos data by cmperpix to convert from pixel units to cm
            % columns 2 and 3: diode 1, columns 4 and 5: diode 2
            pos.cmperpixel = cmperpix;
            pos.data(:,2:3) = pos.data(:,2:3); %*cmperpix;
            
            
            
            % now run functions to add additional behavior parameters to pos struct
            pos = addvelocity(pos, posfilt);
        else
            pos.data = [];
            pos.cmperpixel = cmperpix;
        end
        
        if (strncmp(epochname, 'run', 3))
            runepochs = [runepochs epoch];
            % assign the environment and task variables
        end
    end
end

Pos={};
Pos.rawpos=rawpos;
Pos.pos=pos;
if ~isempty(SKpos)
    Pos.pos.data=SKpos{1,1}.data(:,1:5);
    Pos.pos.direction=headdirection;
end




end


%cd(animdirect);
%eval([varprefix,'rawpos = rawpos']);
%eval(['save ',fileprefix,'rawpos',daystring,' ', varprefix,'rawpos']);
%eval([varprefix,'pos = pos']);
%eval(['save ',fileprefix,'pos',daystring,' ',varprefix,'pos']);


