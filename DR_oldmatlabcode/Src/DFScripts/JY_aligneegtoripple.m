function [out] = JY_aligneegtoripple(index, excludetimes,data,eeg,meaneegspectrograms,ripple,ripples, varargin)

%
% get the times of ripples as specified by filtering
% excludetime -- times of ripples
% does spectral analysis time locked to ripples
% subtracts mean of the epoch


% index - [day epoch tetrode ]



appendindex = 0;
std = 1;
binsize = 1;
for option = 1:2:length(varargin)-1
    if isstr(varargin{option})
        switch(varargin{option})
            case 'appendindex'
                appendindex = varargin{option+1};
            case 'std'
                std = varargin{option+1};
            case 'binsize'
                binsize = varargin{option+1};
            case 'tetrodereference'
                ripref=varargin{option+1};
                
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end
    else
        error('Options must be strings, followed by the variable');
    end
end
warning('OFF','MATLAB:divideByZero');

% get times of ripples
excludeperiods=excludetimes;

%two ways to define ripple times
% use get ripples
%rippletimes = JY_getripples(index, ripples,'minstd',std,'tetrodereference',ripref)*10000;
% or manually go through ripples
% find ripples with thresholds above std
getripindex=find(ripples{index(1)}{index(2)}{ripref}.maxthresh>=std);
rippletimes=ripples{index(1)}{index(2)}{ripref}.midtime(getripindex)*10000;
ripplesize=ripples{index(1)}{index(2)}{ripref}.maxthresh(getripindex);
%out.ripplestart=rippletimes;

% sets up parameters

days=index(1);
epoch=index(2);

epocheeg=meaneegspectrograms{index(1)}{index(2)};
originaleeg=eeg{index(1)}{index(2)};
epochripple=ripple{index(1)}{index(2)}{ripref};

rippleeeg=cell(1,size(rippletimes,1));

halfwindow=5000;

for i=1:length(rippletimes)
    % go through each tetrode
    tstart=rippletimes(i,1)-halfwindow;
    tend=rippletimes(i,1)+halfwindow;
    %tend=rippletimes(i,2);+halfwindow;
    
    for jj=1:size(epocheeg,2)
        %eeg data
        eegsamprate=1500;
        eegstart=epocheeg{jj}.starttime*10000;
        originalstart=originaleeg{jj}.starttime*10000;
        tstartindex=ceil((tstart-eegstart)/(1*10000/eegsamprate));
        if tstartindex <1 % check for mismatch between epoch start times defined by times.mat and position reconstruction
            tstartindex=1;
        end
        tendindex=ceil((tend-eegstart)/(1*10000/eegsamprate));
        if tendindex>size(epocheeg{jj}.data,1);
            tendindex=size(epocheeg{jj}.data,1);
        end
        rippleeeg{i}(jj,:)=double(epocheeg{jj}.data(tstartindex:tendindex,1));
    end
      %ripple data
        ripplestart=epochripple.starttime*10000;
         
        tstartindex=ceil((tstart-ripplestart)/(1*10000/eegsamprate));
        if tstartindex <1 % check for mismatch between epoch start times defined by times.mat and position reconstruction
            tstartindex=1;
        end
        tendindex=ceil((tend-ripplestart)/(1*10000/eegsamprate));
        
        if tendindex>size(epochripple.data,1);
            tendindex=size(epochripple.data,1);
        end
        
        rippleripple{i}=double(epochripple.data(tstartindex:tendindex,3));
        
    vstartindex=lookup(tstart,data{index(1)}{index(2)}.Pos.correcteddata(:,1)*10000);
    vendindex=lookup(tend,data{index(1)}{index(2)}.Pos.correcteddata(:,1)*10000);
    ripplespeed{i}=data{index(1)}{index(2)}.Pos.correcteddata(vstartindex:vendindex,5);
    
end

%% output

out.rippleeeg=rippleeeg;
out.ripplespeed=ripplespeed;
out.rippleripple=rippleripple;
out.ripplesize=ripplesize;


warning('ON','MATLAB:divideByZero');