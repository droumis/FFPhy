function [out] = JY_trialeeg(index, excludetimes,data,meaneegspectrograms,ripple, varargin)

% gets the eeg for each trial
% get the times of reward events
% excludetime -- times want to exclude from analysis

% index - [day epoch tetrode cell]



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
                  case 'ripplereference'
                ripref = varargin{option+1};   
                
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end
    else
        error('Options must be strings, followed by the variable');
    end
end
warning('OFF','MATLAB:divideByZero');

% get the start/end times of each trial

trialtimes=data{index(1)}{index(2)}.Run(:,3:4);

% get the eeg data for each trial


epocheeg=meaneegspectrograms{index(1)}{index(2)};
epochripple=ripple{index(1)}{index(2)}{ripref};

eegstart=epocheeg{1,1}.starttime*10000;
eegsamprate=epocheeg{1,1}.samprate;

% go through each trial and get eeg from all channels

trialeeg=cell(1,size(trialtimes,1));


for i=1:length(trialtimes)
    % go through each tetrode
    
    tstart=trialtimes(i,1);
        tend=trialtimes(i,2);
        
        
    for jj=1:size(epocheeg,2)
        %eeg data
        eegstart=epocheeg{jj}.starttime*10000;
        
      
        
        tstartindex=ceil((tstart-eegstart)/(1*10000/eegsamprate));
        if tstartindex <1 % check for mismatch between epoch start times defined by times.mat and position reconstruction
            tstartindex=1;
        end
        tendindex=ceil((tend-eegstart)/(1*10000/eegsamprate));
        
        if tendindex>size(epocheeg{jj}.data,1);
            tendindex=size(epocheeg{jj}.data,1);
        end
        
        trialeeg{i}(jj,:)=double(epocheeg{jj}.data(tstartindex:tendindex,1));
        
        
        
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
        
        trialripple{i}=double(epochripple.data(tstartindex:tendindex,1));
    
    
   vstartindex=lookup(tstart,data{index(1)}{index(2)}.Pos.correcteddata(:,1)*10000);
   vendindex=lookup(tend,data{index(1)}{index(2)}.Pos.correcteddata(:,1)*10000);
    trialspeed{i}=data{index(1)}{index(2)}.Pos.correcteddata(vstartindex:vendindex,5);    

end

% go through each inter trial and get eeg from all channels


intertrialeeg=cell(1,size(trialtimes,1)-1);

for i=1:length(trialtimes)-1
      tstart=trialtimes(i,2);
        tend=trialtimes(i+1,1);
    for jj=1:size(epocheeg,2)
        eegstart=epocheeg{jj}.starttime*10000;

        tstartindex=ceil((tstart-eegstart)/(1*10000/eegsamprate));
        
        tendindex=ceil((tend-eegstart)/(1*10000/eegsamprate));
        
        intertrialeeg{i}(jj,:)=epocheeg{jj}.data(tstartindex:tendindex,1);
        
       
 
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
        
        intertrialripple{i}=double(epochripple.data(tstartindex:tendindex,1));
    
    
    vstartindex=lookup(tstart,data{index(1)}{index(2)}.Pos.correcteddata(:,1)*10000);
   vendindex=lookup(tend,data{index(1)}{index(2)}.Pos.correcteddata(:,1)*10000);
    intertrialspeed{i}=data{index(1)}{index(2)}.Pos.correcteddata(vstartindex:vendindex,5);   
    
end

% get original eeg
for kk=1:size(epocheeg,2)
raweeg{kk}=epocheeg{kk}.data;
end

%eegstruct=loadeegstruct(out.animal{1,2},animalname,'eeg',days,epochs,tetrodes);

out.trialtimes=trialtimes;
out.eeg=trialeeg;
out.intertrialeeg=intertrialeeg;
out.samplerate=eegsamprate;
out.raweeg=raweeg;
out.trialspeed=trialspeed;
out.intertrialspeed=intertrialspeed;
out.trialripple=trialripple;
out.intertrialripple=intertrialripple;


warning('ON','MATLAB:divideByZero');