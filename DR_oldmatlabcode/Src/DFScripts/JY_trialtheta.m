function [out] = JY_trialtheta(index, excludetimes,data,theta, varargin)

% gets the theta for each trial
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





% get the theta data for each trial


epochtheta=theta{index(1)}{index(2)};


thetastart=epochtheta{1,1}.starttime*10000;
thetasamprate=epochtheta{1,1}.samprate;





% go through each trial and get theta from all channels

trialtheta=cell(1,size(trialtimes,1));
trialthetaphase=cell(1,size(trialtimes,1));

for i=1:length(trialtimes)
    tstart=trialtimes(i,1);
    tend=trialtimes(i,2);
    
    tstartindex=ceil((tstart-thetastart)/(1*10000/thetasamprate));
    if tstartindex <1 % check for mismatch between epoch start times defined by times.mat and position reconstruction
        tstartindex=1;
    end
    tendindex=ceil((tend-thetastart)/(1*10000/thetasamprate));
    
    if tendindex>size(epochtheta{1,1}.data,1);
        tendindex=size(epochtheta{1,1}.data,1);
    end
    
    
    
    trialtheta{i}=rot90(flipud(cell2mat(cellfun(@(x) x.data(tstartindex:tendindex,1), epochtheta,'UniformOutput', false))),3);
    trialthetaphase{i}=rot90(flipud(cell2mat(cellfun(@(x) x.data(tstartindex:tendindex,2), epochtheta,'UniformOutput', false))),3);
    
end


% go through each inter trial and get theta from all channels


intertrialtheta=cell(1,size(trialtimes,1)-1);
intertrialthetaphase=cell(1,size(trialtimes,1)-1);

for i=1:length(trialtimes)-1
    tstart=trialtimes(i,2);
    tend=trialtimes(i+1,1);
    tstartindex=ceil((tstart-thetastart)/(1*10000/thetasamprate));
    tendindex=ceil((tend-thetastart)/(1*10000/thetasamprate));
    intertrialtheta{i}=rot90(flipud(cell2mat(cellfun(@(x) x.data(tstartindex:tendindex,1), epochtheta,'UniformOutput', false))),3);
    intertrialthetaphase{i}=rot90(flipud(cell2mat(cellfun(@(x) x.data(tstartindex:tendindex,2), epochtheta,'UniformOutput', false))),3);
    
end

    
%get theta phase for all times, the number of datapoints is not always the
%same between tetrodes, take the minimum
datapoints=min(cell2mat(cellfun(@(x) size(x.data,1), epochtheta,'UniformOutput', false)));
timeindex=[thetastart:(10000/thetasamprate):thetastart+(10000/thetasamprate)*(datapoints-1)];
allthetaphase=rot90(flipud(cell2mat(cellfun(@(x) x.data(1:datapoints,2), epochtheta,'UniformOutput', false))),3);

%thetastruct=loadthetastruct(out.animal{1,2},animalname,'theta',days,epochs,tetrodes);
alltimes=repmat(timeindex,size(allthetaphase,1),1);

out.trialtimes=trialtimes;
out.theta=trialtheta;
out.thetaphase=trialthetaphase;
out.intertrialtheta=intertrialtheta;
out.intertrialthetaphase=intertrialthetaphase;
out.samplerate=thetasamprate;
out.timeindex=alltimes;
out.allthetaphase=allthetaphase;


warning('ON','MATLAB:divideByZero');