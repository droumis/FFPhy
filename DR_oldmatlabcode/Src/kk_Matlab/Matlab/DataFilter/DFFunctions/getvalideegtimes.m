function out = getvalideegtimes(index,phase,time,varargin)

%% REMOVED 'INDEX' argument on 3.13.12 % KK

%This function looks through a trace of phases and determines which times
%are valid based on the desired frequency band. Default is to work as a
%theta filter.

%Set Defaults
mintime = 0.1;
maxtime = 0.15;

for option = 1:2:length(varargin)-1   
    if isstr(varargin{option})       
        switch(varargin{option})
            case 'min'
                mintime = varargin{option+1};
            case 'max'
                maxtime = varargin{option+1};
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end        
    else
        error('Options must be strings, followed by the variable');
    end
end

%Detect local troughs
[max min] = peakdet(phase,0.1,time);

%Determine if troughs are between minHz and maxHz
valid = diff([min(:,1); 0])>mintime & diff([min(:,1); 0])<maxtime & min(:,2)<-3;

%Construct a time vector of valid and invalid times
goodtimes = zeros(length(time),1);
for i=1:length(valid)-1
    if valid(i)
        starttime = lookup(min(i,1),time);
        endtime = lookup(min(i+1,1),time)-1;
        goodtimes(starttime:endtime)=1;
    end
end

out = goodtimes;
end