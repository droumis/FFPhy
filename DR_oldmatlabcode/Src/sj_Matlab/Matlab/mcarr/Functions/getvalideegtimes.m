function out = getvalideegtimes(phase,time,varargin)
%This function looks through a trace of phases and determines which times
%are valid based on the desired frequency band. Default is to work as a
%theta filter.

%Set Defaults
minHz = 6;
maxHz = 10;

for option = 1:2:length(varargin)-1   
    if isstr(varargin{option})       
        switch(varargin{option})
            case 'minHz'
                minHz = varargin{option+1};
            case 'maxHz'
                maxHz = varargin{option+1};
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end        
    else
        error('Options must be strings, followed by the variable');
    end
end

%Detect local troughs
[maxtab mintab] = peakdet(phase,0.1,time);

%Determine if troughs are between minHz and maxHz
valid = diff([mintab(:,1); 0])>(1/maxHz) & diff([mintab(:,1); 0])<(1/minHz) & mintab(:,2)<-2;

%Construct a time vector of valid and invalid times
goodtimes = zeros(length(time),1);
for i=1:length(valid)-1
    if valid(i)
        starttime = lookup(mintab(i,1),time);
        endtime = lookup(mintab(i+1,1),time);
        goodtimes(starttime:endtime)=1;
    end
end

out = goodtimes;
end