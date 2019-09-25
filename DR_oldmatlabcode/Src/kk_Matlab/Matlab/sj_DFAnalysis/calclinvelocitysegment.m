function out = calclinvelocitysegment(index, excludetimes, linpos, varargin)
% function out = calclinvelocitysegment(index, excludetimes, linpos, varargin)
%
%  Returns the average linear velocity for each valid time segment
%       

% set options

minO = [];
maxO = [];

for option = 1:2:length(varargin)-1   
    if isstr(varargin{option})       
        switch(varargin{option})
            case 'smooth'
                s = varargin{option+1};
            case 'smoothwidth'
                sw = varargin{option+1};
            case 'minO'
                minO = varargin{option+1};
            case 'maxO'
                maxO = varargin{option+1};
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end        
    else
        error('Options must be strings, followed by the variable');
    end
end

time = linpos{index(1)}{index(2)}.statematrix.time;

%Calculate the linear velocity across all time
linv = abs(linpos{index(1)}{index(2)}.statematrix.linearVelocity(:,2));

%apply excludetimes
includetimes = getincludedtimes(excludetimes);

%Find the start and end times for valid times
starttime = includetimes(:,1);
endtime   = includetimes(:,2);

%Only use intervals that are longer than min segment length and shorter
%than max segment length.
temp = endtime-starttime;
if ~isempty(minO)
    temp(temp<minO) = 0;
end
if ~isempty(maxO)
    temp(temp>maxO) = 0;
end
temp(temp>0) = 1;
starttime = starttime(logical(temp));
endtime = endtime(logical(temp));

%Assign the average velocity to each valid time
startind = lookup(starttime,time);
endind   = lookup(endtime,time);

vout = zeros(size(startind));
for i = 1:length(startind)
    vout(i) = mean(linv(startind(i):endind(i)));
end

out.velocity = vout;
out.time = [starttime endtime];

end