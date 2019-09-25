function out = calcvelocitysegment(index, excludetimes, pos, varargin)
% function out = calclinvelocitysegment(index, excludetimes, linpos, varargin)
%
% Computes the mean 2D velocity for each valid time segment. Uses the
% 5th column of the pos structure for the estimate of 2D velocity.
%
% Options:
%   minO: defines a minimum segment length, default is no minimum
%   maxO: defines a maximum segment length, default is no maximum
%

% set options
smooth = 0;
minO = [];
maxO = [];

for option = 1:2:length(varargin)-1   
    if isstr(varargin{option})       
        switch(varargin{option})
            case 'minO'
                minO = varargin{option+1};
            case 'maxO'
                maxO = varargin{option+1};
            case 'smooth'
                smooth = varargin{option+1};
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end        
    else
        error('Options must be strings, followed by the variable');
    end
end

time = pos{index(1)}{index(2)}.data(:,1);
if smooth == 0
    velocity = pos{index(1)}{index(2)}.data(:,5);
else
	velocity = pos{index(1)}{index(2)}.data(:,8);  
end
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
    vout(i) = mean(velocity(startind(i):endind(i)));
end

out.velocity = vout;
out.time = [starttime endtime];

end