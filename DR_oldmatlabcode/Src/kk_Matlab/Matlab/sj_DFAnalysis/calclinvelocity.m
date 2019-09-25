function out = calclinvelocity(index, excludetimes, linpos, varargin)
% out = calclinvelocity(linpos,index, varargin)
% Produces a cell structure with the fields:
% time, velocity (linearized)
%   INDEX - N by 2 vector [day epoch]
%
%   OPTION: 'smooth', default no smoothing
%                   to compute linear speed, can smooth linear position data.
%           'smoothwidth', default 2
%                   It is smoothed with a gaussian of length VSW and std VSW/4.
%                   default for lineardayprocess is 2 seconds


smooth = [];
smoothwidth = 2;
excludetimes = [];

for option = 1:2:length(varargin)-1   
    if isstr(varargin{option})       
        switch(varargin{option})
            case 'smooth'
                smooth = varargin{option+1};
            case 'smoothwidth'
                smoothwidth = varargin{option+1};
            case 'excludetimes'
                excludetimes = varargin{option+1};
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end        
    else
        error('Options must be strings, followed by the variable');
    end
end

if isempty(excludetimes)
    out.time = linpos{index(1)}{index(2)}.statematrix.time;
    goodtimes = 1:length(out.time);
else
    goodtimes = ~isExcluded(linpos{index(1)}{index(2)}.statematrix.time,excludetimes);
    out.time = linpos{index(1)}{index(2)}.statematrix.time(goodtimes);
end

timestep = mean(diff(linpos{index(1)}{index(2)}.statematrix.time));
if isempty(smooth)
    linv = abs(linpos{index(1)}{index(2)}.statematrix.linearVelocity(:,1));
    linv=linv(goodtimes);
elseif ~isempty(smooth)
    %to smooth
    npoints = smoothwidth/timestep;
    filtstd = smoothwidth/(4*timestep);
    % the default filter for smoothing motion direction is a n second long gaussian
    % with a n/4 second stdev
    kernal = gaussian(filtstd, npoints);
    % smooth the linear distances with the filter and then go through the linear
    % distance positions and take all of the differences to get velocities
    smoothdist = smoothvect(linpos{index(1)}{index(2)}.statematrix.lindist, kernal);
    linv = diff([smoothdist; 0; 0]) / timestep;
    linv = linv(goodtimes);
end
    out.velocity=abs(linv);

end