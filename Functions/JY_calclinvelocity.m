function out = JY_calclinvelocity(linpos,index, varargin)
% out = calclinvelocity(linpos,index, varargin)
% Produces a cell structure with the fields:
% time, velocity (linearized)
%   INDEX - N by 2 vector [day epoch]
%
%   OPTION: 'smooth', default no smoothing
%                   to compute linear speed, can smooth linear position data.
%                   It is smoothed with a gaussian of length VSW and std VSW/4.
%                   default for lineardayprocess is 2 seconds
% 20111216: change linearvelocity to reference in
%           Data{day}{epoch}Pos.correcteddata column 5
%           smoothing is not supported since lindist is undefined
% 20130304: get velocity data from Pos.corrected data, it's a better
% version of the smoothed velocity data

smooth = [];
if ~isempty(varargin)
    smooth = varargin{2};
end

%out.time = linpos{index(1)}{index(2)}.statematrix.time; 
out.time = linpos{index(1)}{index(2)}.Pos.correcteddata(:,1); 

%timestep = mean(diff(linpos{index(1)}{index(2)}.statematrix.time));
timestep = mean(diff(out.time));

if isempty(smooth)
    %linv=abs(linpos{index(1)}{index(2)}.statematrix.linearVelocity);
    linv=abs(linpos{index(1)}{index(2)}.Pos.correcteddata(:,5));
    %linv(linpos{index(1)}{index(2)}.statematrix.linearVelocity==-1)=NaN;
elseif ~isempty(smooth)
    %to smooth
    npoints = smoothwidth/timestep;
    filtstd = smoothwidth/(4*timestep);
    % the default filter for smoothing motion direction is a n second long gaussian
    % with a n/4 second stdev
    filt = gaussian(filtstd, npoints);
    % smooth the linear distances with the filter and then go through the linear
    % distance positions and take all of the differences to get velocities
    %smoothdist = smoothvect(linpos{index(1)}{index(2)}.statematrix.lindist, filt);
    smoothdist = smoothvect(linpos{index(1)}{index(2)}.Pos.correcteddata(:,5), filt);
    
    linv = diff(smoothdist) / timestep;
end
    out.velocity=linv;

end