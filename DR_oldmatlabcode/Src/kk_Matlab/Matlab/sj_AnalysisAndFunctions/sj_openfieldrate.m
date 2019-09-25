function [output] = sj_openfieldrate(spikes,statevector,pos,index,binsize,std,minabsvel)

% Shantanu: Adding Abs Vel from "pos" field. Accurate calculation of
% velocity, and independent of minvel inputed/used in sj_getbehavestate
% Similar to sj_calclinfields 
% 10 Jun 2011

%[output] = openfieldratemap(spikes,pos, binsize, std)
%[output] = openfieldratemap(spikes,pos, binsize)
%
%Calculates the 2d occupancy normalized firing rate for the cell
%
%spikes - the spike data structure 
%pos - the position data structure
%binsize- the length of each spatial bin (default 1 cm)
%std - defines the shape of the 2d gaussian used to smooth spikerate.
%              (default 1)
%
%The output is a structure with n matrices. The matrices are: 
% 1) occupancy, 
% 2) bin vector x, 
% 3) bin vector y, 
% 4) bin spike count, 
% 5)occ normailized firing per bin, and 
% 6) smoothed occ normalized firing. 
%

if nargin < 5,
    binsize = 1; % cm square 
end

if nargin < 6,
    std = 1; 
end

if nargin < 7,
    minabsvel = 3; % cm/sec - Most conservative for runs and place fields
end

% Use Calebs speed filter if necessary
defaultfilter = 'velocitydayprocess_filter.mat';
eval(['load ', defaultfilter]);

warning('OFF','MATLAB:divideByZero');

% Get data
spikes = spikes{index(1)}{index(2)}{index(3)}{index(4)}.data;
% Position Data with Fields 'time x y dir vel x-sm y-sm dir-sm vel-sm'
posdata = pos{index(1)}{index(2)}.data;
if size(posdata,2)>5 % already smoothed position and filtered velocity
    absvel = abs(posdata(:,9)); % Should already be absolute value
else
    currvel = abs(pos{day}{epoch}.data(:,5));
    % If vel-sm does not exist, create it
    absvel = [filtfilt(velocityfilter.kernel, 1, currvel)];
end
% Using original velocity
absvel = abs(posdata(:,5));

if length(statevector)~=length(absvel)
    finalind = min([length(statevector) length(absvel)]);
    statevector = statevector(1:finalind);
    absvel = absvel(1:finalind);
end


% spikes.data fields: time x y dir not_used amplitude(highest variance
% channel) posindex x-sm y-sm dir-sm
% if size(spikes,2)>7
%     posx = 8; posy = 9;
% else
%     posx = 2; posy = 3;
% end

posx = 2; posy = 3;
posindexfield = 7;

if ~isempty(spikes)
    spikes = spikes(:,[1 posx posy posindexfield]);    
    spikes(:,5) = statevector(spikes(:,4)); %add the traj state for each spike
    spikes(:,6) = absvel(spikes(:,4)); %add the absvel for each spike
else
    spikes = [0 0 -1];
end

timestep = posdata(2,1) - posdata(1,1);

goodspikes = [];
goodspikes = spikes;
% Non sharp-wave spikes
goodspikeind = (spikes(:,5) ~= -1) & (spikes(:,6) >= minabsvel);
goodspikes = spikes(goodspikeind,:);

%make a cell array, where each cell contains data for one trajectory.
%inside each cell [binlocation occupancy spikecount firingrate]

tmpposition = (posdata(:,[2 3]));
tmpspikes = (goodspikes(:,[2 3]));

if ~isempty(tmpposition)
    minx = floor(min(tmpposition(:,1))) - 1;
    maxx = ceil(max(tmpposition(:,1))) + 1;
    binx = (minx:binsize:maxx);
    miny = floor(min(tmpposition(:,2))) - 1;
    maxy = ceil(max(tmpposition(:,2))) + 1;
    biny = (miny:binsize:maxy);
    %[output.occupancy output.xticks output.yticks] = hist2(tmpposition, binx, biny);
    % Shantanu - I am using hist2 in ~/Src/Matlab/sj_Utility
    [output.occupancy output.xticks output.yticks] = hist2(tmpposition(:,1), tmpposition(:,2), binx, biny);
    %[output.occupancy output.xticks output.yticks] = hist2d(tmpposition, min([length(binx),length(biny)]));
    
    nonzero = find(output.occupancy ~= 0);
    %[output.spikes BX BY] = hist2(tmpspikes, binx, biny);
    [output.spikes BX BY] = hist2(tmpspikes(:,1), tmpspikes(:,2), binx, biny);
    %[output.spikes BX BY] = hist2d(tmpspikes, binsize);
    output.spikerate = zeros(size(output.spikes));
    output.spikerate(nonzero) = output.spikes(nonzero) ./(timestep* output.occupancy(nonzero) );

%smoothed occupancy
    g = gaussian2(std,(6*std));
    output.smoothedspikerate = filter2(g,(output.spikerate)); % is this the right filter?
    smoothedoccupancy = [];
    smoothedoccupancy = zeros(size(output.spikes));
    smoothedoccupancy = filter2(g, output.occupancy);
    zero = find(smoothedoccupancy == 0);

    output.smoothedspikerate(zero) = -1;
end

warning('ON','MATLAB:divideByZero');
