
% function out = DFAsj_getplacefieldparams(index, excludetimes, spikes, linfields, mapfields, varargin)
function out = DFAsj_getplacefieldparams(index, excludetimes, linpos, spikes, varargin) %DR

% Combines calctotalmeanrate and sj_calcoverlap

appendindex = 0;
normalize=0; % redundant - used in sj_calcoverlap
binsize=2; % for already linearized trajectories: in DFAsj_filtercalclinfields_tf
thresh = 3; % threshold peak rate used in sj_calcoverlap 
thresh1 =1; % 1 Hz firing rate for area. Also use 3 Hz, and return both
minbins = 0.5;

if ~isempty(excludetimes)
    excludetimes = [];
end

for option = 1:2:length(varargin)-1
    if isstr(varargin{option})
        switch(varargin{option})
            case 'binsize'
                binsize = varargin{option+1};
            case 'thresh'
                thresh = varargin{option+1};
            case 'minbins'
                minbins = varargin{option+1};
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end
    end
end

% From calctotalmeanrate
% -----------------------

ontime = diff(spikes{index(1)}{index(2)}{index(3)}{index(4)}.timerange)/10000;
if ~isempty(excludetimes)
    totalexclude = sum(excludetimes(:,2) - excludetimes(:,1));
else
    totalexclude = 0;
end

totalontime = ontime-totalexclude;
if ~isempty(spikes{index(1)}{index(2)}{index(3)}{index(4)})
    if ~isempty(spikes{index(1)}{index(2)}{index(3)}{index(4)}.data)
        goodspikes = ~isExcluded(spikes{index(1)}{index(2)}{index(3)}{index(4)}.data(:,1), excludetimes);
    else
        goodspikes = 0;
    end
else
    goodspikes = nan;
end
numgoodspikes = sum(goodspikes);
meanrate = numgoodspikes/totalontime;

% From sj_calcoverlap - dealing with only 1 cell here. simpler
% ---------------------------------------------------------

% calc linearized fields -DR
% ------------------------
statematrix = linpos{index(1)}{index(2)}.statematrix;
lindist = statematrix.lindist;
state = statematrix.traj;
exclind = isExcluded(linpos{index(1,1)}{index(1,2)}.statematrix.time,excludetimes);  %find excluded time indeces
state(find(exclind)) = -1; %set excluded times to -1
trajdata = calclinfields(spikes,state,lindist,linpos, index); %DR
% trajdata = linfields{index(1)}{index(2)}{index(3)}{index(4)}; 

% [overlap, peakcomb, trajpeakcomb, peak1, peak2, trajpeak1, trajpeak2] = sj_calcoverlap(lf1,lf2,...
%     'normalize',normalize,'thresh',thresh,'minbins',minbins);
peaks=[];
for traj = 1:length(trajdata)
        if ~isempty(trajdata{traj}) 
            peaks = [peaks (max(trajdata{traj}(:,5)))];
        end
end
[peakrate,peaktraj] = max(peaks);

% Get area under field
% --------------------
for traj = 1:length(trajdata)
   
    if ~isempty(trajdata{traj})
        trajlth(traj) = length(trajdata{traj}(:,5));
        currtrajrate = trajdata{traj}(:,5);
        % Find number of bins with rate > thresh1(1Hz) and thresh(3Hz)
        lthunder1(traj)=length(find(currtrajrate>thresh1));
        lthunder3(traj)=length(find(currtrajrate>thresh));
        %DR mean thresh within traj for each cell
        threshmean = mean(currtrajrate);
        lthundermean(traj)=length(find(currtrajrate>threshmean));
         end
end

% Get fraction for trajs separately, and also combine for trajs
fracunder1=lthunder1./trajlth; % for each traj
fracunder3=lthunder3./trajlth; % for each traj
fracundermean=lthundermean./trajlth; % for each traj
% For all trajs combined    
total_lthunder1=sum(lthunder1);
total_lthunder3=sum(lthunder3);
total_lthundermean=sum(lthundermean);
total_trajlth=sum(trajlth);
total_fracunder1=total_lthunder1/total_trajlth; % Can also take mean of fracunder1
total_fracunder3=total_lthunder3/total_trajlth;
total_fracundermean=total_lthundermean/total_trajlth;


% Put in structure and return
% -----------------------------

out.index = index;
% Meanrate
out.meanrate=meanrate; % for epoch
% Peakrate
out.peakrate=peakrate;
% Which is the trajectory with the max firing rate. Can be multiple
out.peaktraj=peaktraj;
% For all trajectories combined - scalar no
out.total_fracunder1=total_fracunder1;
out.total_fracunder3=total_fracunder3;
out.total_fracundermean=total_fracundermean;
out.total_trajlth=total_trajlth;
out.total_lthunder1=total_lthunder1;
out.total_lthunder3=total_lthunder3;
out.total_lthundermean=total_lthundermean;
% For individual trajectories separately - Vector with length ntraj, usually 4
out.fracunder1=fracunder1;
out.fracunder3=fracunder3;
out.fracundermean=fracundermean;
out.trajlth=trajlth;
out.lthunder1=lthunder1;
out.lthunder3=lthunder1;
% Return linearized trajectories and 2d maps as well
out.trajdata = trajdata;
% out.mapdata = mapfields{index(1)}{index(2)}{index(3)}{index(4)}; %DR 






