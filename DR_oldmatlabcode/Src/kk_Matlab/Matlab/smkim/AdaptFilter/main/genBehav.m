function data= genBehav(info)

% generate behavioral data, i.e. velocity, position and theta phase, given
% parameters in info
%
% info.
%    startTime,endTime      : [sec]
%    dT                     : size of one timestep [sec]
%    minX, maxX             : [cm]
%    velocity               : [cm/sec]
%    thetaFreq              : frequency of theta oscillation [1/sec]


% in future we need some randomness in trajetory
%    meanVel, sigmaVel      : velocities ~ N(meanVel, sigmaVel^2) [cm/sec]

randn('state',sum(100*clock));

ts= info.startTime;
te= info.endTime;
dT= info.dT;

time= [ts:dT:te-dT]';
TN= length(time);

% generate velocities
%vel= info.meanVel + info.sigmaVel*randn(TN,1);
vel= info.velocity* ones(TN,1);

% integrate velocities to get position information
dL= info.maxX-info.minX;
pos= mod(cumsum(vel * dT), dL);
pos= pos + info.minX*ones(TN,1);
 
% generate theta phase
phase= mod(2*pi*info.thetaFreq*time, 2*pi);

% @@ tmp.
id=zeros(TN,1);

data.time= time;
data.linpos= pos;
data.phase= phase;
data.traj= id;
