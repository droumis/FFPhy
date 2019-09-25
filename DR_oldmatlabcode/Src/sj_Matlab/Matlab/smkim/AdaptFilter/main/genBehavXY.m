function data= genBehavXY(info)

% generate behavioral data, i.e. velocity, position and theta phase, given
% parameters in info
%
% info.
%    startTime,endTime      : [sec]
%    dT                     : size of one timestep [sec]
%    minX, maxX             : [a]
%    minY, maxY             : [a]
%    velX, velY             : [a/sec]


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
vx= info.velX* ones(TN,1);

% integrate velocities to get position information
dLX= info.maxX-info.minX;
pos= mod(cumsum(vx * dT), dLX);
pos= pos + info.minX*ones(TN,1);
 
vy= info.velY* ones(TN,1);
dLY= info.maxY-info.minY;
phase= mod(cumsum(vy * dT), dLY);
phase= phase + info.minY*ones(TN,1);


% @@ tmp.
id=zeros(TN,1);

data.timearray= time;
data.posarray= pos;
data.phasearray= phase;
data.fieldID= id;