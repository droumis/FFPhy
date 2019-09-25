
function [out] = DFAsj_behstats (index,excludeperiods,prefix,pos,linpos,varargin)

% From DFAsj_stim_behstats and sj_stim_behstats. No plotting. Put all output for [day epoch] in out.
% Get Position and Velocity of Stimulations in run epochs for given animal - Borrow from sj_stim_behstats_onlyvelocity
% eg sj_stim_behstats('/data25/sjadhav/RippleInterruption/RCb_direct','RCb',1:8,[2 4],0,0,0);


%% Fixed parameters
Fs=29.97; %video sampling rate
tsamp=1/Fs; % in sec
lockout=0.25; % in sec
dur = 15*60; % Standard epoch duration in secs
% Use Calebs speed filter if necessary
defaultfilter = 'velocitydayprocess_filter.mat';
eval(['load ', defaultfilter]);
L = length(velocityfilter.kernel);
fix=0;

% Variable parameters - You can do these calculations in the script
thrsvel=5;  %<5cm per sec is still on track (3/2 for sleep session in sj_stimresp2_withvel.m)
pastsecs = 1; % For looking at history of Beh, eg. velocity: used in stimresp2
thrsdistwell = 10; %in cm, for defining well boundary
thrsdistint = 10; %in cm, for defining intersection boundary


for option = 1:2:length(varargin)-1
    switch varargin{option}
        case 'thrsvel'
            thrsvel = varargin{option+1};
        case 'pastsecs'
            proptetrodes = varargin{option+1};
        case 'thrsdistwell'
            appendindex = varargin{option+1};
        case 'thrsdistint'
            minenergy = varargin{option+1};
        case 'armthrsdist'
            minthresh = varargin{option+1};
        otherwise
            error(['Option ''', varargin{option}, ''' not defined']);
    end
end

if isempty(excludeperiods), % excludeperiods are empty
    excludeperiods=[];
end

% Current day and epoch
% ---------------------
day = unique(index(:,1));
epoch = unique(index(:,2));


% Get Position and Velocity
% ----------------------------
currtime = pos{day}{epoch}.data(:,1);
totaltime = length(currtime)/Fs; % in sec. This will be inaccurate-limited by Fs. ~100ms
if size(pos{day}{epoch}.data,2)>5 % already smoothed position and filtered velocity
    currpos = pos{day}{epoch}.data(:,6:7);
    currvel = pos{day}{epoch}.data(:,9);
    smooth_currvel=currvel;
else
    currpos = pos{day}{epoch}.data(:,2:3);
    currvel = pos{day}{epoch}.data(:,5);
    % Look at instantaneous smoothed velocity - Smooth velocity
    smooth_currvel = [filtfilt(velocityfilter.kernel, 1, currvel)];
    %smooth_currvel = smoothvect(currvel,filt); 2nd option for convolution
end

% Vel-Moving vs Still
% --------------------
% Divide all velocity into moving and still
moving_vel =  smooth_currvel(find( smooth_currvel>thrsvel));
still_vel =  smooth_currvel(find( smooth_currvel<=thrsvel));
time_moving = length(moving_vel)/Fs;
time_still = length(still_vel)/Fs;

% Get Well and Intersection positions
wellpos = linpos{day}{epoch}.wellSegmentInfo.wellCoord;
intpos(1,:) = linpos{day}{epoch}.segmentInfo.segmentCoords(1,3:4);
intpos(2,:) = linpos{day}{epoch}.segmentInfo.segmentCoords(2,3:4);
intpos(3,:) = linpos{day}{epoch}.segmentInfo.segmentCoords(4,3:4);


% After review - Get Time Spent at Wells
% --------------------------------------

% For all 3 wells - get no. of time posn is within thrsdist to well and convert to time
for i=1:3
    currwellpos = wellpos(i,:);
    nposwell(i) = length(find(dist(currpos,repmat(currwellpos,length(currpos),1)) <= thrsdistwell));
    timewell(i) = nposwell(i)*tsamp;
end

% Sum time across all 3 wells
totaltimewells = sum(timewell); 

% Now, you also need number of trials. Use yout behavior plotting function for this
[outbound_logic] = sj_day_findoutbound (linpos,day,epoch);
[inbound_logic] = sj_day_findinbound (linpos,day,epoch);

Ntrials = length(outbound_logic) + length(inbound_logic);
Ntrials_out = length(outbound_logic);
Ntrials_in = length(inbound_logic);
avgwelltime = totaltimewells/Ntrials;
avgwelltime1 = timewell(1)/Ntrials_in; %Inbounds end at vtr well
avgwelltime2 = timewell(2)/(0.5*Ntrials_out); %Outbounds end at side wells
avgwelltime3 = timewell(3)/(0.5*Ntrials_out);

% --------------------
% Save 
% --------------------

% Vel and pos in entire epoch
out.vel = smooth_currvel;
out.pos = currpos;
% Well and intersection positions for current epoch - if you need to calculate pos_stim
out.wellpos = wellpos;
out.intpos = intpos;
out.time_moving = time_moving; % in sec
out.time_still = time_still; % in sec
out.vel_moving = moving_vel;
out.vel_still = still_vel;
%Time at reward wells
out.totaltimewells = totaltimewells;
out.nposwell = nposwell;
out.timewell = timewell;
out.Ntrials = Ntrials;
out.avgwelltime = avgwelltime;
out.avgwelltime1 = avgwelltime1;
out.avgwelltime2 = avgwelltime2;
out.avgwelltime3 = avgwelltime3;
        
        
        
        
        
        


     
     
       
    
  


