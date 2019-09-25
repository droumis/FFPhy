
function [out] = DFAsj_stim_behstats (index,excludeperiods,prefix,DIO,pos,linpos,varargin)

% From sj_stim_behstats. No plotting. Put all output for [day epoch] in out.
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
thrsdistwell = 12; % used to be 20 %in cm, for defining well boundary
thrsdistint = 12; %used to be 15  %in cm, for defining intersection boundary
armthrsdist = 20;

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

% Get Stimln information
% ------------------------
stim = DIO{day}{epoch}{16};
if isempty(stim),
    stim = DIO{day}{epoch}{15};
end
stim_starttime = stim.pulsetimes(:,1); %in tens of ms
stim_endtime = stim.pulsetimes(:,2);
stim_length = stim.pulselength;
stim_isi = stim.timesincelast(2:end);
% Get rid of errors
stim_isi=diff(stim_starttime);
rem=find(stim_isi<=lockout*10000);
stim_starttime(rem+1)=[];
stim_isi=diff(stim_starttime)./10000; % in sec
pt = stim.pulsetimes ./ 10000;    % in sec

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

% Stimulation indices in pos and vel vector
% ------------------------------------------
% Get indices for stimulation time in the pos vector
posind = lookup(pt(:,1), currtime);
pos_stim = currpos(posind,:);
% Skip first few and last few if needed
posind=posind(1:end);
pos_stim = currpos(posind,:);
% look at instantaneous velocity at stimulation times
velindk = posind; npre=1;npost=1;

% After revision - also add getting speed before and after stimulation in the following loop
% Get 5s before and after. So skip the first few and last stim in session

%velstim_befaft=[];
for i=1:length(velindk)
    vel_stim(i) = smooth_currvel(velindk(i));
       
    % 5 seconds Before and after stim
    if (velindk(i)>149) && (velindk(i)+150<=length(smooth_currvel))
        velstim_befaft(i,:) = smooth_currvel(velindk(i)-149:velindk(i)+150);
    end            
end


% Vel-Moving vs Still
% --------------------
% Divide into moving stimulations and still stimulations
moving_stimidx = find(vel_stim>thrsvel); moving_stimvel = vel_stim(find(vel_stim>thrsvel));
still_stimidx = find(vel_stim<=thrsvel); still_stimvel = vel_stim(find(vel_stim<=thrsvel));
% Divide all velocity into moving and still
moving_vel =  smooth_currvel(find( smooth_currvel>thrsvel));
still_vel =  smooth_currvel(find( smooth_currvel<=thrsvel));
time_moving = length(moving_vel)/Fs;
time_still = length(still_vel)/Fs;

% Position at stimulation calculations
% -------------------------------------
% Get Well and Intersection positions
wellpos = linpos{day}{epoch}.wellSegmentInfo.wellCoord;
intpos(1,:) = linpos{day}{epoch}.segmentInfo.segmentCoords(1,3:4);
intpos(2,:) = linpos{day}{epoch}.segmentInfo.segmentCoords(2,3:4);
intpos(3,:) = linpos{day}{epoch}.segmentInfo.segmentCoords(4,3:4);
% Find number of stimulations at well, at intersection and center-arm run
sstimarm=[]; pstimarm=[];
for i=1:3
    currwellpos = wellpos(i,:);
    currintpos = intpos(i,:);
    nstimwell(i) = length(find(dist(pos_stim,repmat(currwellpos,length(pos_stim),1)) <= thrsdistwell));
    nstimint(i) = length(find(dist(pos_stim,repmat(currintpos,length(pos_stim),1)) <= thrsdistint));
    
    switch i
        case 2
            stimarm = pos_stim(find(pos_stim(:,1)<currwellpos(1)+armthrsdist),:); % Find curr arm by using current well X-position
            sstimarm{i}=stimarm;
            
        case 1
            stimarm = pos_stim(find((pos_stim(:,1)<currwellpos(1)+armthrsdist) & (pos_stim(:,1)>currwellpos(1)-armthrsdist)),:);
            sstimarm{i}=stimarm;
            
        case 3
            stimarm = pos_stim(find(pos_stim(:,1)>currwellpos(1)-armthrsdist),:);
            sstimarm{i}=stimarm;
    end
    
    nstimarm(i) = length(find( (dist(stimarm,repmat(currintpos,length(stimarm),1))>thrsdistint) & (dist(stimarm,repmat(currwellpos,length(stimarm),1))>thrsdistwell) ));
    pstimarm{i} = stimarm(find( (dist(stimarm,repmat(currintpos,length(stimarm),1))>thrsdistint) & (dist(stimarm,repmat(currwellpos,length(stimarm),1))>thrsdistwell) ),:);
    
end
edgewell=[2,3]; centerwell=1;
nstimedgewell = nstimwell(edgewell);
nstimcenterwell = nstimwell(centerwell);
nstimcenterarm = nstimarm(centerwell);

% --------------------
% Save 
% --------------------

% Stim-time and rate. Can also get stim-isi from this
out.stimtime = pt; % stimtimes in current day and epoch, in secs
out.totalstim = length(pt);
out.totaltime = totaltime; % in sec
out.stimrate = length(pt)./(totaltime); % in Hz
% Vel and position at stimulation
out.velstim = vel_stim;
out.posstim = pos_stim;
% Vel and pos in entire epoch
out.vel = smooth_currvel;
out.pos = currpos;
% Well and intersection positions for current epoch - if you need to calculate pos_stim
out.wellpos = wellpos;
out.intpos = intpos;
% Old - Stuff you have calculated above about velocity
out.velstim_moving = moving_stimvel;
out.velstim_still = still_stimvel;
out.time_moving = time_moving; % in sec
out.time_still = time_still; % in sec
out.vel_moving = moving_vel;
out.vel_still = still_vel;
% Posn at stimln calculations
out.nstimwell = nstimwell;
out.nstimint = nstimint;
out.nstimarm = nstimarm;
out.pstimarm = pstimarm;
% After review
out.velstim_befaft = velstim_befaft;



        
        
        
        
        


     
     
       
    
  


