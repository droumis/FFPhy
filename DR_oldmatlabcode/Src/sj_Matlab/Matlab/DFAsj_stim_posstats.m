
function [out] = DFAsj_stim_posstats (index,excludeperiods,prefix,DIO,pos,linpos,varargin)

% Shantanu - Apr 17 2012. From DFAsj_stim_behstats. Just keeping position and velocity here - for David, etc.
% To be used with DFAsj_getstimpos_summ 

% From sj_stim_behstats. No plotting. Put all output for [day epoch] in out.
% Get Position and Velocity of Stimulations in run epochs for given animal - Borrow from sj_stim_behstats_onlyvelocity
% eg sj_stim_behstats('/data25/sjadhav/RippleInterruption/RCb_direct','RCb',1:8,[2 4],0,0,0);


%% Fixed parameters
Fs=29.97; %video sampling rate
tsamp=1/Fs; % in sec

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
velindk = posind;

%velstim_befaft=[];
for i=1:length(velindk)
    vel_stim(i) = smooth_currvel(velindk(i));
end

% --------------------
% Save 
% --------------------

% Stim-time and rate. Can also get stim-isi from this
out.stimtime = pt; % stimtimes in current day and epoch, in secs
out.totalstim = length(pt);
out.totaltime = totaltime; % in sec
% Vel and position at stimulation
out.velstim = vel_stim;
out.posstim = pos_stim;
% Vel and pos in entire epoch
out.vel = smooth_currvel;
out.pos = currpos;




        
        
        
        
        


     
     
       
    
  


