function [out] = dfakk_getrewardtrigspiking(index, excludeperiods, spikes, rewardinfo, ripples, pos, varargin)
% [out] = dfakk_getrewardtrigspiking(index, excludeperiods, spikes, rewardinfo, varargin)
%
%   This function transcribes reward triggered spiking for every reward trigger. The ripples
%   are defined by the first group of cells and the spikes are from the
%   second group of cells.
%
%   use singlecelleeganal
%
%   index [day epoch tetrode cell]
%
%   options are 
%   'window', 1x2 vector specifies the window before and after each included ripple.
%                   Default is 100 mseconds before and 15 seconds after
%                   ripple start time.
%   out = out.out   An R x C sized matrix where each entry is the number of
%                   spikes cell C fired during ripple R
%         out.index [D E T C], gives the identity of the cells for each
%                   column in out.out
%         out.times [starttime endtime], givest the starttime and endtime
%                   of the ripples for each row in out.out

% default options
posflag = 0;
pseudodelay = 0;
minthresh = 0;
window = [5 5];  % in sec
binsize = 0.001;  % in sec

for option = 1:2:length(varargin)-1
    switch varargin{option}
        case 'posflag'
            posflag = varargin{option+1};
        case 'window'
            window = varargin{option+1};
        case 'binsize'
            binsize = varargin{option+1};
        case 'minthresh'                          % ripple SD minimum threshold
            minthresh = varargin{option+1};    
        case 'pseudo'
            pseudodelay = varargin{option+1};     % impose offset (in sec) for error trials, to mimic delay in error. e.g. Egypt had 500 ms delay in reward 
        otherwise
            error(['Option ''', varargin{option}, ''' not defined']);
    end
end

% establish times vector (suited for histc, where center bin straddles time 0)
out.index = index;
times = (-window(1)-0.5*binsize):binsize:(window(2)+0.5*binsize);
    timescenters = times(1:(end-1))+binsize/2;      
out.times = times;

reward = rewardinfo{index(1)}{index(2)};
spiketimes = spikes{index(1)}{index(2)}{index(3)}{index(4)}.data(:,1);


%collect error and reward triggers
        % first gather unique reward triggers (sweep through all epochs)
    wells = [];
            for d=1:20
                for e=1:20
                    try
                        wells = [wells ; rewardinfo{d}{e}(:,1)];
                    catch
                    end
                end
            end
     wells = unique(wells(~isnan(wells)));

     % construct trigmatrix --  error (0) and reward (1) triggers and inbound (10)-outbound (11)
     % also initialize output
     trigmatrix = [];
     for w=1:length(wells)
         trigmatrix = [trigmatrix ; wells(w) 0 10];  % well error inbound
         trigmatrix = [trigmatrix ; wells(w) 1 10];  % well reward inbound
         trigmatrix = [trigmatrix ; wells(w) 0 11];  % well error outbound
         trigmatrix = [trigmatrix ; wells(w) 1 11];  % well reward outbound
     end
     
     % initialize output matrix     
     out.trigmatrix = trigmatrix;
     out.trigmatrixdescript = 'well #  //  error (0) or reward (1) // inbound (10) or outbound (11)';
     out.output = {};
     out.ripples = {};
     out.rewardtimes = {};
     out.rewardoutcomes = {};
     out.velocity = {};
     for t=1:size(trigmatrix,1)
         out.output{t} = [];
         out.ripples{t} = [];
         out.rewardtimes{t} = [];
         out.rewardoutcomes{t} = [];
         out.velocity{t} = [];
     end
     
    % for this epoch, construct array where each element represents whether
    % ripple is present on cell's parent tetrode (1 ms timesteps).
	r = ripples{index(1)}{index(2)}{index(3)};
	timevec = r.timerange(1):0.001:r.timerange(end);
	nrip = zeros(size(timevec));
	    % apply the minthresh threhsold
	    rvalid = find(r.maxthresh > minthresh);
	    rtimes = [r.starttime(rvalid) r.endtime(rvalid)];
	    % create another parallel vector with bordering times for zeros
	    nrtimes = [(rtimes(:,1) - 0.00001) (rtimes(:,2) + 0.00001)];
	    rtimes = reshape(rtimes', length(rtimes(:)), 1); 
	    rtimes(:,2) = 1;
	    nrtimes = [r.timerange(1) ; reshape(nrtimes', ...
	    	length(nrtimes(:)), 1) ; r.timerange(2)];
	    nrtimes(:,2) = 0;
	    % create a new list with all of the times in it
	    tlist = sortrows([rtimes ; nrtimes]);
	    % use interp to create a set of ones and zeros for each time
	    % and add to nrip to get a cumulative count of the number of
	    % ripples per timestep
		nrip = interp1(tlist(:,1), tlist(:,2), timevec, 'nearest');
        ripstamps = nrip .* timevec;

     
     % collect spikes in window around reward times
     % iterate over each reward trigger, iterating by trigger type
     for t=1:size(trigmatrix,1)
         
         % if pseudo specified, then check if error trial
         %   - if so, institute the known offset (e.g. 0.5 sec)
         if pseudodelay ~= 0
             if trigmatrix(t,2)==0
                 offset = pseudodelay;
             elseif trigmatrix(t,2)==1
                 offset = 0;
             end            
         end
         
         row = rowfind(trigmatrix(t,:),reward(:,[1 3 4]));
         while row ~= 0
             % collect spike times
             spikehist = histc(spiketimes,reward(row,2)/10000 + offset + times);
             out.output{t} = [out.output{t} ; spikehist'];
             % collect velocity corresponding to the window
             if posflag
             vel = nan(size(timescenters));
             bincenters = reward(row,2)/10000 + offset + timescenters;
             for ii=1:length(timescenters)
                posind = lookup(bincenters(ii),pos{index(1)}{index(2)}.data(:,1));
                vel(ii)=pos{index(1)}{index(2)}.data(posind,9);
             end
             out.velocity{t} = [out.velocity{t} ; vel];
             disp('done w/ pos!')
             end
             % collect ripple times
             riphist = histc(ripstamps,reward(row,2)/10000 + offset + times);
             out.ripples{t} = [out.ripples{t} ; riphist];             
             % transcribe reward time in .rewardtimes
             out.rewardtimes{t} = [out.rewardtimes{t} ; reward(row,2) + offset];
             % transcribe reward outcome in .rewardoutcomes
             out.rewardoutcomes{t} = [out.rewardoutcomes{t} ; reward(row,3)];
             % strike off trigger
             reward(row,[1 3 4])=[NaN NaN NaN];
             % look for another
             row = rowfind(trigmatrix(t,:),reward(:,[1 3 4]));
         end
         
     end
         
     
     
% exclude any reward triggers that are in excludeperiods







    
end











