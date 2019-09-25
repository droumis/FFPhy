function [out] = dfakk_getrewardtrigspiking(index, excludeperiods, spikes, rewardinfo, varargin)
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

% assign the options
minthresh = 0;
window = [2 2];  % in sec
binsize = 0.001;  % in sec

for option = 1:2:length(varargin)-1
    switch varargin{option}
        case 'window'
            window = varargin{option+1};
        case 'binsize'
            binsize = varargin{option+1};
        otherwise
            error(['Option ''', varargin{option}, ''' not defined']);
    end
end

% establish times vector (suited for histc, where center bin straddles time 0)
out.index = index;
times = (-window(1)-0.5*binsize):binsize:(window(2)+0.5*binsize);      
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
     wells = unique(wells);

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
     for t=1:size(trigmatrix,1)
         out.output{t} = [];
         out.rewardtimes{t} = [];
     end
     
     % collect spikes in windows centered on trigger type
     for t=1:size(trigmatrix,1)
         row = rowfind(trigmatrix(t,:),reward(:,[1 3 4]));
         while row ~= 0
             % retrieve spike times
             spikehist = histc(spiketimes,reward(row,5)/10000+times);
             out.output{t} = [out.output{t} ; spikehist'];
             % transcribe reward time in .rewardtimes
             out.rewardtimes{t} = [out.rewardtimes{t} ; reward(row,5)];
             % strike off trigger
             reward(row,[1 3 4])=[NaN NaN NaN];
             % look for another
             row = rowfind(trigmatrix(t,:),reward(:,[1 3 4]));
         end
     end
     
     
% exclude any reward triggers that are in excludeperiods







    
end











