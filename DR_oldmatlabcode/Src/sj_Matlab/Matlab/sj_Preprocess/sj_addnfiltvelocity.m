
function [vel_sm] = sj_addnfiltvelocity(pos_sm, time, filt)
% ESTIMATE VELOCITY given x,y position. Called by sj_velocitydayprocess
% assign default values to varargin options
if (nargin < 3)
    filt = gaussSmooth(30,0.5/sqrt(2),2); %Maggie's filter for position. 
    
    % Option: Filtering velocity with real-time filter
    filt = gaussian(15,30); % Non-causal part of this filter is used in NSpike2 for online speed estimation
                            % Npts = 30, and sigma is half length of filter.  
                            % 15 pts at 29.97 fps is 15*33.33 = 500ms. So
                            % filter is 0.5 second on each side
    filt_rt = 2*filt(1:15);   
end

% Use Calebs filter
defaultfilter = 'velocitydayprocess_filter.mat';
eval(['load ', defaultfilter]);
L = length(velocityfilter.kernel);

if size(pos_sm,2) ~= 2
    error(['input vector must be N by 2']);
end
if size(pos_sm,1) ~= size(time)
    error(['time and pos vectors must be same length']);
end
frametime = diff(time);

% get all invalid positions
indzeropos = find((pos_sm(:,2) == 0));

if indzeropos
    % find all chunks of valid pos which flank contiguous regions of 0
    % one thing we know is that the first and last pos are not 0
    goodpos = [];   % this will be a 2XN array, with each column contain the
    % start and end indeces of a good chunk of pos data
    lastzeroind = 0;
    for i = indzeropos'    % cycle through all the indeces of zero pos
        % we know that we've just gotten to the end of a chunk of good pos if
        % two things happen:
        %    1) the difference between the last zero index and the current zero
        %       index is greater than 1
        %    2) or, if we've gotten to the end of the array of indeces pointing
        %       to invalid pos'sength(a)-1).
        
        % if either of these two things are true, then save the starting and
        % ending indeces of the chunk of good pos in the array 'goodpos'
        if ((i-lastzeroind)>1)
            goodpos = [goodpos [(lastzeroind+1); (i-1)]];
        end
        if (i==indzeropos(end))
            goodpos = [goodpos [(i+1); size(pos_sm,1)]];
        end
        % now, save the current index to invalid pos as the previous index to
        % invalid pos
        lastzeroind = i;
    end
else
    goodpos = [1; size(pos_sm,1)];
end

% now iterate through every 'good' chunk of position and smooth it
filtlength = length(filt);
for i = goodpos
    if ((i(2)-i(1)) > 3*filtlength)
        pos_sm(i(1):i(2),:) = [filtfilt(velocityfilter.kernel, 1, pos_sm(i(1):i(2),1)) ...
                                filtfilt(velocityfilter.kernel, 1, pos_sm(i(1):i(2),2))];
        %pos_sm(i(1):i(2),:) = [smoothvect(pos_sm(i(1):i(2),1),filt) smoothvect(pos_sm(i(1):i(2),2),filt)]; %2nd option for convolution
    end
end

% compute the velocity as the distance between adjoining points divided by the
% frame time
vel_sm = zeros(length(time), 1);
vel_sm(2:end) = dist(pos_sm(1:end-1,:), pos_sm(2:end,:)) / mean(frametime);
vel_sm(1) = vel_sm(2);


% zero out the velocities that correpond to invalid pos's.  zero out the
% velocity at both the index that contains the invalid pos and the
% subsequent index (although it may contain a valid pos)

indzerovel = unique([indzeropos; (indzeropos+1)]);
vel_sm(indzerovel) = 0;



end % end function addnfilt_velocity