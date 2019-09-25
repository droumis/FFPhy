function velocitydayprocess(directoryname,fileprefix,days,varargin)
% function velocitydayprocess(directoryname,fileprefix,days, varargin)
%
%DIRECTORYNAME -- example '/data99/user/animaldatafolder/', a folder 
%                 containing processed matlab data for the animal
%
%FILEPREFIX --    animal specific prefix for each datafile (e.g. 'fre')
%
%DAYS --          a vector of experiment day numbers 
%
%OPTIONS:   CMPERPIX-- the size of each pixel in centimeters (default 1)
%
%
cmperpix = 1;

[otherOptions] = procOptions(varargin);

for day = days
   fname = fullfile(directoryname,sprintf('%srawpos%02d.mat',fileprefix,day));
   load(fname);
   fname = fullfile(directoryname,sprintf('%spos%02d.mat',fileprefix,day));
   load(fname);
   fprintf('(%d)',day);

   nepochs = length(rawpos{day});
   
   cmperpix = pos{day}{2}.cmperpixel;
   tmppos = estimate_position(rawpos{day},'centimeters_per_pixel',cmperpix);
   for e = 1:nepochs
       if ~isempty(rawpos{day}{e})
           position = nan(size(rawpos{day}{e}.data,1),9);
           speedind = lookup(pos{day}{e}.data(:,1),tmppos{e}.data(:,1));
           %keep non-smoothed x and y position and original estimates of
           %head position and velocity
           position(speedind,2:5) = pos{day}{e}.data(:,2:5);
           %add denoised estimate of x-position, y-position, and head
           %position from estimate_position.m
           position(:,[1 6 7 9]) = tmppos{e}.data(:,[1 2 3 4]);
           %add smoothed velocity based on denoised estimate of x and y
           %position using add_velocity function below
           position(:,8) = add_velocity(position(:,6:7),position(:,1));
           pos{day}{e}.data = position;
           pos{day}{e}.descript{3} = sprintf('estimated velocity and head direction using esitmate_position.m');
           while length(pos{day}{e}.descript)>3
               pos{day}{e}.descript = pos{day}{e}.descript(1:end-1);
           end
           pos{day}{e}.fields = 'time x y dir vel smooth-x smooth-y smooth-v smooth-dir';
           fprintf('.');
       end
   end
   % save the resulting file
   save(fname, 'pos');
   clear pos rawpos
end
fprintf('\n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ESTIMATE VELOCITY
function vel = add_velocity(pos, time, filt)
    % assign default values to varargin options
    if (nargin == 2)
        filt = gaussSmooth(30,0.5/sqrt(2),2);
    end

    if size(pos,2) ~= 2
         error(['input vector must be N by 2']);
    end
    if size(pos,1) ~= size(time)
        error(['time and pos vectors must be same length']);
    end
    frametime = diff(time);

    % get all invalid positions 
    indzeropos = find((pos(:,2) == 0));

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
                goodpos = [goodpos [(i+1); size(pos,1)]];
            end
            % now, save the current index to invalid pos as the previous index to
            % invalid pos
            lastzeroind = i;
        end
    else
        goodpos = [1; size(pos,1)];
    end

    filtlength = length(filt);
    % now iterate through every 'good' chunk of pos and smooth it
    for i = goodpos
        if ((i(2)-i(1)) > 3*filtlength)
            pos(i(1):i(2),:) = [filtfilt(filt, 1, pos(i(1):i(2),1)) ...
                                filtfilt(filt, 1, pos(i(1):i(2),2))];
        end
    end

    % compute the velocity as the distance between adjoining points divided by the
    % frame time
    vel = zeros(length(time), 1);
    vel(2:end) = dist(pos(1:end-1,:), pos(2:end,:)) ./ frametime;
    vel(1) = vel(2);

    % zero out the velocities that correpond to invalid pos's.  zero out the
    % velocity at both the index that contains the invalid pos and the
    % subsequent index (although it may contain a valid pos)

    indzerovel = unique([indzeropos; (indzeropos+1)]);
    vel(indzerovel) = 0;
end % end function add_velocity