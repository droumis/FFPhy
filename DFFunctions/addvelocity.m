% velocity = GETVELOCITY(posstruct, filt)
%          takes a pos structure and appends velocity info
%               filt (optional) -- specifies a filter to use for position 
%                           smoothing, (default no smoothing)

function [newpos] = addvelocity(posdata, filt)
%AS changed funtion call to match function name 4-30-09

% assign default values to varargin options
if (nargin == 1)
    filt = [];
end

% copy posdata into newpos
newpos = posdata;

pos = newpos.data(:,2:3);
time = newpos.data(:,1);
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

% if a filter was provided, then smooth the positions with that filter
if (~isempty(filt))
    filtlength = length(filt);
    % now iterate through every 'good' chunk of pos and smooth it
    for i = goodpos
        if ((i(2)-i(1)) > 3*filtlength)
            pos(i(1):i(2),:) = [filtfilt(filt, 1, pos(i(1):i(2),1)) ...
                                filtfilt(filt, 1, pos(i(1):i(2),2))];
        end
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

% now save the velocity data to the newpos structure.  either it will be
% appended or replace the velocity data that is already in the structure.
% first, see if 'vel' exists in 'field'
% in order to do that, find all the tokens in 'fields' and see if any of
% them match 'vel'
toks = {};
i = 1;
fields = newpos.fields;
while (fields)
    [toks{i} fields] = strtok(fields);
    i = i+1;
end

fieldcol = find(strcmp(toks, 'vel'));

if fieldcol     % this pos struct already contains 'vel', so replace the existing vel
    newpos.data(:,fieldcol) = vel;
    % then replace filt
    filtrow = find(strcmp(newpos.arg(:,1), 'filt'));
    newpos.arg{filtrow, 2} = filt;
else    % this is the first time vel has been computed for this pos structure, so add
    % put value of 'filt' in the newpos.arg field
    newpos.arg{end+1,1} = 'filt';
    newpos.arg{end,2} = filt;

    % write out the command parameters to a string
    descript1 = sprintf('Smoothed velocity data', inputname(1));
    newpos.descript{end+1,1} = descript1;
    newpos.fields = [newpos.fields ' vel'];
    
    % append the velocity data to the data structure
    newpos.data = [newpos.data vel];
end
