% velocity = GETVELOCITY(posstruct, index, options)
%          returns an array of times and velocities
%               options take the form 'option', value, 'option', value, ... 
%                   'filt', filt specifies a filter to use for smoothing,
%                            (default no smoothing)

function [velocity] = getvelocity(posstruct, in, varargin)

ps = posstruct{in(1)}{in(2)};

arglist = char('filt');

% assign default values to varargin options
filt = [];

% process varargin
if (~isempty(varargin))
	assign(varargin{:});
end

% put all of the arguments in the newpos.arg field
velocity.arg = cell(size(arglist,1) * 2, 1);
for i = 1:size(arglist,1)
	velocity.arg{i*2-1} = arglist(i,:);
	velocity.arg{i*2} = eval([arglist(i,:)]);
end

% write out the command parameters to a string
descript1 = sprintf('Smoothed velocity data', inputname(1));
descript2 = ps.descript;
velocity.descript = char(descript1, descript2);
fields = 'time velocity';
velocity.fields = fields;

pos = ps.data(:,2:3);
frametime = diff(ps.data(:,1));

% get all invalid positions 
invalidpos = find((pos(:,2) == 0));

%set the invalid positions to infinity so that they can be removed later
pos(invalidpos,:) = Inf;

if (~isempty(filt))
    newpos = [smoothvect(pos(:,1), filt) smoothvect(pos(:,2), filt)];
else
    newpos = pos;
end

times = ps.data(1:(length(newpos)-1),1) + frametime / 2;
%compute the velocity as the distance between adjoining points divided by the
%frame time
velocity.data = [times dist(newpos(1:length(pos)-1,:), newpos(2:end,:)) ./ frametime];

