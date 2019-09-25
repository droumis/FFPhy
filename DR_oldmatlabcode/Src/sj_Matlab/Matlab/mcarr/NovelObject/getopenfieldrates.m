function output = getopenfieldrates(index, excludetimes, spikes,pos, task, varargin)
%
%For each time determines:
    % x, y position
    % valid spiketimes
    % quad position

for option = 1:2:length(varargin)-1
    if isstr(varargin{option})
        switch(varargin{option})
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end
    else
        error('Options must be strings, followed by the variable');
    end
end

%Define variables and filter out exclude times
pos = pos{index(1,1)}{index(1,2)}.data;
task = task{index(1,1)}{index(1,2)};
valid =  ~isExcluded(pos(:,1), excludetimes);

quadtypes = task.quadrants;
pos = pos(valid,:);

type = zeros(4,1);
%Determine the type of each quadrant
for q = 1:4
    %Determine what type of quadrant q is
    if task.objects(q) && task.novel(q);
        type(q) = 1;
        quadtypes(task.quadrants==q) = 1;
    elseif task.objects(q) && ~task.novel(q);
        type(q) = 2;
        quadtypes(task.quadrants==q) = 2;
    elseif ~task.objects(q) && task.novel(q);
        type(q) = 3;
        quadtypes(task.quadrants==q) = 3;
    elseif ~task.objects(q) && ~task.novel(q);
        type(q) = 4;
        quadtypes(task.quadrants==q) = 4;
    else
        type(q) = NaN;
    end   
end

output.x = pos(:,2);
output.y = pos(:,3);
output.time = pos(:,1);
output.quad = quadtypes(valid);
output.types = [task.objects; task.novel];
output.quadcenters = task.linearcoords;

%Go through each cell
output.spikes = cell(size(index,1),1);
output.index = index;
for cells = 1:size(index,1)
        
    s = spikes{index(cells,1)}{index(cells,2)}{index(cells,3)}{index(cells,4)}.data;
    
	indgoodspikes = ~isExcluded(s(:,1), excludetimes);
    spiketimes = s(indgoodspikes,1);
    output.spikes{cells} =spiketimes;
end

end