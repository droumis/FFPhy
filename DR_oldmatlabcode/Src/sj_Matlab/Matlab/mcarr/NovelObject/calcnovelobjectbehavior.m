function out = calcnovelobjectbehavior(index, excludetimes, task,pos,varargin)
% function out = calcnovelobjectbehavior(index, excludetimes, task,pos,varargin)
%
%  Loads the task and pos structure for each epoch and computes
%   - percent time spent in each quadrant type for all time and for just
%   the first minute
%   - the preference for novel objects as compared to familiar objects
%   - the preference for novel objects as a function of time.
%
%   out is a structure with the following fields
%       totalpercenttime: 4 X 1 vector with % time for all minutes spent in: 
%               novel object quadrants
%               familiar object quadrants
%               novel empty quadrants
%               familiar empty quadrants
%       if there was no novel quadrants (familiar session) novel object and
%       novel empty are set to -1.

%       totalpreference: preference for novel object (or one of the
%       familiar objects if a familiar session) over the entire session

%       percenttime: 4 X 1 vector with % time for first minute spent in: 
%               novel object quadrants
%               familiar object quadrants
%               novel empty quadrants
%               familiar empty quadrants
%       if there was no novel quadrants (familiar session) novel object and
%       novel empty are set to -1.

%       preference: preference for novel object (or one of the familiar
%       objects if a familiar session) as a function of session minute

% set options
appendindex = 1;
out = [];

for option = 1:2:length(varargin)-1   
    if isstr(varargin{option})       
        switch(varargin{option})
            case 'appendindex'
                appendindex = varargin{option+1};
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end        
    else
        error('Options must be strings, followed by the variable');
    end
end

%Initialize variables
pos = pos{index(1)}{index(2)};

objects = find(task{index(1)}{index(2)}.objects);
task = task{index(1)}{index(2)};
quadrants = task.quadrants;

time = pos.data(:,1);

%Apply excludetimes
valid = logical(isExcluded(time, excludetimes));
quadrants(valid) = -1;
totaltime = sum(~valid);

%Compute percent time in each quadrant type
out.totalpercenttime = -1*ones(4,1);
out.percenttime = -1*ones(4,1);
for q = 1:4
    %Determine time spent in quadrant q for all time
    tmp = sum(quadrants == q)./totaltime;

    %Determine time spent in quadrant q for first minute
    tmp1 = sum(quadrants(1:1800) == q)./sum(~valid(1:1800));
    
    %Determine what type of quadrant q is
    if task.objects(q) && task.novel(q);
        type = 1;
    elseif task.objects(q) && ~task.novel(q);
        type = 2;
    elseif ~task.objects(q) && task.novel(q);
        type = 3;
    elseif ~task.objects(q) && ~task.novel(q);
        type = 4;
    else
        type = NaN;
    end
    
    if out.percenttime(type) == -1
        out.totalpercenttime(type) = tmp;
        out.percenttime(type) = tmp1;
        clear tmp tmp1
    else
        
        out.totalpercenttime(type) = out.totalpercenttime(type) + tmp;
        out.percenttime(type) = out.percenttime(type) + tmp1;
        clear tmp tmp1
    end
end

%Compute preference score for all time and for each 30 seconds
obj1 = quadrants == objects(1);
obj2 = quadrants == objects(2);
if any(find(task.novel) == objects(2))
    obj1 = quadrants == objects(2);
    obj2 = quadrants == objects(1);
end
out.totalpreference = (sum(obj1)./totaltime-sum(obj2)./totaltime)./(sum(obj1)./totaltime+sum(obj2)./totaltime);

%Compute preference score over each minute
t = lookup(time,time(1):60:time(end),-1);

%Get rid of any bins that are shorter than 15 seconds
invalid = find(hist(t,1:max(t))<450);
if ~isempty(invalid)
    for i = 1:length(invalid)
        t(t==invalid(i)) = -1;
    end
end

out.preference = zeros(max(t),1);
for q = 1:max(t)
    out.preference(q) = (sum(obj1(t<=q))./sum(~valid(t<=q))-sum(obj2(t<=q))./sum(~valid(t<=q)))./(sum(obj1(t<=q))./sum(~valid(t<=q))+sum(obj2(t<=q))./sum(~valid(t<=q)));
end

end