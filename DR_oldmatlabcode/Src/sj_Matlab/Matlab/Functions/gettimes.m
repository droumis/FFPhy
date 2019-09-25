function out = gettimes(timesfile)

load(timesfile)

timelist = names;
count = 0;
for i = 2:length(names)
    dash = strfind(timelist{i},'-');
    if ~isempty(dash)
        count = count+1;
        spaces = strfind(timelist{i},' ');
        lastspace = spaces(end);
        firstspace = spaces(1);
        name = timelist{i}(firstspace:lastspace);
        name(strfind(name,' ')) = '';
        name = lower(name);
        range = ranges(i,:);
        starttime = timelist{i}(lastspace+1:dash-1);
        endtime = timelist{i}(dash+1:end);
        out(count).name = name;
        out(count).range = range;
        out(count).starttime = starttime;
        out(count).endtime = endtime;        
    else
        count = count+1;
        out(count).name = [];
        out(count).range = [];
        out(count).starttime = [];
        out(count).endtime = [];
    end
end

