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
        starttime = timelist{i}(lastspace+1:dash-1);
        endtime = timelist{i}(dash+1:end);
        out(count).name = name;
        out(count).starttime = starttime;
        out(count).endtime = endtime;
        out(count).range = ranges(i,:);
    else
        count = count+1;
        out(count).name = [];
        out(count).starttime = [];
        out(count).endtime = [];
        out(count).range = [];
    end
end

