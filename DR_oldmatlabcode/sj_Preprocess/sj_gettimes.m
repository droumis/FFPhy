function out = sj_gettimes(timesfile)
%% I am changing it to get times from the aranges field instead of names:
%% Output will be in 1/10000 s units 


load(timesfile)
timelist = names;
count = 0;
for i = 2:size(ranges,1)
    count=count+1;
    out(count).starttime = ranges(i,1);
    out(count).endtime = ranges(i,2);
    % To get name
     spaces = strfind(timelist{i},' ');
     lastspace = spaces(end);
     firstspace = spaces(1);
     name = timelist{i}(firstspace:lastspace);
     name(strfind(name,' ')) = '';
     out.name = lower(name);
end
    
% for i = 2:size(names,1) 
%     dash = strfind(timelist{i},'-');
%     if ~isempty(dash)
%         count = count+1;
%         spaces = strfind(timelist{i},' ');
%         lastspace = spaces(end);
%         firstspace = spaces(1);
%         name = timelist{i}(firstspace:lastspace);
%         name(strfind(name,' ')) = '';
%         name = lower(name);
%         starttime = timelist{i}(lastspace+1:dash-1);
%         endtime = timelist{i}(dash+1:end);
%         out(count).name = name;
%         out(count).starttime = starttime;
%         out(count).endtime = endtime;        
%     else
%         count = count+1;
%         out(count).name = [];
%         out(count).starttime = [];
%         out(count).endtime = [];
%     end
% end

