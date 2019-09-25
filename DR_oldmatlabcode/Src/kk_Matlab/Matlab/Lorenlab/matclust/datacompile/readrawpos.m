function posstruct = readrawpos(posfiles,starttime, endtime)


startnum = timetrans({starttime},10000,2);
endnum = timetrans({endtime},10000,2);
found = 0;
times = [];
pos = [];
for i = 1:length(posfiles)
   posfilename = posfiles(i).name;
   [tmptimes, tmppos] = getposinfo(posfilename);
   times = [times; tmptimes];
   pos = [pos; tmppos];
end
times = double(times);
pos = double(pos);
[times, sortindex] = unique(times);
pos = pos(sortindex,:);


index = find((times >= startnum) & (times <= endnum));
times = times/10000;
posstruct.data = [times(index) pos(index,:)];
posstruct.fields = 'timestamp x1 y1 x2 y2';
posstruct.descript = {['position data from ', starttime, ' to ', endtime]};



