function posstruct = readrawpos_steve(posfiles, timestamp_range)


startnum = timestamp_range(1);
endnum = timestamp_range(2);
found = 0;
times = [];
pos = [];
for i = 1:length(posfiles)
   posfilename = posfiles(i).name;
   rawdata = load(posfilename);
   tmptimes = rawdata.rawpos.timestamp;
   tmppos = [rawdata.rawpos.xfront rawdata.rawpos.yfront rawdata.rawpos.xback rawdata.rawpos.yback];
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
posstruct.descript = {['position data from ', ...
  ts2str(uint32(startnum)), ' to ', ...
  ts2str(uint32(endnum))]};



