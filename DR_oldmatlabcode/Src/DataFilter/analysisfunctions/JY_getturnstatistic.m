function [out] = JY_getturnstatistic(index, excludetimes, data, linpos)
% gets the number of segments used in each trajectory

if ~isempty(linpos{index(1)}{index(2)})
    
warning('OFF','MATLAB:divideByZero');
segment1=[];
segment2=[];


%load turn table

load('/home/jai/Src/NSpikeProcess/JY_turntable.mat');

segmenttransition=cell2mat(cellfun(@(x) x(:,2:end),linpos{index(1)}{index(2)}.trialsegments(2,:),'UniformOutput', false)');

change=[segmenttransition(1:end-1,3) segmenttransition(2:end,3)];  

idx=sub2ind([12 12],change(:,1),change(:,2));

turntablerep=reshape(turntable,[],1);

turn=turntablerep(idx);

% turns at center
% must be coming from segments 5 6 7 8

centerindex=find(segmenttransition(1:end-1,2)==5);

out.centerturnsall=size(turn(centerindex),1);
out.centerturnleft=size(find(turn(centerindex,1)==1),1);
out.centerturnright=size(find(turn(centerindex,1)==2),1);
out.centerturnstraight=size(find(turn(centerindex,1)==3),1);



% turns on corners

cornerindex=find(ismember(segmenttransition(1:end-1,2),[1 2 3 4])==1);

out.cornerturnsall=size(turn(cornerindex),1);
out.cornerturnleft=size(find(turn(cornerindex,1)==1),1);
out.cornerturnright=size(find(turn(cornerindex,1)==2),1);
out.cornerturnstraight=size(find(turn(cornerindex,1)==3),1);



% count turns
out.duration=data{index(1)}{index(2)}.Stats.duration/10000;
out.leftturn=size(find(turn==1),1);
out.rightturn=size(find(turn==2),1);
out.straight=size(find(turn==3),1);
out.all=[change turn];





warning('ON','MATLAB:divideByZero');
end
end