function [out] = JY_gettrajectorysegments(index, excludetimes, data, linpos)
% gets the number of segments used in each trajectory

if ~isempty(linpos{index(1)}{index(2)})
    
warning('OFF','MATLAB:divideByZero');
segment1=[];
segment2=[];

trialsegment=cellfun(@(x) size(x,1),linpos{index(1)}{index(2)}.trialsegments(2,:))';




% print out error messages if there are potential missegmentation, ie
% missed one segment. 
% 2 segment trajectories are not possible
 % only print error if trajectories with less than 2 segments happen
    % after 1st trial or before last trial
    
nofirstorlast=trialsegment(2:end-1);
twosegmenttrajs=find(nofirstorlast<=2);

if ~isempty(twosegmenttrajs)
for i=1:size(twosegmenttrajs,1)

fprintf('Segment violation: %s day %s epoch %s traj %s - %s segs \n', ...
    data{index(1)}{index(2)}.Stats.animal,num2str(index(1)),num2str(index(2)),...
    num2str(twosegmenttrajs(i)+1),num2str(nofirstorlast(twosegmenttrajs(i))));
end
end

%work out standardised segment length

%segmentlength=linpos{index(1)}{index(2)}.segmentInfo.segmentLength/100;

segmentlength=[1 1 1 1 0.7 0.7 0.7 0.7  0.3 0.3 0.3 0.3];

standardlengths=cellfun(@(x) sum(segmentlength(x(:,4))),linpos{index(1)}{index(2)}.trialsegments(2,:));


% trials with more than one segment
greaterthanoneseg=find(trialsegment>1);

% 1st segment the rat travelled on
segment1=cellfun(@(x) x(1,4),linpos{index(1)}{index(2)}.trialsegments(2,:))';

% 2nd segment the rat travelled on
segment2=zeros(size(segment1));
segment2=cellfun(@(x) x(2,4),linpos{index(1)}{index(2)}.trialsegments(2,greaterthanoneseg))';
segment1=segment1(greaterthanoneseg,:);

segment=[segment1 segment2];






out.segments=trialsegment;
out.ntrial=size(trialsegment,1);
out.barrier=data{index(1)}{index(2)}.Run(:,6);
out.trajectory=linpos{index(1)}{index(2)}.trialsegments(2,:);
out.segment12=segment;
out.standardisedsegments=standardlengths;


else
out.segments=[];
out.ntrial=0;
out.barrier=[];  
out.segment12=[];

warning('ON','MATLAB:divideByZero');
end
end