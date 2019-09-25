function [out] = JY_gettrajectorydistance(index, excludetimes, data, linpos)
% gets the change in linear distance between every time point in an epoch


warning('OFF','MATLAB:divideByZero');

%linpos = linpos{index(1)}{index(2)}.statematrix.linearDistanceToWells;
%data = data{index(1)}{index(2)}.Pos.correcteddata;

if ~isempty(linpos{index(1)}{index(2)})
    
    % 04/17/2013
    % see if processed distance has been scaled from original incorrectly
       maxxcorrected=max(data{index(1)}{index(2)}.Pos.correcteddata(:,2));
    maxxraw=max(data{index(1)}{index(2)}.Pos.rawpos(:,2));
    if maxxraw/maxxcorrected<1.2
        correctionfactor=1;
    else correctionfactor=2;
    end
  
    
    % calculate linear distance
    % get linear distance to wells
    lineardistancetowells=linpos{index(1)}{index(2)}.statematrix.linearDistanceToWells*correctionfactor;
    % get segment index
    segmentindex=linpos{index(1)}{index(2)}.statematrix.segmentIndex;
    segmentInfo=linpos{index(1)}{index(2)}.segmentInfo;
    
    % count the difference in linear distance between each time point but reset
    % everytime the segment index changes
    
    %index of transitions between segments
    Ind(:,1) = [1 find(diff(segmentindex))'+1];
    Ind(:,2) = [find(diff(segmentindex))' length(segmentindex)];
    
    D=0;
    
    for i=1:size(Ind,1);
        % get change in distance between each point
        dD=abs(diff(lineardistancetowells(Ind(i,1):Ind(i,2),1)));
        
        if lineardistancetowells(Ind(i,2))<0.5*segmentInfo.segmentLength(segmentindex(Ind(i,2)))*correctionfactor;
            lastdistancevalue=lineardistancetowells(Ind(i,2));
        else
            lastdistancevalue=segmentInfo.segmentLength(segmentindex(Ind(i,2)))*correctionfactor-lineardistancetowells(Ind(i,2));
        end
        
        
        
        dD=[dD; lastdistancevalue];
        D=[D;dD];
        
    end
    displacement=0;
    % calculate 2D displacement
    for i=1:length(data{index(1)}{index(2)}.Pos.correcteddata)-1
        displacement=[displacement;dist(data{index(1)}{index(2)}.Pos.correcteddata(i+1,2:3),data{index(1)}{index(2)}.Pos.correcteddata(i,2:3))*correctionfactor];
    end
    
    
    
    
    
    % get linear trajectory distance by calculating change in distance between the
    % start and end of each run
    timeind=linpos{index(1)}{index(2)}.statematrix.time;
    
    % trials
    runstart=lookup(data{index(1)}{index(2)}.Run(:,3),timeind*10000);
    runend=lookup(data{index(1)}{index(2)}.Run(:,4),timeind*10000);
    tD=[];
    aD=[];
    for i=1:size(runstart,1);
        d=sum(D(runstart(i):runend(i)));
        ad=sum(displacement(runstart(i):runend(i)));
        tD=[tD;d];
        aD=[aD;ad];
    end
    
    % find number of segments for each trial
    
    trialsegments=cellfun(@(x) size(x,1),linpos{index(1)}{index(2)}.trialsegments(2,:))';
    trialsegments_details=linpos{index(1)}{index(2)}.trialsegments(2,:);
    
    
    
    
    % find optimal distance
    
    % load connectivity table to find the most direct route between reward wells
    
    load('/home/jai/Src/NSpikeProcess/optdisttable.mat');
    activewells=data{index(1)}{index(2)}.Wellinfo.rewardedwells;
    
    
    optdist=sum(linpos{index(1)}{index(2)}.segmentInfo.segmentLength(1,optdisttable{activewells(1),activewells(2)}...
        (optdisttable{activewells(1),activewells(2)}(:,1)>0,1)))*correctionfactor;
    optdist_barrier=sum(linpos{index(1)}{index(2)}.segmentInfo.segmentLength(1,optdisttable{activewells(1),activewells(2)}(:,2)))*correctionfactor;
    tD2(data{index(1)}{index(2)}.Run(:,6)==0)=tD(data{index(1)}{index(2)}.Run(:,6)==0)./optdist;
    tD2(data{index(1)}{index(2)}.Run(:,6)==1)=tD(data{index(1)}{index(2)}.Run(:,6)==1)./optdist_barrier;
    
    % find optimal distance from actual rat distance
    % trialsegmentused=cellfun(@(x) size(x,1),linpos{index(1)}{index(2)}.trialsegments(2,:))';
    % seg=size(optdisttable{activewells(1),activewells(2)}(:,1)>0,1);
    % seg_barrier=size(optdisttable{activewells(1),activewells(2)}(:,2)>0,1);
    % optdist=mean(tD(trialsegmentused==seg));
    % optdist_barrier=mean(tD(trialsegmentused==seg_barrier));
    
    tD2(data{index(1)}{index(2)}.Run(:,6)==0)=tD(data{index(1)}{index(2)}.Run(:,6)==0)./optdist;
    tD2(data{index(1)}{index(2)}.Run(:,6)==1)=tD(data{index(1)}{index(2)}.Run(:,6)==1)./optdist_barrier;
    
    out.normdist=tD2';
    out.optdist=optdist;
    out.optdist_barrier=optdist_barrier;
    out.distance=D;
    out.time=timeind;
    out.trajectorydistance=tD;
    out.ntrial=size(data{index(1)}{index(2)}.Run(:,5),1);
    out.barrier=data{index(1)}{index(2)}.Run(:,6);
    out.trialsegments=trialsegments;
    out.trialsegments_details=trialsegments_details;
    out.trajectorydisplacement=aD;
    out.meanepochtrajectorydisplacement=mean(aD);
    out.meantrialsegments=mean(trialsegments);
    
    
else
    out.normdist=[];
    out.optdist=[];
    out.optdist_barrier=[];
    out.distance=[];
    out.time=[];
    out.trajectorydistance=[];
    out.ntrial=0;
    out.barrier=[];
    out.trialsegments=[];
    out.trialsegments_details=[];
    out.trajectorydisplacement=[];
end

warning('ON','MATLAB:divideByZero');
