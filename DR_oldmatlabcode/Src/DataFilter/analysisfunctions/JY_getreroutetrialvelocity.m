function [out] = JY_getreroutetrialvelocity(index, excludetimes, data, linpos)
% gets the linear velocity


warning('OFF','MATLAB:divideByZero');

%linpos = linpos{index(1)}{index(2)}.statematrix.linearDistanceToWells;
%data = data{index(1)}{index(2)}.Pos.correcteddata
if ~isempty(linpos{index(1)}{index(2)})
    
    
    % 04/17/2013
    % see if processed distance has been scaled from original incorrectly
    
    maxxcorrected=max(data{index(1)}{index(2)}.Pos.correcteddata(:,2));
    maxxraw=max(data{index(1)}{index(2)}.Pos.rawpos(:,2));
    if maxxraw/maxxcorrected<1.2
        correctionfactor=1;
    else correctionfactor=2;
    end
    
    %correctionfactor=1;
    
    
    % get linear distance to wells
    lineardistancetowells=linpos{index(1)}{index(2)}.statematrix.linearDistanceToWells*correctionfactor;
    % get segment index
    segmentindex=linpos{index(1)}{index(2)}.statematrix.segmentIndex;
    
    % count the difference in linear distance between each time point but reset
    % everytime the segment index changes
    
    %index of transitions between segments
    Ind(:,1) = [1 find(diff(segmentindex))'+1];
    Ind(:,2) = [find(diff(segmentindex))' length(segmentindex)];
    
    D=0;
    
    for i=1:size(Ind,1);
        % get change in distance between each point
        dD=abs(diff(lineardistancetowells(Ind(i,1):Ind(i,2),1)));
        dD=[dD; 0];
        D=[D;dD];
        
    end
    
    vD=abs(D)/0.0333;
    
 
    velocity_original=data{index(1)}{index(2)}.Pos.correcteddata(:,5)*correctionfactor;
    % get trajectory velocity by calculating change in distance between the
    % start and end of each run
    timeind=linpos{index(1)}{index(2)}.statematrix.time;
    
    % runs
    runstart=lookup(data{index(1)}{index(2)}.Run(:,3),timeind*10000);
    runend=lookup(data{index(1)}{index(2)}.Run(:,4),timeind*10000);
    
    tD=[];
    tDc=[];
    
    for i=1:size(runstart,1);
        d=mean(abs(velocity_original(runstart(i):runend(i))));
        %     dc=sum(D(runstart(i):runend(i)))/(timeind(runend(i))-timeind(runstart(i)));
        dc=mean(vD(runstart(i):runend(i)));
        tD=[tD;d];
        tDc=[tDc;dc];
    end
    
    
    
    
    % find optimal distance
    
    % load connectivity table to find the most direct route between reward wells
    
    % load('/home/jai/Src/NSpikeProcess/optdisttable.mat');
    % activewells=data{index(1)}{index(2)}.Wellinfo.rewardedwells;
    %
    % optdist=sum(linpos{index(1)}{index(2)}.segmentInfo.segmentLength(1,optdisttable{activewells(1),activewells(2)}));
    % out.normdist=tD./optdist;
    %
    % out.optidist=optdist;
    
    
    out.velocity=velocity_original;
    out.velocity_orig=tD;
    out.velocity_calc=tDc;
    out.time=timeind;
    out.ntrial=size(data{index(1)}{index(2)}.Run(:,5),1);
    out.barrier=data{index(1)}{index(2)}.Run(:,6);
    
else
    out.velocity=[];
    out.velocity_orig=[];
    out.velocity_calc=[];
    out.time=[];
    out.ntrial=0;
    out.barrier=[];
    
end


warning('ON','MATLAB:divideByZero');