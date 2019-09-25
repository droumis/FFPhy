function [out] = JY_getturnstatistic(index, excludetimes, data, linpos)
% gets the number of segments used in each trajectory

if ~isempty(linpos{index(1)}{index(2)})
    
    warning('OFF','MATLAB:divideByZero');
    segment1=[];
    segment2=[];
    
    
    %load turn table
    
    load('/home/jai/Src/NSpikeProcess/JY_turntable.mat');
    
    segmenttransition=cell2mat(cellfun(@(x) x(:,1:end),linpos{index(1)}{index(2)}.trialsegments(2,:),'UniformOutput', false)');
    
    % times exiting intersections
    %[start end from to intersection];
    
    change=[ segmenttransition(1:end-1,1) segmenttransition(2:end,1)...
        segmenttransition(1:end-1,4) segmenttransition(2:end,4) segmenttransition(1:end-1,3)];
    
    
    % get coordinates of intersections
    
    intcoord=linpos{index(1)}{index(2)}.segmentInfo.segmentCoords;
    coordindx=[1 1; 2 1; 3 1; 4 1; 5 3; 9 3; 10 3; 11 3; 12 3];
    coordindy=[1 2; 2 2; 3 2; 4 2; 5 4; 9 4; 10 4; 11 4; 12 4];
    
    lincoordindx=sub2ind([12 4],coordindx(:,1),coordindx(:,2));
    lincoordindy=sub2ind([12 4],coordindy(:,1),coordindy(:,2));
    intcoordr=reshape(intcoord,[],1);
    intersections=[intcoordr(lincoordindx) intcoordr(lincoordindy)];
    
    % loop through each choice point
    % get position indices for all points x distance away from intersection
    % get velocity indices for all points x distance away from intersection
    
    timeind=linpos{index(1)}{index(2)}.statematrix.time;
    pos=data{index(1)}{index(2)}.Pos.correcteddata;
    velocity_lin=linpos{index(1)}{index(2)}.statematrix.linearVelocity;
    velocity=data{index(1)}{index(2)}.Pos.correcteddata(:,5);
    
    out.pos=pos;
    
    % convert times into indices
    
    startind=lookup(change(:,1),timeind);
    endind=lookup(change(:,2),timeind);
    
    % generate array for each choice event at cutoffdistance from intersection
    cutoffdist=10;
    
    out.distancecutoff=cutoffdist;
    
    choicearray=cell(1,size(change,1));
    
    for ii=1:size(change,1)
        tmppos=pos(max(1,startind(ii)-100):min(size(pos,1),endind(ii)+100),2:3);
        
        if isempty(tmppos)
            index(1)
            index(2)
            ii 

        end

  
        
        
        
        tmpindx=[max(1,startind(ii)-100):min(size(pos,1),endind(ii)+100)];
        currintersectioncoord=intersections(change(ii,5),:);
        
        
        
        tmppos2=mat2cell(tmppos, size(tmppos,1),[1]);
        
        % calculate distance to intersection
        distancetointersect=cellfun(@(x) pdist([currintersectioncoord; x(1) x(2)]),tmppos2);
        distanceind=tmpindx(find(distancetointersect<cutoffdist));
        choicearray{ii}.disttointersection=distancetointersect;
        choicearray{ii}.postartend=[distanceind(1) distanceind(end)];
        
       
        
        %closest point to intersection
        [C, minIndx]=min(distancetointersect);
        choicearray{ii}.intersectionpointindex=tmpindx(minIndx);
        % mean velocity during period
        
        choicearray{ii}.meanvelo=mean(abs(velocity(max(1,distanceind(1)):min(size(pos,1),distanceind(end)))));
        
%         
%         plot(tmppos(:,1),tmppos(:,2),'.');
%         hold on;
%         plot(currintersectioncoord(1),currintersectioncoord(2),'+r');   
%         plot(pos(distanceind(1):distanceind(end),2),pos(distanceind(1):distanceind(end),3),'r.');
%         hold on;
%         
%         figure;
%         plot(distancetointersect,velocity(max(1,startind(ii)-100):min(size(pos,1),endind(ii)+100)));
%         hold on;
%         plot(distancetointersect,velocity_lin(max(1,startind(ii)-100):min(size(pos,1),endind(ii)+100)),'g');
%         xlabel('Distance to intersection');
%         ylabel('Speed');
%         figure;
%         plot(timeind(max(1,startind(ii)-100):min(size(pos,1),endind(ii)+100)),...
%             velocity(max(1,startind(ii)-100):min(size(pos,1),endind(ii)+100)));
%         
%         hold on;
%          plot(timeind(max(1,startind(ii)-100):min(size(pos,1),endind(ii)+100)),...
%             velocity_lin(max(1,startind(ii)-100):min(size(pos,1),endind(ii)+100)),'g');
%           xlabel('Time');
%         ylabel('Speed');
%         
%         close all;
       % end
    end
    
    %get data
    
    out.choicevelocity=cellfun(@(x) x.meanvelo,choicearray)';
    out.choicestartendindx=cell2mat(cellfun(@(x) x.postartend,choicearray,  'UniformOutput', false)');
    out.choicecenterindx=cellfun(@(x) x.intersectionpointindex,choicearray)';
    
    
    % get linear index of positions in turntable reference to classify turn
    % into left right or center
    idx=sub2ind([12 12],change(:,3),change(:,4));
    
    turntablerep=reshape(turntable,[],1);
    
    turn=turntablerep(idx);
    
    % turns at center
    % must be coming from segments 5 6 7 8
     centerindex=find(change(:,5)==5);
     
    out.centerturnsall=size(turn(centerindex),1);
    out.centerturnleft=size(find(turn(centerindex,1)==1),1);
    out.centerturnright=size(find(turn(centerindex,1)==2),1);
    out.centerturnstraight=size(find(turn(centerindex,1)==3),1);
    
    out.centerindex=centerindex;
    out.centerturnleftindx=find(turn(centerindex,1)==1);
    out.centerturnrightindx=find(turn(centerindex,1)==2);
    out.centerturnstraightindx=find(turn(centerindex,1)==3);
    
    out.centerleftvelocity=mean(out.choicevelocity(out.centerturnleftindx));
    out.centerrightvelocity=mean(out.choicevelocity(out.centerturnrightindx));
    out.centerstraightvelocity=mean(out.choicevelocity(out.centerturnstraightindx));
    
    
    % turns on corners
    
    cornerindex=find(ismember(change(:,5),[1 2 3 4])==1);
    
    out.cornerturnsall=size(turn(cornerindex),1);
    out.cornerturnleft=size(find(turn(cornerindex,1)==1),1);
    out.cornerturnright=size(find(turn(cornerindex,1)==2),1);
    out.cornerturnstraight=size(find(turn(cornerindex,1)==3),1);
    
    out.cornerindex=cornerindex;
    out.cornerturnleftindx=find(turn(cornerindex,1)==1);
    out.cornerturnrightindx=find(turn(cornerindex,1)==2);
    out.cornerturnstraightindx=find(turn(cornerindex,1)==3);
    
    out.cornerleftvelocity=mean(out.choicevelocity(out.cornerturnleftindx));
    out.cornerrightvelocity=mean(out.choicevelocity(out.cornerturnrightindx));
    out.cornerstraightvelocity=mean(out.choicevelocity(out.cornerturnstraightindx));
    
    
    
    
    % count turns
    out.duration=data{index(1)}{index(2)}.Stats.duration/10000;
    out.leftturn=size(find(turn==1),1);
    out.rightturn=size(find(turn==2),1);
    out.straight=size(find(turn==3),1);
    out.all=[change turn];
    
    
    
    
    
    warning('ON','MATLAB:divideByZero');
end
end