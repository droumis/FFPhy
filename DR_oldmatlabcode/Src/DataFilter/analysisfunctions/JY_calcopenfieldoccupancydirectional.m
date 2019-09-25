function [out] = JY_calcopenfieldoccupancydirectional(index, excludetimes, spikes,data,linpos, varargin)
%[out] = openfieldoccupancy(index, excludetimes, spikes,pos, options)
%
%Calculates the 2d occupancy normalized firing rate for the cell.
% Each segment is divided up into each direction, based on wellexitenter
% values in linpos
%
% excludetime -- times want to exclude from analysis
% spikes - the 'spikes' cell array for the day you are analyzing
% pos - the output of nspike_fixpos
% index - [day epoch tetrode cell]
% options:
%  binsize- the length of each spatial bin (default 1 cm)
%  std - defines the shape of the 2d gaussian used to smooth spikerate.
%              (default 1)
%  appendindex - 0 or 1, 1 appends index infront of output (default 0 )
%
%The output is a structure with fields: occupancy, bin vector x (.xticks), bin vector
%y (.yticks), bin spike count (.spikes), occ normailized firing per bin (.spikerate), and smoothed occ
% normalized firing (.smoothedspikerate).


appendindex = 0;
std = 1;
binsize = 1;
for option = 1:2:length(varargin)-1
    if isstr(varargin{option})
        switch(varargin{option})
            case 'appendindex'
                appendindex = varargin{option+1};
            case 'std'
                std = varargin{option+1};
            case 'binsize'
                binsize = varargin{option+1};
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end
    else
        error('Options must be strings, followed by the variable');
    end
end
warning('OFF','MATLAB:divideByZero');

spikesfields = spikes{index(1)}{index(2)}{index(3)}{index(4)}.fields;
spikes = spikes{index(1)}{index(2)}{index(3)}{index(4)}.data;
pos = data{index(1)}{index(2)}.Pos.correcteddata;
%direction=linpos{index(1)}{index(2)}.statematrix.wellExitEnter;
currentsegment=linpos{index(1)}{index(2)}.statematrix.segmentIndex;
headdirection=data{index(1)}{index(2)}.Pos.headdirection;

if (nargin < 6)
    std = 1;
else
    %std = user defined;
end

if (nargin < 5)
    binsize = 1;
else
    %binsize = user defined;
end

posx = 2;
posy = 3;
posindexfield = 7;

if ~isempty(spikes)
    spikes = spikes(:,[1 posx posy posindexfield]); %columns: time, x y, posindex
    
else
    spikes = [0 0 -1];
end

%load segment lookup table;

load('/home/jai/Src/NSpikeProcess/segmentlookuptable.mat');

output={};

timestep = mean(diff(pos(:,1)));;

%% filter out exclude times

indgoodpos =  ~isExcluded(pos(:,1), excludetimes);

%original posindex of included times
indgoodpos=find(indgoodpos==1);

goodpos = pos(indgoodpos,:);  %select spikes not excluded by exclude times

goodheaddirection=headdirection(indgoodpos,:);

goodcurrentsegment=currentsegment(indgoodpos,:);

goodspikes = goodpos(lookup( spikes(:,1),goodpos(:,1)),:); %select spikes not excluded by exclude times
%goodspikes=[ goodspikes gooddirection(lookup(spikes(:,1),goodpos(:,1)),:)];


%% calculate head direction relative to segment plane

if ~isempty(goodpos)
    %     minx = floor(min(tmpposition(:,1)));
    %     maxx = ceil(max(tmpposition(:,1)));
    %     binx = (minx:binsize:maxx);
    %     miny = floor(min(tmpposition(:,2)));
    %     maxy = ceil(max(tmpposition(:,2)));
    %     biny = (miny:binsize:maxy);
    
    
    
    minx = -10;
    maxx = 400*0.5;
    binx = (minx:binsize:maxx);
    miny = -10;
    maxy = 300*0.5;
    biny = (miny:binsize:maxy);
    
    % find unique directions
    
    % go through each of the segments and get head direction
    tmpoutput={};
    
    segmentlist=[];
    
    for jj=1:12;
        
        %find all positions on current segment
        currentsegmentpos=goodpos(find(goodcurrentsegment==jj),:);

        currentsegmentheaddirection=goodheaddirection(find(goodcurrentsegment==jj),:);

        % use lookuptable to find the start/end combinations for current
        % segment, this determines the direction
        
        directions=segmentlookuptable(find(segmentlookuptable(:,2)==jj),:);
        directionone=num2str(directions(directions(:,3)==1,1));
        directiontwo=num2str(directions(directions(:,3)==0,1));
        directionone=[str2num(directionone(1:1)) str2num(directionone(2:2))];
        directiontwo=[str2num(directiontwo(1:1)) str2num(directiontwo(2:2))];
        
        segmentdirections=[directionone;directiontwo];
        
        % process data for each direction
        segmentheaddirection=[];
        
        segmentends=[linpos{index(1)}{index(2)}.segmentInfo.segmentCoords(jj,1:2); linpos{index(1)}{index(2)}.segmentInfo.segmentCoords(jj,3:4)];
        
        for kk=1:2
            
            [C,minindex]=min(segmentdirections(kk,:));
            [C,maxindex]=max(segmentdirections(kk,:));
            segstart=segmentends(minindex,1:2);
            segend=segmentends(maxindex,1:2);
            segment_complex=complex(segstart(1)-segend(1),segstart(2)-segend(2));
            segment_angle=angle(segment_complex);
            
            % workout angle difference between segment and head direction
            segmentheaddirection(:,kk)=currentsegmentheaddirection-segment_angle;
        end
        
        % find the smallest difference to 0 radians and that's the
        % direction
        
        [C,actualdirection]=min(abs(segmentheaddirection),[],2);

        
        
        
        
        
        
        %directionlist=unique(gooddirection,'rows');
        
        directionlist=[1 2;2 1;2 3;3 2; 3 4;4 3; 4 1;1 4; 1 5; 5 1; 3 5;5 3; 2 5;5 2; 4 5;5 4; 1 6;6 1; 2 7; 7 2; 3 8; 8 3; 4 9;9 4];
        
        for kk=1:2;
            
            % start with first direction
            
            currentdirection=segmentdirections(kk,:);
            
            % find the matching direction from directionlist
            directionlistindex=rowfind(currentdirection,directionlist);
            
            
            
            % calculate occupancy for each direction
            
            
            
            % 2D occupancy
            % index of current directon
            
            currentdirectionind=find(actualdirection==kk);
            
            
            selectposition = currentsegmentpos(currentdirectionind,[2 3]); %xy position
            currentdirectionspiketimes=currentsegmentpos(currentdirectionind,1);
            selectspikes = goodspikes(find(rowfind(goodspikes(:,1),currentdirectionspiketimes(:,1))>0),[2 3]);
            [histout.twoDoccupancy histout.twoDxticks histout.twoDyticks] = hist2(selectposition(:,1), selectposition(:,2), binx, biny);
            
            nonzero = find(histout.twoDoccupancy ~= 0);
            
            
            g = gaussian2(std,(6*std));
            smoothedoccupancy = filter2(g, histout.twoDoccupancy);
            
            zero = find(smoothedoccupancy == 0);
            
            % linear occupancy
            % project each point to segment to get % of segment travelled
            onsegment=jj;
            % find coordinates of segment
            onsegmentends=[linpos{index(1)}{index(2)}.segmentInfo.segmentCoords(onsegment,1:2); linpos{index(1)}{index(2)}.segmentInfo.segmentCoords(onsegment,3:4)];
            % project each point of occupancy onto the segment
            onsegmentout=[];
            
            [C,minindex]=min(segmentdirections(kk,:));
            [C,maxindex]=max(segmentdirections(kk,:));
            segmentstart=segmentends(minindex,1:2);
            segmentend=segmentends(maxindex,1:2);
            
            
            segmentlength=dist(segmentstart,segmentend);
            for pp=1:size(selectposition,1)
                tmponsegmentout=JY_projectpoint(selectposition(pp,:),onsegmentends);
                %find distance to each end of the segment
                startdist=dist(tmponsegmentout(1,1:2),segmentstart);
                enddist=dist(tmponsegmentout(1,1:2),segmentend);
                onsegmentout=[onsegmentout; [startdist enddist startdist/segmentlength enddist/segmentlength ]];
                
                
            end
            
            
            
            onsegmenthistout=zeros(1,51);
            
            % histoutogram
            if ~isempty(onsegmentout)
                onsegmenthistout=hist(onsegmentout(:,3),0:0.02:1);
            end
            
            % segment of travel
            
            if ~isempty(selectspikes)
                
                %2D version
                [histout.twoDspikes BX BY] = hist2(selectspikes(:,1), selectspikes(:,2), binx, biny);
                
                % linear version
                
                onsegmentspikeout=[];
                
                for pp=1:size(selectspikes,1)
                    tmponsegmentspikeout=JY_projectpoint(selectspikes(pp,:),onsegmentends);
                    %find distance to each end of the segment
                    startdist=dist(tmponsegmentspikeout(1,1:2),segmentstart);
                    enddist=dist(tmponsegmentspikeout(1,1:2),segmentend);
                    onsegmentspikeout=[onsegmentspikeout; [startdist enddist startdist/segmentlength enddist/segmentlength]];
                    
                    
                end
                
                % histoutogram
                %find the smallest value of the two intersections, this is the
                %start
                
                [C, smallintersectionindex]=min(currentdirection);
                onsegmentspikehistout=hist(onsegmentspikeout(:,smallintersectionindex+2),0:0.02:1);
                
                % spike rate
                histout.linearspikerate=onsegmentspikehistout./(timestep*onsegmenthistout);
                
                % correct for Nan and Inf
                 histout.linearspikerate(isnan(histout.linearspikerate))=0;
                 
                histout.linearspikerate(isinf(histout.linearspikerate))=0;
                
                
                histout.twoDspikerate = zeros(size(histout.twoDspikes));
                histout.twoDsmoothedspikerate=zeros(size(histout.twoDoccupancy));
                
                histout.linearoccupancy=onsegmenthistout;
                histout.linearspikes=onsegmentspikehistout;
                
                histout.linearspacestep=segmentlength/50;
                
                %smooth linear rate
                
                
                % gaussian smoothing
                %sigma=2*output.bin_length(seg_num);
                sigma=1;
                width = round((6*sigma - 1)/2);
                support = (-width:width);
                gaussFilter = exp( -(support).^2 ./ (2*sigma^2) );
                gaussFilter = gaussFilter/ sum(gaussFilter);
                %             win=gausswin(64,1);
                %             win=win./sum(win);
                
                smoothed_placefield=conv(histout.linearspikerate,gaussFilter,'same');
                %smoothed_placefield=medfilt1(linearspikerate,3,5);
                
                histout.linearsmoothedrate=smoothed_placefield;
                 
                 
                if ~isempty(find(histout.twoDspikes>0))
                    
                    histout.twoDspikerate(nonzero) = histout.twoDspikes(nonzero) ./(timestep* histout.twoDoccupancy(nonzero));
                    %smoothed occupancy
                    
                    histout.twoDsmoothedspikerate = filter2(g,(histout.twoDspikerate)); % is this the right filter?
                    
                    %smoothedoccupancy = zeros(size(histout.twoDspikes));
                    
                    
                    histout.twoDnormsmoothedspikerate=histout.twoDsmoothedspikerate./max(histout.twoDsmoothedspikerate(:));
                    
                    histout.twoDsmoothedspikerate(zero) = -2; %no occupancy is negative and diff/darker from occupancy but no spikes
                    
                    
                    
                end
                
                
            else
                histout.twoDspikes=zeros(size(histout.twoDoccupancy));
                
                
                histout.spikerate=zeros(size(histout.twoDoccupancy));
                histout.twoDsmoothedspikerate=zeros(size(histout.twoDoccupancy));
                %histout.twoDsmoothedspikerate(zero)=-2;
                
                histout.twoDnormsmoothedspikerate=zeros(size(histout.twoDoccupancy));
                histout.linearspacestep=segmentlength/50;
                histout.linearspikes=zeros(1,51);
                histout.linearsmoothedrate=zeros(1,51);
            end
            
            tmpoutput{directionlistindex}=histout;
            segmentlist=[segmentlist;jj];
            
            
        end
    end
end



output.data=tmpoutput;
output.directionlist=directionlist;
output.segmentlist=segmentlist;


if appendindex == 0
    out = output;
elseif appendindex == 1
    output.index = index
    out = output;
end


warning('ON','MATLAB:divideByZero');