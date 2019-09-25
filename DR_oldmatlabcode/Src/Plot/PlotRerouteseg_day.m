function [ output_args ] = PlotRerouteseg_day( animaldir,animalname,day,barrier)
% Plots the segment occupied at each time point for each run
% animaldir  - directory where the output from batchDIOextract is no need to
%              enter full path eg '/F12_/'
% animalname - name of the animal, eg F1
% days       - vector of day or days to be analysed, eg [1], [1:3]
% barrier    - vector containing epoch and path on which barriers are placed
%              first number indicates epoch, second and third specify the
%              intersections between which a barrier is placed
%              eg. [2 1 2]
%

% get the current directory to go back to after the function is done
currentdir = pwd;
% ---- File loading ----
% See if day number needs 0
dsz = '';
if (day < 10)
    dsz = '0';
end
% Specify data directory and load the file
datadir = '/data14/jai/';
%datadir = '/home/daliu/tmp/';
% Converts the day and epoch into text
dayt = num2str(day);
sfilename = strcat(datadir,animaldir,'/',animaldir(1:end-1),'data',dsz,dayt,'.mat');
posfilename = strcat(datadir,animaldir,'/',animaldir(1:end-1),'linpos',dsz,dayt,'.mat');
%posfilename = strcat(datadir,animaldir,'/',animaldir(1:end-1),'linpos',dsz,dayt,'.mat');

% Load the mat file
load(sfilename);
load(posfilename);
% plot trajectory between times
for k=1:size(data{1,day},2);
    if isempty(data{1,day}{1,k});
        k=k+1;
    elseif strcmp(data{1,day}{1,k}.Stats.epochtype,'Sleep')==1;
        k=k+1;
        
    else
        seg=[];
        if strcmp(data{1,day}{1,k}.Stats.animal,animalname)==1;
            figure;
            epoch = k;
            epocht = num2str(epoch);
            %rewardedevents=data{1,day}{1,epoch}.Events.Wellevents((data{1,day}{1,epoch}.Events.Wellevents(:,3)~=-50),:);
            rewardedevents=data{1,day}{1,epoch}.Events.Welltriggers; % reward events as reward release rather than reward trigger
            nrows=ceil(size(rewardedevents,1)/5);
            % plotting paramters
            columns=5;
            gap_h=0.01;
            gap_w=0.01;
            marg_h=[0.01 0.1];
            marg_w=[0.01 0.01];
            ha = tight_subplot(nrows, columns, [gap_h gap_w], marg_h,marg_w);
            
            %             for i=1:size(data{1,day}{1,epoch}.Run,1);
            %                 tstart=data{1,day}{1,epoch}.Run(i,3);
            %                 tend=data{1,day}{1,epoch}.Run(i,4);
            
            for i=1:size(rewardedevents,1)-1;
                timeind=rewardedevents(i,1)+rewardedevents(i,4);
                tstartref=min(find(data{1,day}{1,epoch}.Pos.correcteddata(:,1)*10000>(timeind))); % for some reason +10000 seems to be accurate rather than no offset
                timeind2=rewardedevents(i+1,1)+10000; % when rat first triggers reward well, location is not within defined well region, so take location at 10000 later
                tendref=max(find(data{1,day}{1,epoch}.Pos.correcteddata(:,1)*10000<(timeind2)));
                
                
                % path identity matrix
                pathid=[];
                %                 % get segment info from linpos file
                %                 startind=find(linpos{1,day}{1,epoch}.statematrix.time(:,1)*10000>tstart);
                %                 tstartref=min(startind);
                %                 endind=find(linpos{1,day}{1,epoch}.statematrix.time(:,1)*10000<tend);
                %                 tendref=max(endind);
                % get segment info between times
                seg=linpos{1,day}{1,epoch}.statematrix.segmentIndex(tstartref:tendref,1);
                if ~isempty(seg);
                    
                    % get exit times during period
                    estartref=find(linpos{1,day}{1,epoch}.trajmatrix(:,1)*10000>timeind+10000);
                    [val eendref]=min(abs(linpos{1,day}{1,epoch}.trajmatrix(:,1)*10000-timeind2));
                    if ~isempty(estartref);
                        exitref=linpos{1,day}{1,epoch}.trajmatrix(estartref(1,1):eendref,:);
                    else i=size(data{1,day}{1,epoch}.Run,1);
                    end
                    % concatenate columns to generate integer representation of
                    % path
                    for q=1:size(exitref,1);
                        exitref(q,4)=str2num(sprintf('%d%d',exitref(q,2),exitref(q,3)));
                    end
                    % remove segments containing barrier
%                     for b=1:size(barrier,1);
%                         if barrier(b,1)==epoch;
%                             fdirect=str2num(sprintf('%d%d',barrier(b,3),barrier(b,3)));
%                             forward=find(exitref(:,4)~=fdirect);
%                             rdirect=str2num(sprintf('%d%d',barrier(b,2),barrier(b,2)));
%                             reverse=find(exitref(forward,4)~=rdirect);
%                             exitref=exitref(reverse,:);
%                             
%                             
%                         end
%                         b=b+1;
%                     end
                    
                    for b=size(exitref,1);
                        currpathid=exitref(b,2:3);
                        pathid=[pathid currpathid];
                        b=b+1;
                    end
                    
                    % remove incomplete path crossings
                    
                    exitref=exitref(find(exitref(:,2)-exitref(:,3)~=0),:);
                    
                    
                    axes(ha(i));
                    plot(seg,'.','MarkerEdgeColor','r','MarkerSize',2);
                    ylim([0 size(linpos{1,day}{1,epoch}.segmentInfo.segmentLength,2)]);
                    xlim([0 size(seg,1)]);
                    xlimval=get(gca,'xlim');
                    ylimval=get(gca,'ylim');
                    seglabel=[0.05*(xlimval(2)-xlimval(1))+xlimval(1) 0.85*(ylimval(2)-ylimval(1))+ylimval(1)];
                    % counts number of segments using size of existref, or
                    % times rat leaves the intersection zones
                    text(seglabel(1),seglabel(2),sprintf('%s',num2str(size(exitref,1))),'FontSize',10);
                    set(gca,'xtick',[]);
                    i=i+1;
                end
            end
            % Print title for whole figure
            ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 ...
                1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
            text(0.5, 0.98,sprintf('%s day %s epoch %s',animalname,dayt,epocht),'HorizontalAlignment'...
                ,'center','VerticalAlignment', 'top');
            k=k+1;
        else display(sprintf('No %s data in epoch %s ',animalname,num2str(k)));
            k=k+1;
        end
        
    end
end

end
