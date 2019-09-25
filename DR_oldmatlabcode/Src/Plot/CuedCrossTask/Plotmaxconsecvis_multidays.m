function [] = Plotconsecvis_multidays(datadirectory,animalname,days,runtype)
% Plots the number of consecutive correct visits for all epochs of days
%
% datadirectory - directory where the output from batchDIOextract is
%                 located, eg '/data14/jai/F12_/'
% animalname    -  name of the animal, eg 'F1'
% days          - days to be analysed, eg 1
% runtype       - 'home' homebound runs 
%                 'out'  outbundruns 
%                 'all'  all runs

% 

% load file for days

% all visits across days
allmax=[];

for i=1:size(days,2);

% open files
dayt = num2str(days(1,i));
datafilepref = datadirectory(end-5:end-1);
filename=strcat(datafilepref,num2str(days(1,i)),'.mat');
datafile = strcat(datadirectory,filename);
load(datafile);

% find no. epochs

j =  size(Data{1,days(1,i)},2);

% find data corresponding to animal

%epochs for each day
epochmax=[];

%epochs for each day
for k = 1:j; 
    cnt=[];
    if isempty(Data{1,days(1,i)}{1,k});
        k=k+1;
    else
        runmat =[];
        if strcmp(Data{1,days(1,i)}{1,k}.Stats.animal(1,1),animalname)==1;
            % get wells
            wells=[];
            for l=1:size(Data{1,days(1,i)}{1,k}.Config,2);
                wells(1,l)=Data{1,days(1,i)}{1,k}.Config{1,l}(2,1);
            end
            % which type of run to analyse
            switch runtype;
                case 'out';
            % get all outbound runs
                runmat = Data{1,days(1,i)}{1,k}.Run(find(Data{1,days(1,i)}{1,k}.Run(:,1)==wells(1,1)),7);
                percentcorrect = Data{1,days(1,i)}{1,k}.Stats.outboundcorrect(1,1);
                case 'home';
            % get all homebound runs
                runmat = Data{1,days(1,i)}{1,k}.Run(find(Data{1,days(1,i)}{1,k}.Run(:,1)~=wells(1,1)),7);
                percentcorrect = Data{1,days(1,i)}{1,k}.Stats.inboundcorrect(1,1);
            % get all runs
                case 'all';
                runmat = Data{1,days(1,i)}{1,k}.Run(:,7);
                percentcorrect = mean(Data{1,days(1,i)}{1,k}.Run(:,7));
            end       
            % count consecutive corrects
            r=1;
            while r<=size(runmat,1);
                 currcnt=[];
                 currmax=[];
                if runmat(r,1)==1;
                    % find next different value
                    nextlist=runmat(r+1:size(runmat,1),1);
                    m =find(nextlist==0);
                    if size(m,1)==0;
                        r=size(runmat,1)+1; 
                    else
                        if m(1,1)==1;
                            currcnt=1;
                            r=r+2;
                        else currcnt=m(1,1);
                        r=m(1,1)+r;
                        end
                    end
                
                else r=r+1;
                   currcnt=[]; 
                end
                cnt=[cnt;currcnt];
                currmax=[max(cnt)];
                
            end
            % add epochs to another
            epochmax=[epochmax;currmax];
        else k=k+1;
        end   
    end 
end
daymax=[max(epochmax)];
allmax=[allmax;daymax];
i=i+1;
end

% plot histogram

figure;
bar(allmax,'k');

%figure;

%plot(allmax,'s',...
%    'MarkerEdgeColor','r',...
%    'MarkerFaceColor','r',...
%    'MarkerSize',10);

title({sprintf('Maximum no. consecutive correct %s runs for %s', runtype, animalname);...
    sprintf('days %s to %s',num2str(min(days)),num2str(max(days)))});

set(gca,'xtick',1:size(days,2));
set(gca,'ytick',1:1:max(allmax));
set(gca,'XTickLabel',days);
%xlim([1 size(days,2)]);
ylim([0 max(allmax)]);


xlabel('Day');
ylabel('No. consecutive correct runs');

end







