function [] = Plot1stvisit(datadirectory,animalname,days)
% Plots the 1st outbound visit of an animal over days
% datadirectory - directory where the output from batchDIOextract is
% located, eg '/data14/jai/F12_/'
% animalname -  name of the animal, eg F1
% days - vector of day or days to be analysed, eg [1], [1:3]
% 

% load file for each day

visit =[];

for i=days;

datafilepref = datadirectory(end-5:end-1);

filename=strcat(datafilepref,num2str(i),'.mat');

datafile = strcat(datadirectory,filename);

load(datafile);

% find no. epochs

j =  size(Data{1,i},2);

% find data corresponding to animal

dayvisit=[];
for k = 1:j;
    if isempty(Data{1,i}{1,k});
        k=k+1;
    else
        if strcmp(Data{1,i}{1,k}.Stats.animal(1,1),animalname)==1;
            % get wells
            wells=[];
            for l=1:size(Data{1,i}{1,k}.Config,2);
                wells(1,l)=Data{1,i}{1,k}.Config{1,l}(2,1);
            end
           
            % get all outbound runs
            outbound = Data{1,i}{1,k}.Run(find(Data{1,i}{1,k}.Run(:,1)==wells(1,1)),:);
            
            % find no. runs to each reward well
            currvisit =[];
            for m=2:size(wells,2);
            currvisit(1,1)=i;
            currvisit(1,2)=k;
            currvisit(1,m+1)=size(find(outbound(:,3)==wells(1,m)),1);
            if isempty(find(Data{1,i}{1,k}.Wellinfo.rewardedwells(:,1)==wells(1,m)));
               currvisit(1,m+(size(wells,2)-1)+1)=0;
            else currvisit(1,m+(size(wells,2)-1)+1)=1;
            end
            m=m+1;
            end
            dayvisit=[dayvisit;currvisit]; 
        else k=k+1;
        end
        
    end
    
end
         
 visit=[visit;dayvisit];           
        
end

% convert to %

for i=1:size(visit,1);
    s=sum(visit(i,3:3+(size(wells,2)-2)));
    for j=3:3+(size(wells,2)-2);
        visit(i,j)=visit(i,j)/s;
        j=j+1;
    end
    i=i+1;
end


% rewarded well labels

rwell={};
for i=1:size(visit,1);
    for j=(3+(size(wells,2)-2)):size(visit,2);
    if visit(i,j)==1;
    rwell{i,(j-2-(size(wells,2)-2))}=['Well' num2str(j-2-(size(wells,2)-1))];
    else rwell{i,j-2-(size(wells,2)-2)}='';
    j=j+1;
    end
    
    end
    i=i+1;
end

rwelltxt=[];

for i=1:size(rwell,1);
    txt=[];
    for j=2:size(rwell,2);
    currtxt=rwell(i,j);
    txt=strcat(txt,currtxt);
    j=j+1;
    end
    rwelltxt=[rwelltxt;txt];
    i=i+1;
end

% Generate X tick labels

xtick={};
for i=1:size(visit,1);
    xtick{i,1}=[num2str(visit(i,1)) '-' num2str(visit(i,2)) '-' rwelltxt{i,1} ];
    i=i+1;
end


    
    
% generate legend

for i=2:size(wells,2);
    legendtxt{1,i-1}=['Well',num2str(i-1)];
    i=i+1;
end

figure;
p=bar('v6',visit(:,3:5),'stacked');
colormap Spring;
set(gca,'XtickLabel',xtick);
%set(get(p(1)),'EdgeColor','black','LineWidth',2);
h = legend(legendtxt,'Location','NorthEastOutside');
hold


end


