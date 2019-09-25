function [] = Plotconsecvis(datadirectory,animalname,day,runtype)
% Plots the number of consecutive correct visits for all epochs of a day
%
% datadirectory - directory where the output from batchDIOextract is
%                 located, eg '/data14/jai/F12_/'
% animalname    -  name of the animal, eg 'F1'
% days          - day to be analysed, eg 1
% runtype       - 'home' homebound runs 
%                 'out'  outbundruns 
%                 'all'  all runs

% 

% load file for each day

visit =[];
dayt = num2str(day);


datafilepref = datadirectory(end-5:end-1);

filename=strcat(datafilepref,num2str(day),'.mat');

datafile = strcat(datadirectory,filename);

load(datafile);

% find no. epochs

j =  size(Data{1,day},2);

% find data corresponding to animal

dayvisit=[];
for k = 1:j;
    cnt=[];
    if isempty(Data{1,day}{1,k});
        k=k+1;
    else
       
        if strcmp(Data{1,day}{1,k}.Stats.animal(1,1),animalname)==1;
            % get wells
            wells=[];
            for l=1:size(Data{1,day}{1,k}.Config,2);
                wells(1,l)=Data{1,day}{1,k}.Config{1,l}(2,1);
            end
            
            runmat =[];
            
            switch runtype;
            
                case 'out';
            % get all outbound runs
            runmat = Data{1,day}{1,k}.Run(find(Data{1,day}{1,k}.Run(:,1)==wells(1,1)),7);
                case 'home'; 
            % get all homebound runs
            runmat = Data{1,day}{1,k}.Run(find(Data{1,day}{1,k}.Run(:,1)~=wells(1,1)),7);
                case 'all'; 
            % get all homebound runs
            runmat = Data{1,day}{1,k}.Run(:,7);
            end
            
            % count consecutive corrects
            i=1;
            
            while i<=size(runmat,1);
                 currcnt=[];
                if runmat(i,1)==1;
                
                    % find next different value
                    nextlist=runmat(i+1:size(runmat,1),1);
                    m =find(nextlist==0);

                    if size(m,1)==0;
                        i=size(runmat,1)+1; 
                    else
                        if m(1,1)==1;
                            currcnt=1;
                            i=i+2;
                        else currcnt=m(1,1);
                        i=m(1,1)+i;
                        end
                    end
                
                else i=i+1;
                   currcnt=[]; 
                end
                cnt=[cnt;currcnt];
               
            end
            
            % plot histogram
            
            figure;
            
            hist(cnt,[1:1:max(cnt)]);
            
            %set(gca,'xtick',1:max(cnt));
            
            ylim([0 10]);
            
            xlabel('No. consecutive correct runs');
            ylabel('Count');
            
            epocht = num2str(k);
            
            title([animalname,' ','day ',dayt,' epoch ',epocht, ' ', runtype, ' bound runs']);

        else k=k+1;
        end
        
    end
    
end
         
   
        
end








