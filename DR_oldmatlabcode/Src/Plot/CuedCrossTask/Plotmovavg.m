function [] = Plotmovavg(datadirectory,animalname,day,runtype,step)
% Plots the moving average for all epochs of a day
% 
% datadirectory - directory where the output from batchDIOextract is
%                 located, eg '/data14/jai/F12_/'
% animalname    -  name of the animal, eg 'F1'
% days          - vector of day or days to be analysed, eg 1
% runtype       -  'home' homebound runs 
%                  'out' outbund runs
%                  'all' all runs
% step          - number of runs for average eg 5
% 

% load file for each day

visit =[];

for i=day;
    
dayt = num2str(day);


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
            
            runmat =[];
            
            switch runtype;
            
                case 'out';
            % get all outbound runs
                runmat = Data{1,i}{1,k}.Run(find(Data{1,i}{1,k}.Run(:,1)==wells(1,1)),7);
                percentcorrect = Data{1,i}{1,k}.Stats.outboundcorrect(1,1);
                
                case 'home';
                  
            % get all homebound runs
                runmat = Data{1,i}{1,k}.Run(find(Data{1,i}{1,k}.Run(:,1)~=wells(1,1)),7);
                percentcorrect = Data{1,i}{1,k}.Stats.inboundcorrect(1,1);
                
            % get all runs
                case 'all';
                runmat = Data{1,i}{1,k}.Run(:,7);
                percentcorrect = mean(Data{1,i}{1,k}.Run(:,7));
            
            
            end
            epocht = num2str(k);
            avgrunmat = moving(runmat,step);
            
            figure;
            
            plot(avgrunmat,'Color','red','LineWidth',4);
            
            title([animalname,' ','day ',dayt,' epoch ',epocht,' ', runtype, ' bound runs', char(10), 'moving average of ',num2str(step),' runs']);
            
            text(1,0.95,['fraction total runs correct: ',num2str(percentcorrect,2)]);
            
            xlabel('runs');
            ylabel('correct runs');
            
            ylim([0 1]);
            

        else k=k+1;
        end
        
    end
    
end
         
   
        
end


end


