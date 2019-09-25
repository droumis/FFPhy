



function [linpos]= updatevelocity(animaldir,animal,day);

datadir = '/data14/jai/';
dsz = '';
if (day < 10)
    dsz = '0';
end
dayt = num2str(day);
posfilename = strcat(datadir,animaldir,'/',animal,'linpos',dsz,dayt,'.mat');
load(posfilename);


for epoch=find(~cellfun('isempty',linpos{1,day}));
    Ind=[];
    
    index(1)=day;
    index(2)=epoch;
    
    segment=linpos{index(1)}{index(2)}.statematrix.segmentIndex;
    
    timestep=(1299.02270000000-1298.98940000000);
    
    smoothwidth=1;
    
    lindist=linpos{index(1)}{index(2)}.statematrix.linearDistanceToWells;

    Ind(:,1) = [1 find(diff(segment))'+1];
    Ind(:,2) = [find(diff(segment))' length(segment)];
    
    D=0;
    npoints = smoothwidth/timestep;
    filtstd = smoothwidth/(4*timestep);
    filt = gaussian(filtstd, npoints);
    
    % smoothlindist=smoothvect(lindist(:,1), filt);
    %     plot(lindist,'k');
    %     hold on;
    %     plot(smoothlindist,'r');
    
    
    for i=1:size(Ind,1);
        % get change in distance between each point
        % smooth the lindist within the segment
        
        % the default filter for smoothing motion direction is a n second long gaussian
        %with a n/4 second stdev
        
        %smoothdist = smoothvect(lindist(Ind(i,1):Ind(i,2),1), filt);
        %dD=diff(smoothdist);
        
        dD= diff(lindist(Ind(i,1):Ind(i,2),1));
        % see if the last transition is close to the start of the reference
        % intersection or away from the reference intersection
        
        if lindist(Ind(i,2))<0.5*linpos{index(1)}{index(2)}.segmentInfo.segmentLength(segment(Ind(i,2)));
            lastdistancevalue=lindist(Ind(i,2));
        else
            lastdistancevalue=linpos{index(1)}{index(2)}.segmentInfo.segmentLength(segment(Ind(i,2)))-lindist(Ind(i,2));
        end
        
        dD=[dD; lastdistancevalue];
        D=[D; dD];
    end
    
    smoothlindist=smoothvect(D, filt);
    
    %figure;
    
    
    Dv=D./timestep;
    velocity=Dv;
    smvelocity=smoothlindist./timestep;
    
    linpos{index(1)}{index(2)}.statematrix.linearVelocity=smvelocity;
end


    filename=strcat(datadir,animaldir,animal,'linpos',dsz,num2str(day));
    save(filename,'linpos');
    clear all;
end
    

% plot(velocity,'k');hold on;
% plot(smvelocity,'r');