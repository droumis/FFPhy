function [] = fellowsgen(iter,repeat,well1,well2,filename)
% Generates a sequence of iter repeats with 0 1 0 2 etc
% Specifies the proportion of 1s and 2s in each sequence
% Generates 2 files - .txt contains the sequence, _stat.txt has statistics
% specified in the reference for each sequence
% iter -  number of iterations, 0 1 is considered one iteration
% repeat - maximum number of repeats of one well, currently only 2 or 3
%          repeats are possible
% well1 - well number
% well2 - well number
% filename - file name of sequence
% Based on Fellows series, Fellows 1967, Psych Bull
currdir = pwd;

ones = iter;
twos = iter-ones;
PAreinf1=iter*33;
PAreinf2=iter*33;
PPreinf1=iter*33;
PPreinf2=iter*33;
WSHreinf1=iter*33;
WSHreinf2=iter*33;
WSTreinf1=iter*33;
WSTreinf2=iter*33;

% Set proportion of 

factor = 0.2;  

% keyboard
 
% while abs(ones-twos) >=0.05*iter && max(abs(PAreinf1),abs(PAreinf2))>=factor*iter && max(abs(PPreinf1),abs(PPreinf2))>=factor*iter && max(abs(WSHreinf1),abs(WSHreinf2))>=factor*iter && max(abs(WSTreinf1),abs(WSTreinf2))>=factor*iter;
% while (abs(ones-twos) >=0.05*iter) && (max(abs(PAreinf1),abs(PAreinf2))>=factor*iter) && (max(abs(PPreinf1),abs(PPreinf2))>=factor*iter) && (max(abs(WSHreinf1),abs(WSHreinf2))>=factor*iter) && (max(abs(WSTreinf1),abs(WSTreinf2))>=factor*iter);  
% while mean([(abs(ones-twos)) (max(abs(PAreinf1),abs(PAreinf2))) (max(abs(PPreinf1),abs(PPreinf2))) (max(abs(WSHreinf1),abs(WSHreinf2))) (max(abs(WSTreinf1),abs(WSTreinf2)))] > [0.05*iter factor*iter factor*0*iter factor*30*iter factor*30*iter])<1;
    while mean([(abs(ones-twos)) (max(abs(PAreinf1),abs(PAreinf2))) (max(abs(PPreinf1),abs(PPreinf2)))] > [0.05*iter factor*iter factor*iter])>=0.3;
   
   seq = [];

    for num = 1:iter;
    
    currseq =[];

    currseq(1,1)=0;

    r=randi(2,[1,1]);
    
    nr = nonzeros(seq);
    
    if size(nr,1)<=3;
       
        currseq(2,1)=r;
        
    else
    
    last1=nr(size(nr,1),1);
    last2=nr(size(nr,1)-1,1);
    last3=nr(size(nr,1)-2,1);
    
        
        if repeat==2;
             repeatn=[last1,last2];
        elseif repeat==3;
            repeatn=[last1,last2,last3];
        end
    
            if repeatn(1,1)==mean(repeatn);
       
                if last1 == 2;
            
                currseq(2,1)=1;
                
                else
                
                currseq(2,1)=2;
                end
       
            else currseq(2,1)=r;
            end
  
        
    end
    
seq = [seq; currseq];
    
    end

    %seq

num = num +1;

ones =size(find(seq==1),1);
twos =size(find(seq==2),1);


%end
% Remove 0s

seqnonezeros = seq((seq~=0),1);

% Position perseveration

% Position perserveration scores

PPperf1 = size(find(seqnonezeros==1),1)
PPperf2 = size(find(seqnonezeros==2),1)

PP = zeros(size(seqnonezeros,1),4);
PP(:,1)=1;
PP(:,3)=2;
PP(:,2)=1;
PP(:,4)=2;

% Compare sequence with PP expected sequence

for PPn=1:size(PP,1);
    if isequal(seqnonezeros(PPn,1),PP(PPn,1));
        PP(PPn,2)=1;
    else PP(PPn,2)=-1;
    end
    PPn= PPn+1;
end

for PPn=1:size(PP,1);
    if isequal(seqnonezeros(PPn,1),PP(PPn,3));
        PP(PPn,4)=1;
    else PP(PPn,4)=-1;
    end
    PPn= PPn+1;
end

% Score sequence for PP

for n= [2 4];
    for PPn = 2:size(PP,1);
    
    if PP(PPn-1,n)==1 && PP(PPn,n)==1;
        PP(PPn,n)=2;
    elseif PP(PPn-1,n)==2 && PP(PPn,n)==1 ;
        PP(PPn,n)=3;
    elseif PP(PPn-1,n)==-1 && PP(PPn,n)==-1;
        PP(PPn,n)=-2;
    elseif PP(PPn-1,n)==-2 && PP(PPn,n)==-1;
        PP(PPn,n)=-3;
    end
    PPn = PPn+1;
    end
    n = 4;
end


% Perseveration reinforcement score

PPreinf1= sum(PP(:,2))
PPreinf2= sum(PP(:,4))


% Position alternation
PA = zeros(size(seqnonezeros,1),4);

PA(1,1)=1;
PA(1,3)=2;

% Generate PA expected sequence

for n=[1 3];
    for PAn = 2:size(PA,1);
    if PA(PAn-1,n)==1;
        PA(PAn,n)=2;
    else PA(PAn,n)=1;
    end
    PAn= PAn+1;
    end
    n=3;
end

% Compare sequence with PA expectation

for n=[2 4];
    for PAn=1:size(PA,1);
    if isequal(seqnonezeros(PAn,1),PA(PAn,1));
        PA(PAn,2)=1;
    else PA(PAn,2)=-1;
    end
    PAn= PAn+1;
    end
    n=4;
end

for PAn=1:size(PA,1);
    if isequal(seqnonezeros(PAn,1),PA(PAn,3));
        PA(PAn,4)=1;
    else PA(PAn,4)=-1;
    end
    PAn= PAn+1;
end

% Score PA expectation

for n= [2 4];
    for PAn = 2:size(PA,1);
    if PA(PAn-1,n)==1 && PA(PAn,n)==1;
        PA(PAn,n)=2;
    elseif PA(PAn-1,n)==2 && PA(PAn,n)==1 ;
        PA(PAn,n)=3;
    elseif PA(PAn-1,n)==-1 && PA(PAn,n)==-1;
        PA(PAn,n)=-2;
    elseif PA(PAn-1,n)==-2 && PA(PAn,n)==-1;
        PA(PAn,n)=-3;
    end
    PAn = PAn+1;
    end
    n = 4;
end


% Position alternation reinforcement score

PAreinf1= sum(PA(:,2))
PAreinf2= sum(PA(:,4))

% Position alternation reinforcement score

PAperf1 = size(find(PA(:,2)>0),1)
PAperf2 = size(find(PA(:,4)>2),1)

% Win-stay, lose-shift
WST = zeros(size(seqnonezeros,1),4);

WST(1,1)=1;
WST(1,3)=2;

% Generate WST expectation
 
for n =[1 3];
    for WSTn = 2:size(WST,1);
    if isequal(WST(WSTn-1,n),seqnonezeros(WSTn-1,1));
        WST(WSTn,n)=WST(WSTn-1,n);
    elseif WST(WSTn-1,n)==1;
        WST(WSTn,n)=2;
    elseif WST(WSTn-1,n)==2
        WST(WSTn,n)=1;
    end
    WSTn= WSTn+1; 
    end
    n = 3;
end   
 
% Compare sequence with WST expectation

for n=[2 4];
    for WSTn=1:size(WST,1);
    if isequal(seqnonezeros(WSTn,1),WST(WSTn,n-1));
        WST(WSTn,n)=1;
    else WST(WSTn,n)=-1;
    end
    WSTn= WSTn+1;
    end
    n=4;
end 

for n= [2 4];
    for WSTn = 2:size(WST,1);
    if WST(WSTn-1,n)==1 && WST(WSTn,n)==1;
        WST(WSTn,n)=2;
    elseif WST(WSTn-1,n)==2 && WST(WSTn,n)==1 ;
        WST(WSTn,n)=3;
    elseif WST(WSTn-1,n)==-1 && WST(WSTn,n)==-1;
        WST(WSTn,n)=-2;
    elseif WST(WSTn-1,n)==-2 && WST(WSTn,n)==-1;
        WST(WSTn,n)=-3;
    end
    WSTn = WSTn+1;
    end
    n = 4;
end


% Win-stay reinforcement score

WSTreinf1= sum(WST(:,2))
WSTreinf2= sum(WST(:,4))

% Win-stay performance score

WSTperf1 = size(find(WST(:,2)>0),1)
WSTperf2 = size(find(WST(:,4)>0),1)

% Win-shift, lose-stay
WSH =zeros(size(seqnonezeros,1),4);

WSH(1,1)=1;
WSH(1,3)=2;

% Generate WSH expectation

for n=[1 3];
    for WSHn = 2:size(WSH,1);
        if isequal(WSH(WSHn-1,n),seqnonezeros(WSHn-1,1));
            if WSH(WSHn-1,1)==1;
               WSH(WSHn,n)=2;
            elseif WSH(WSHn-1,1)==2;
               WSH(WSHn,n)=1;
            end
        else WSH(WSHn,n)=WSH(WSHn-1,n);
        end
        WSHn= WSHn+1;
    end
    n=3;
end

% Score WSH expectation

for n=[2 4];
    for WSHn=1:size(WSH,1);
    if isequal(seqnonezeros(WSHn,1),WSH(WSHn,n-1));
        WSH(WSHn,n)=1;
    else WSH(WSHn,n)=-1;
    end
    WSHn= WSHn+1;
    end
    n=4;
end

for n= [2 4];
    for WSHn = 2:size(WSH,1);
    if WSH(WSHn-1,n)==1 && WSH(WSHn,n)==1;
        WSH(WSHn,n)=2;
    elseif WSH(WSHn-1,n)==2 && WSH(WSHn,n)==1 ;
        WSH(WSHn,n)=3;
    elseif WSH(WSHn-1,n)==-1 && WSH(WSHn,n)==-1;
        WSH(WSHn,n)=-2;
    elseif WSH(WSHn-1,n)==-2 && WSH(WSHn,n)==-1;
        WSH(WSHn,n)=-3;
    end
    WSHn = WSHn+1;
    end
    n = 4;
end 

% Win-shift reinforcement score

WSHreinf1= sum(WSH(:,2))
WSHreinf2= sum(WSH(:,4))

% Win-shift performance score

WSHperf1 = size(find(WSH(:,2)>0),1)
WSHperf2 = size(find(WSH(:,4)>0),1)

    end  

    

listones =find(seq==1);
listtwos =find(seq==2);
 
seq(listones,:)=well1;
seq(listtwos,:)=well2;
 
well1 =size(find(seq==well1),1)
well2 =size(find(seq==well2),1)

datafile = {'PAperf1', PAperf1; 'PAperf2', PAperf2; 'PAreinf1', PAreinf1; 'PAreinf2', PAreinf2; 'PPperf1', PPperf1; 'PPperf2', PPperf2; 'WSTperf1', WSTperf1; 'WSTperf2', WSTperf2; 'WSTreinf1', WSTreinf1; 'WSTreinf2', WSTreinf2; 'WSHperf1', WSHperf1; 'WSHperf2', WSHperf2; 'WSHreinf1', WSHreinf1; 'WSHreinf2', WSHreinf2}; 

stat = [(abs(well1-well2)) (max(abs(PAreinf1),abs(PAreinf2))) (max(abs(PPreinf1),abs(PPreinf2))) (max(abs(WSHreinf1),abs(WSHreinf2))) (max(abs(WSTreinf1),abs(WSTreinf2)))]
temp = [0.05*iter factor*iter factor*iter factor*iter factor*iter]
meanst = mean(stat>temp)
%seq = num2str(seq);
         
cd('/data14/jai/seqfiles/');

filename1 = ([filename,'.txt'])
dataname = ([filename,'_stat_','.txt'])

dlmwrite(filename1,seq,'precision',0);
cell2csv(dataname,datafile);
 
    
cd(currdir);
 
end

