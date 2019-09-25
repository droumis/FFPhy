function [ output_args ] = twochoicegen(iter,well1,well2,filename)
% Generates a sequence of iter repeats with 0 1 0 2 etc
% Specifies the proportion of 1s and 2s in each sequence
currdir = pwd;

ones = iter/2-5;
twos = iter-ones;

% Set proportion of 

while sqrt((ones-twos)^2) >=0.1*iter;

    seq=[];

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
    
            if isequal(last1,last2,last3);
            
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

num = num +1;

ones =size(find(seq==1),1);
twos =size(find(seq==2),1);
 
end

ones =find(seq==1);
twos =find(seq==2);

seq(ones,:)=well1;
seq(twos,:)=well2;

 
well1 =size(find(seq==well1),1)
well2 =size(find(seq==well2),1)

%seq = num2str(seq);

cd('/data14/jai/seqfiles/');

filename = ([filename,'.txt'])

dlmwrite(filename,seq,'precision',0);
    
cd(currdir);  