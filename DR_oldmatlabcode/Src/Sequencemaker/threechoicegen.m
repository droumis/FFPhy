function [ output_args ] = threechoicegen(iter,filename)
% Generates a sequence of iter repeats with 0 1 0 2 0 3 etc
% 
currdir = pwd;

seq=[];

for num = 1:iter;
    
    currseq =[];

    currseq(1,1)=0;

    r=randi(3,[1,1]);
    
    currseq(2,1)=r;
    
    seq = [seq; currseq];
    
    num = num +1;
    
end

%seq = num2str(seq);

cd('/data14/jai/seqfiles/');

filename = ([filename,'.txt'])

dlmwrite(filename,seq,'precision',0);
    
cd(currdir); 