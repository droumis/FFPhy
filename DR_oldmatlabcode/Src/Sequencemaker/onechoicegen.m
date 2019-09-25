function [ output_args ] = onechoicegen(iter,well,filename)
% Generates a sequence of iter repeats.
% Well specifies the number repeated, 1, 2, 3 etc.
% Specifies the proportion of 1s and 2s in each sequence
currdir = pwd;



seq=[];

for num = 1:iter;
    
currseq =[];

currseq(1,1)=0;
       
currseq(2,1)=well;
    
seq = [seq; currseq];
    
end 

num = num +1;

cd('/data14/jai/seqfiles/');

filename = ([filename,'.txt'])

dlmwrite(filename,seq,'precision',0);
    
cd(currdir);  