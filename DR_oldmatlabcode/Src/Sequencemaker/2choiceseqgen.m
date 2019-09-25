function [ output_args ] = twochoicegen(iter)
% Generates a sequence of iter repeats with 0 1 0 2 etc
% 

seq=[];

for num = 1:iter;
    
    currseq =[];

    0 = currseq(1,1);

    randi(2,[1,1])=currseq(2,1);
    
    seq = seq + currseq;
    
    num = num +1;
    
end

    
