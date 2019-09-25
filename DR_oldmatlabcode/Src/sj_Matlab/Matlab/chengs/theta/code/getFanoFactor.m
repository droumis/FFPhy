function F=getFanoFactor(counts)
% function F=getFanoFactor(c)
%
% calculate Fano Factor from 9 adjacent bins 
% periodic boundary conditions in first dimension (theta)

n= size(counts);
F= zeros(n);

% periodic boundary conditions in y (1st dim)
c= zeros(n(1)+2, n(2));
c(2:end-1,:)= counts;
c(1,:)= counts(end,:);
c(end,:)= counts(1,:);

% first and last bin in x don't have 9 adjacent bins
F(:,1)= nan; F(:,end)= nan;

for i=2:n(1)+1
    for j=2:n(2)-1
        block= reshape(c(i-1:i+1, j-1:j+1),9,1);
        if mean(block)==0
            F(i-1,j)= 1;
        else
            F(i-1,j)= var(block)/mean(block);
        end
    end
end

