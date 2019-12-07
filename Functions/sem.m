function out = sem(x,dim)
% out = sem(x,dim)
% Calculate standard error of the mean using n-1
% x is the data array
% dim is the dimension over which to calculate std and n (1 = rows, 2 = columns).

out = std(x,0,dim)/sqrt(size(x,dim)-1);

end