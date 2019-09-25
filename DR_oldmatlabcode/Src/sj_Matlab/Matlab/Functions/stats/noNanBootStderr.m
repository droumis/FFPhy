function out = noNanBootStderr(input)
%[meanout, stderrout] = noNanMean(input)
%Computes the mean and stderr of a matrix, ignoring the NaN's

if ~isnumeric(input)
    error('The input must be numeric');
end
input = shiftdim(input);
out = [];
for i = 1:size(input,2)
    tmp = input(find(~isnan(input(:,i))),i);  
    tmp2 = bootstrp(1000,@mean,tmp);
    out(i) = std(tmp2);
end