function out = overlap(vector1, vector2)

%out = overlap(vector1, vector2)
%Computes the overlap between two vectors of equal length.
%Overlap is defined as 2*(overlapping area)/(area1+area2)


if (length(vector1) ~= length(vector2))
    error('The length of the two inputs vectors must be the same.');
end

ratediff = abs(vector1 - vector2);
totalrates = vector1 + vector2;

ratediff(find(isnan(ratediff))) = 0;
totalrates(find(isnan(totalrates))) = 0;
if (sum(totalrates) == 0)
    out = nan;
else
    out = (sum(totalrates)-sum(ratediff))/(sum(totalrates));
end