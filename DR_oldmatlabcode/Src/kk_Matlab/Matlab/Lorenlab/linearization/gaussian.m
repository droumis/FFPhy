% g = GAUSSIAN(stdev, size)
%     creates a 1D gaussian matrix with SIZE elements with the given stdev
%     centered in the middle of g

function [g] = gaussian(stdev, size)

size = round(size);
g = zeros(size,1);

m = 1 + (size - 1) / 2;
mean = [m];


for i = 1:size
	x = i;
	g(i) = 1 / (stdev * sqrt(2 * pi)) * exp(-.5 * (mean - x)^2 / stdev^2);
end

% renormalize 
g = g ./ sum(g);
