% [pval] = WATSONWILLIAMS(a1, a2)
%         Given two distributions on a circle, uses the Watson-Williams F test to
%         determine the significance of the difference of the means.
%         The test requires a sufficiently tight distribution of angles and returns -1 if
%         that condition is not met.
%         Note that it is important to examine the distributions by eye as well.
%         From _Circular Statistics in Biology_ by Edward Batshelet, 1981, Academic Press
function [pval] = watsonwilliams(a1, a2)

% load the table of the estimates of kappa
load /home/loren/matlab/stats/kappa.mat;

% get the values for R1 and R2
c1 = sum(cos(a1));
s1 = sum(sin(a1));
c2 = sum(cos(a2));
s2 = sum(sin(a2));
R1 = sqrt(c1^2 + s1^2);
R2 = sqrt(c2^2 + s2^2);
R = sqrt((c1 + c2)^2 + (s1 + s2)^2);
n = length(a1) + length(a2);
  
rtilde = (R1 + R2) / n;
ntilde = n / 2;

% get the value for kappa
kval = interp2(kappar, kappan, kappa, rtilde, ntilde);

%keyboard
if (kval < 2)
	% the assuptions of the test are not met, so set sig to -1
	pval = -1;
	return
end

g = 1 + 3 / (8 * kval);

F = g * (n - 2) * (R1 + R2 - R) / (n - (R1 + R2));

pval = 1 - fcdf(F, 1, n-2)
