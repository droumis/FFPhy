% to gain intuition about  mutual information as measure of correlation
% mutual information

addpath('../aux');

disp(['mutual information vs. linear correlation coefficient r for Gaussian' ...
       ' PDF ...'])
MutualInfoTest
pause

disp(['mutual information as function of ''bin width''...'])
MutualInfoTest2
pause

disp(['mutual information vs. linear correlation coefficient for simple' ...
      ' distributions  ...'])
MutualInfoTest3

