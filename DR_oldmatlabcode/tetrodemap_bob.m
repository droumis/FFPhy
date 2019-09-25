% Note: Columns correspond to chips.
% The final 3 rows of a column correspond to the triode at the Ground pin.

% Here is an example tetrodemap where all 3 chips have their lowest #
% tetrodes (#1, #8, and #15) begin at the ground pin (i.e. they're triodes)

tetrodemap = zeros(27,2);

% Chip 1 hippocampus

tetrodemap([23 25 27],2) = 14;
tetrodemap([15 17 19 21],2) = 13;
tetrodemap([7 9 11 13],2) = 12;
tetrodemap([3 5 24 26],2) = 11;
tetrodemap([16 18 20 22],2) = 10;
tetrodemap([8 10 12 14],2) = 9;
tetrodemap([2 4 6],2) = 8;


% Chip 2 medial septum

tetrodemap([23 25 27],1) = 1;
tetrodemap([15 17 19 21],1) = 2;
tetrodemap([7 9 11 13],1) = 3;
tetrodemap([3 5 24 26],1) = 4;
tetrodemap([16 18 20 22],1) = 5;
tetrodemap([8 10 12 14],1) = 6;
tetrodemap([2 4 6],1) = 7;


% % Chip 3
% 
% tetrodemap(1:4,3) = 21;
% tetrodemap(5:8,3) = 20;
% tetrodemap(9:12,3) = 19;
% tetrodemap(13:16,3) = 18;
% tetrodemap(17:20,3) = 17;
% tetrodemap(21:24,3) = 16;
% tetrodemap(25:27,3) = 15;

% And this is how to call GenerateDSPmap.

GenerateDSPmap(tetrodemap,'master','hot','slave1','humid','datadir','/data/kkay');
