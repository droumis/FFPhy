% Note: Columns correspond to chips.
% The final 3 rows of a column correspond to the triode at the Ground pin.

% Here is an example tetrodemap where all 3 chips have their lowest #
% tetrodes (#1, #8, and #15) begin at the ground pin (i.e. they're triodes)

tetrodemap = zeros(27,2);

% Chip 1 (MS)

tetrodemap(1:4,1) = 1;
tetrodemap(5:8,1) = 2;
tetrodemap(9:12,1) = 3;
tetrodemap(13:16,1) = 4;
tetrodemap(17:20,1) = 5;
tetrodemap(21:24,1) = 6;


% Chip 2

tetrodemap(1:4,2) = 7;
tetrodemap(5:8,2) = 8;
tetrodemap(9:12,2) = 9;
tetrodemap(13:16,2) = 10;
tetrodemap(17:20,2) = 11;
tetrodemap(21:24,2) = 12;
tetrodemap(25:27,2) = 13;



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

GenerateDSPmap(tetrodemap,'master','global','slave1','warming','datadir','/data');
