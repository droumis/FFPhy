

tetrodemap = zeros(27,2);

% Partial

% chip 1
% 
% tetrodemap(1:4,1) = 1;
% tetrodemap(5:8,1) = 2;
% tetrodemap(9:12,1) = 3;
% tetrodemap(13:16,1) = 4;
% tetrodemap(17:20,1) = 5;
% tetrodemap(21:24,1) = 6;
% tetrodemap(25:27,1) = 7;
% 
% %chip 2
% 
% tetrodemap(1:4,2) = 8;
% tetrodemap(5:8,2) = 9;
% tetrodemap(9:12,2) = 10;
% tetrodemap(13:16,2) = 11;
% tetrodemap(17:20,2) = 12;
% tetrodemap(21:24,2) = 13;
% tetrodemap(25:27,2) = 14;

% % chip 3
% 
% tetrodemap(1:4,3) = 15;
% tetrodemap(5:8,3) = 16;
% % tetrodemap(9:12,2) = 17;
% tetrodemap(13:16,3) = 18;
% tetrodemap(17:20,3) = 19;
% tetrodemap(21:24,3) = 20;
% tetrodemap(25:27,3) = 21;
% 
% % chip 4
% 
% tetrodemap(1:4,4) = 22;
% tetrodemap(5:8,4) = 23;
% tetrodemap(9:12,4) = 24;
% tetrodemap(13:16,4) = 25;
% tetrodemap(17:20,4) = 26;
% tetrodemap(21:24,4) = 27;
% tetrodemap(25:27,4) = 28;
% ----------------------------------------------------------------------
% FULL (post adjust)
% chip 1

tetrodemap(1:4,1) = 1;
tetrodemap(5:8,1) = 2;
tetrodemap(9:12,1) = 3;
tetrodemap(13:16,1) = 4;
tetrodemap(17:20,1) = 5;
tetrodemap(21:24,1) = 6;
tetrodemap(25:27,1) = 7;

%chip 2

tetrodemap(1:4,2) = 8;
tetrodemap(5:8,2) = 9;
tetrodemap(9:12,2) = 10;
% % tetrodemap(13:16,2) = 11;
% % tetrodemap(17:20,2) = 12;
tetrodemap(21:24,2) = 13;
tetrodemap(25:27,2) = 14;

% chip 3

% % tetrodemap(1:4,3) = 15;
tetrodemap(5:8,3) = 16;
% % tetrodemap(9:12,2) = 17;
tetrodemap(13:16,3) = 18;
tetrodemap(17:20,3) = 19;
tetrodemap(21:24,3) = 20;
tetrodemap(25:27,3) = 21; %bad but i need another tetrode to make nspike not fuckup

% chip 4

tetrodemap(1:4,4) = 22;
tetrodemap(5:8,4) = 23;
tetrodemap(9:12,4) = 24;
tetrodemap(13:16,4) = 25;
tetrodemap(17:20,4) = 26;
tetrodemap(21:24,4) = 27;
tetrodemap(25:27,4) = 28;

% GenerateDSPmap(tetrodemap,'master','hail','slave1','snow','datadir','/data');
GenerateDSPmap(tetrodemap,'master','snow','slave1','hail','datadir','/data/demetris');




