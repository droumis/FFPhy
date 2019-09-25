
%for adjusting, i'm making an mec/por config and an HC config.. later I'll combine into 1 config


tetrodemap = zeros(27,2);

%MEC and POR

% tetrodemap(1:4,1) = 8;
% tetrodemap(5:8,1) = 9; %por
% tetrodemap(9:12,1) = 10; %por
% tetrodemap(13:16,1) = 11; %por
% tetrodemap(17:20,1) = 12; %por
% tetrodemap(21:24,1) = 13; %por
% tetrodemap(25:27,1) = 14; %por
% 
% tetrodemap(1:4,2) = 1; %por
% tetrodemap(5:8,2) = 2;
% tetrodemap(9:12,2) = 3;
% tetrodemap(13:16,2) = 4;
% tetrodemap(17:20,2) = 5;
% tetrodemap(21:24,2) = 6;
% tetrodemap(25:27,2) = 7;

% HC

tetrodemap(1:4,3) = 22;
tetrodemap(5:8,3) = 23;
tetrodemap(9:12,3) = 24;
tetrodemap(13:16,3) = 25;
tetrodemap(17:20,3) = 26;
tetrodemap(21:24,3) = 27;
tetrodemap(25:27,3) = 28;

tetrodemap(1:4,4) = 15;
tetrodemap(5:8,4) = 16;
tetrodemap(9:12,4) = 17;
tetrodemap(13:16,4) = 18;
tetrodemap(17:20,4) = 19;
tetrodemap(21:24,4) = 20;
tetrodemap(25:27,4) = 21;


% GenerateDSPmap(tetrodemap,'master','hail','slave1','snow','datadir','/data');
GenerateDSPmap(tetrodemap,'master','snow','slave1','hail','datadir','/data/demetris');

% 
% % LLBL 2 mm LFP probe   kk ea jc 6.4.13
% 
% tetrodemap(1,1) = 37;
% tetrodemap(2,1) = 38;
% tetrodemap(3,1) = 3;
% tetrodemap(4,1) = 2;
% tetrodemap(5,1) = 1;
% tetrodemap(6,1) = 6;
% tetrodemap(7,1) = 5;
% tetrodemap(8,1) = 4;
% tetrodemap(9,1) = 9;
% tetrodemap(10,1) = 8;
% tetrodemap(11,1) = 7;
% tetrodemap(12,1) = 12;
% tetrodemap(13,1) = 11;
% tetrodemap(14,1) = 10;
% tetrodemap(15,1) = 15;
% tetrodemap(16,1) = 14;
% tetrodemap(17,1) = 13;
% tetrodemap(18,1) = 18;
% tetrodemap(19,1) = 17;
% tetrodemap(20,1) = 16;
% tetrodemap(21,1) = 39;
% tetrodemap(22,1) = 40;
% tetrodemap(23,1) = 41;
% tetrodemap(24,1) = 42;
% tetrodemap(25,1) = 43;
% tetrodemap(26,1) = 44;
% tetrodemap(27,1) = 45;
% 
% tetrodemap(1,2) = 46;
% tetrodemap(2,2) = 47;
% tetrodemap(3,2) = 20;
% tetrodemap(4,2) = 21;
% tetrodemap(5,2) = 24;
% tetrodemap(6,2) = 19;
% tetrodemap(7,2) = 22;
% tetrodemap(8,2) = 23;
% tetrodemap(9,2) = 26;
% tetrodemap(10,2) = 27;
% tetrodemap(11,2) = 30;
% tetrodemap(12,2) = 25;
% tetrodemap(13,2) = 28;
% tetrodemap(14,2) = 29;
% tetrodemap(15,2) = 32;
% tetrodemap(16,2) = 33;
% tetrodemap(17,2) = 36;
% tetrodemap(18,2) = 31;
% tetrodemap(19,2) = 48;
% tetrodemap(20,2) = 49;
% tetrodemap(21,2) = 50;
% tetrodemap(22,2) = 51;
% tetrodemap(23,2) = 52;
% tetrodemap(24,2) = 53;
% tetrodemap(25,2) = 54;
% tetrodemap(26,2) = 55;
% tetrodemap(27,2) = 56;
% 
% GenerateDSPmap(tetrodemap,'eegexclude',37:56,'master','global','slave1','warming','datadir','/data/emily')
% 
