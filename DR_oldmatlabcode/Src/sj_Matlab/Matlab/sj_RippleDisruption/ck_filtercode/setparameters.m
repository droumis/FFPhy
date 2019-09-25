
fprintf('Setting parameters');

pulsewin = [-0.15, 0.25];

FS = 30000;

% {[Epochs] [CA1 CA3]}
animalTetrodeTable = {...
  'Neb',          {[1:10],[12 13]}};  % 10 has the largest pop spike

animals = {'Neb'}; 

% Normalization parameters
movingThresh = 10; % cm/s
quieThresh = 5; % still for at least 5 s

% normType = 'z-score';
normType = 'average';
normCondition = sprintf('> %f',movingThresh);

% Helpful for plotting
%-----------------------------------------------------
twin = (round(pulsewin(1)*FS):round(pulsewin(2)*FS))/FS;
twin = twin(1:end-1);
[closetozero,t0] = min(abs(twin));
%-----------------------------------------------------

CheckSignificance = 1;
BalanceVels = 0;
BootstrapSize = 1000;
DoMedians = 1;

% For plotting vs. velocity
% C = 20; vMax = 40; velCats = linspace(0,vMax,C+1);
velCats = [0.125 0.25 0.5 1 2 4 8 16 32];
veltick = [0.125 0.25 0.5 1 2 4 8 16 32]; 
minV = 0.125;
maxV = 32;


fprintf('\n');
