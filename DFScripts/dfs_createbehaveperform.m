% get the behaveperform (list of 1's and 0's for correct and incorrect
% trials for each epoch) and run state space analysis
%-----------------------------------------------------

savefile = 0;

Ltype = 'S1Acquisition'; %specify learning type: 'S1Acquisition' or 'Switch' - if running a W track animal, choose 'S1Acquisition'
runstatespace = 1; %set to 1 to run state space analysis
% if runscript
    
    % Animal Selection
    animals = {'Bond'}; {'Chapati'}; {'Higgs'}; {'Government'};
    
    % Day Filter
    dayfilter = 10; %Fabio days 6-13 for S1 acquisition, 14-21 for switch
    
    % Epoch Filter
    epochfilter = 'isequal($type, ''run'') && (isequal($environment, ''WTrackA''))'; %'isequal($type, ''run'') && isequal($environment, ''WTrackA'') && (($exposure>=1) && ($exposure<=10))';  %

iterator = 'singleepochanal';

timefilter = {}; % {'getlinvelocity', ['((abs($velocity) > -1))']} };

F = createfilter('animal', animals,'epochs', epochfilter, 'excludetime', timefilter, 'iterator', iterator);
% F = createfilter('animal', animals,'epochs', epochfilter, 'excludetime', timefilter, 'iterator', iterator);

if strcmp(Ltype,'S1Acquisition')
    F = setfilterfunction(F, 'dfa_createbehaveperform', {'rewardinfo','task'},'learningtype',Ltype); 
elseif strcmp(Ltype,'Switch')
    F = setfilterfunction(F, 'dfa_createbehaveperform_switch', {'rewardinfo','oppseqrewardinfo','task'},'learningtype',Ltype); % Include 'oppseqrewardinfo' as an input variable if Ltype == Switch
end
F = runfilter(F);

%% Plotting options
combineplots = 0; %combine Acquisition and Switch plots

%% Concatenate epochs into behaveperform

if length(animals)==1
    f = F.output{1};
else
    f = struct([]);
    for an=1:length(animals)
    f = ([f F(an).output{1}]);
    end
end

%% for S1Acquisition

if strcmp(Ltype,'S1Acquisition') 
%initialize behaveperform
behaveperform = struct;
behaveperform.outbound = [];
behaveperform.outtrials = 0;
behaveperform.inbound = [];
behaveperform.intrials = 0;
behaveperform.all = [];
behaveperform.alltrials = 0;
% Collect epoch boundaries
AcEpochEdges = zeros(length(f),1);

for ep = 1:length(f)
    if ep == 1
        AcEpochEdges(ep) = numel(f(ep).combinedperf);
    else
         AcEpochEdges(ep) = AcEpochEdges(ep-1)+numel(f(ep).combinedperf);
    end
    behaveperform.outbound = [behaveperform.outbound; f(ep).outboundperf];
    behaveperform.outtrials = [behaveperform.outtrials + length(f(ep).outboundperf)];
    behaveperform.inbound = [behaveperform.inbound; f(ep).inboundperf];
    behaveperform.intrials = [behaveperform.intrials + length(f(ep).inboundperf)];
    behaveperform.all = [behaveperform.all; f(ep).combinedperf];
    behaveperform.alltrials = [behaveperform.alltrials + length(f(ep).combinedperf)];
    behaveperform.learningtype = Ltype;   
end

%save 
if savefile
filename = sprintf('%s%sbehaveperform_%s.mat',F.animal{2},F.animal{3},Ltype);
save(filename,'behaveperform');
end

% run state state algorithm
if runstatespace
    if strcmp(animals,'Eliot') || strcmp(animals,'Fabio')
        chance = 0.166667;
    elseif strcmp(animals,'Egypt') || strcmp(animals,'Government') || strcmp(animals,'Higgs') || strcmp(animals,'Chapati')
        chance = 0.5;
    end
%     [pcAcOut, ltAcOut] = getestprobcorrect(behaveperform.outbound,chance,1);
%     [pcAcIn, ltAcIn] = getestprobcorrect(behaveperform.inbound,chance,1);
    [pcAcAll, ltAcAll] = getestprobcorrect(behaveperform.all,chance,1,0);  %setting plot to 0 to hide individual days
end

% plot behavior curve
Actrials = [1:(size(pcAcAll,1)-1)]';

figure
hold on
plot(Actrials,pcAcAll(2:end,1),'-r','LineWidth',2)
plot(Actrials,pcAcAll(2:end,2),'--r','LineWidth',1)
plot(Actrials,pcAcAll(2:end,3),'--r','LineWidth',1)
plot(Actrials,repmat(chance,length(Actrials),1),'--','Color',[0.3 0.3 0.3],'LineWidth',1.5)
for i = 1:length(AcEpochEdges)
line(AcEpochEdges(i),[0:0.001:1],'Color',[0.5 0.5 0.5])
end
% line(AcEpochEdges(end),[0:0.001:1],'Color','k')
xlabel('Trials')
ylabel('Probability of a Correct Response')
end

%% for Switch

if strcmp(Ltype,'Switch')
%initialize behaveperform
behaveperform = struct;
behaveperform.outbound = cell(1,length(f));
behaveperform.outtrials =0;
behaveperform.inbound = cell(1,length(f));
behaveperform.intrials = 0;
behaveperform.all = cell(1,length(f));
behaveperform.alltrials = 0;

   for ep = 1:length(f)
    if strcmp(f(ep).sequence,'S1')
    behaveperform.outbound{ep}{1} = f(ep).outboundperf; % S1 is rewarded performance
    behaveperform.inbound{ep}{1} = f(ep).inboundperf; % S1 is rewarded performance
    behaveperform.all{ep}{1} = f(ep).combinedperf; % S1 is rewarded performance
    behaveperform.outbound{ep}{2} = f(ep).OPPoutboundperf; % S2 is opposite, unrewarded performance
    behaveperform.inbound{ep}{2} = f(ep).OPPinboundperf; % S2 is opposite, unrewarded performance
    behaveperform.all{ep}{2} = f(ep).OPPcombinedperf; % S2 is opposite, unrewarded performance
    elseif strcmp(f(ep).sequence,'S2')
    behaveperform.outbound{ep}{1} = f(ep).OPPoutboundperf; % S1 is opposite, unrewarded performance 
    behaveperform.inbound{ep}{1} = f(ep).OPPinboundperf; % S1 is opposite, unrewarded performance 
    behaveperform.all{ep}{1} = f(ep).OPPcombinedperf; % S1 is opposite, unrewarded performance 
    behaveperform.outbound{ep}{2} = f(ep).outboundperf; % S2 is rewarded performance
    behaveperform.inbound{ep}{2} = f(ep).inboundperf; % S2 is rewarded performance
    behaveperform.all{ep}{2} = f(ep).combinedperf; % S2 is rewarded performance
    end
    behaveperform.outtrials = [behaveperform.outtrials + length(f(ep).outboundperf)];
    behaveperform.intrials = [behaveperform.intrials + length(f(ep).inboundperf)];
    behaveperform.alltrials = [behaveperform.alltrials + length(f(ep).combinedperf)];
   end

   %save 
   if savefile
       filename = sprintf('%s/%sbehaveperform_%s.mat',F.animal{2},F.animal{3},Ltype);
       save(filename,'behaveperform');
   end

% Run state space algorithm
if runstatespace
    S1outProbcorrect = zeros(behaveperform.outtrials(1),3);
    S2outProbcorrect = zeros(behaveperform.outtrials(1),3);
    S1inProbcorrect = zeros(behaveperform.intrials(1),3);
    S2inProbcorrect = zeros(behaveperform.intrials(1),3);
    S1allProbcorrect = zeros(behaveperform.alltrials(1),3);
    S2allProbcorrect = zeros(behaveperform.alltrials(1),3);
    backprob1 = pcAcAll(end,1);
    chance = 0.166667;
    startindS1 = 1;
    startindS2 = 1;
    endindS1 = 0;
    endindS2 = 0;
    for ep = 1:length(behaveperform.outbound)
        if ep == 1
            %        [pcS1, ltS1] = getestprobcorrect(behaveperform.outbound{ep}{1},backprob1,1);
            %        [pcS2, ltS2] = getestprobcorrect(behaveperform.outbound{ep}{2},chance,1);
            [pcS1, ltS1] = getestprobcorrect(behaveperform.all{ep}{1},backprob1,1,0); %setting plot to 0 to hide individual days
            [pcS2, ltS2] = getestprobcorrect(behaveperform.all{ep}{2},chance,1,0);
        else
            %               [pcS1, ltS1] = getestprobcorrect(behaveperform.outbound{ep}{1},backprobS1,1);
            %        [pcS2, ltS2] = getestprobcorrect(behaveperform.outbound{ep}{2},backprobS2,1);
            [pcS1, ltS1] = getestprobcorrect(behaveperform.all{ep}{1},backprobS1,1,0);
            [pcS2, ltS2] = getestprobcorrect(behaveperform.all{ep}{2},backprobS2,1,0);
        end
        %        endindS1 = endindS1 + numel(behaveperform.outbound{ep}{1});
        %        endindS2 = endindS2 + numel(behaveperform.outbound{ep}{2});
        endindS1 = endindS1 + numel(behaveperform.all{ep}{1});
        endindS2 = endindS2 + numel(behaveperform.all{ep}{2});
        %            S1outProbcorrect(startindS1:endindS1,:) = pcS1(2:end,:);
        %            S2outProbcorrect(startindS2:endindS2,:) = pcS2(2:end,:);
        S1allProbcorrect(startindS1:endindS1,:) = pcS1(2:end,:);
        S2allProbcorrect(startindS2:endindS2,:) = pcS2(2:end,:);
        backprobS1 = pcS1(end,1); %background prob for the next epoch = mode of probability at end of this epoch
        backprobS2 = pcS2(end,1);
        startindS1 = endindS1+1;
        startindS2 = endindS2+1;
    end
end

% Collect epoch boundaries
SwEpochEdges = zeros(length(behaveperform.all),1);
for i= 1:length(behaveperform.all)
    if i==1
    SwEpochEdges(i) = numel(behaveperform.all{i}{1});
    else
    SwEpochEdges(i) = SwEpochEdges(i-1)+numel(behaveperform.all{i}{1});
    end
end

% Plot behavior curves

Swtrials = [1:size(S1allProbcorrect,1)]';

%plot only first 8 switch epochs for clarity
figure
hold on
plot(Swtrials(1:SwEpochEdges(16)),S1allProbcorrect(1:SwEpochEdges(16),1),'-r','LineWidth',2)
plot(Swtrials(1:SwEpochEdges(16)),S1allProbcorrect(1:SwEpochEdges(16),2),'--r','LineWidth',1)
plot(Swtrials(1:SwEpochEdges(16)),S1allProbcorrect(1:SwEpochEdges(16),3),'--r','LineWidth',1)
plot(Swtrials(1:SwEpochEdges(16)),S2allProbcorrect(1:SwEpochEdges(16),1),'-b','LineWidth',2)
plot(Swtrials(1:SwEpochEdges(16)),S2allProbcorrect(1:SwEpochEdges(16),2),'--b','LineWidth',1)
plot(Swtrials(1:SwEpochEdges(16)),S2allProbcorrect(1:SwEpochEdges(16),3),'--b','LineWidth',1)
plot(Swtrials(1:SwEpochEdges(16)),repmat(chance,(SwEpochEdges(16)),1),'Color',[0.5 0.5 0.5])
for i = 1:16
line(SwEpochEdges(i),[0:0.001:1],'Color','k')
end
xlabel('Trials')
ylabel('Probability of a Correct Response')

%plot all switch epochs
figure
hold on
plot(Swtrials,S1allProbcorrect(:,1),'-r','LineWidth',2)
plot(Swtrials,S1allProbcorrect(:,2),'--r','LineWidth',1)
plot(Swtrials,S1allProbcorrect(:,3),'--r','LineWidth',1)
plot(Swtrials,S2allProbcorrect(:,1),'-b','LineWidth',2)
plot(Swtrials,S2allProbcorrect(:,2),'--b','LineWidth',1)
plot(Swtrials,S2allProbcorrect(:,3),'--b','LineWidth',1)
plot(Swtrials,repmat(chance,length(Swtrials),1),'Color',[0.5 0.5 0.5])
for i = 1:length(SwEpochEdges)
line(SwEpochEdges(i),[0:0.001:1],'Color','k')
end

end

%% Plot S1Acquisition and Switch together

if combineplots
    newswitchtrials = Swtrials+Actrials(end);
    newSwEpochEdges = SwEpochEdges+Actrials(end);
    Swstart = Actrials(end)+1;
    alltrials = [Actrials; newswitchtrials];

% plot all acquisitoin + 8 switch epochs
figure
hold on
plot(Actrials,pcAcAll(2:end,1),'-b','LineWidth',2)
plot(Actrials,pcAcAll(2:end,2),'--b','LineWidth',1)
plot(Actrials,pcAcAll(2:end,3),'--b','LineWidth',1)
plot(alltrials(1:newSwEpochEdges(8)),repmat(chance,newSwEpochEdges(8),1),'--','Color',[0.3 0.3 0.3],'LineWidth',1.5) %chance line
% for i = 1:length(AcEpochEdges)
% line(AcEpochEdges(i),[0:0.001:1],'Color',[0.5 0.5 0.5])
% end
line(AcEpochEdges(end),[0:0.001:1],'Color','k')
plot(Swstart:newSwEpochEdges(8),S1allProbcorrect(1:SwEpochEdges(8),1),'-b','LineWidth',2)
plot(Swstart:newSwEpochEdges(8),S1allProbcorrect(1:SwEpochEdges(8),2),'--b','LineWidth',1)
plot(Swstart:newSwEpochEdges(8),S1allProbcorrect(1:SwEpochEdges(8),3),'--b','LineWidth',1)
plot(Swstart:newSwEpochEdges(8),S2allProbcorrect(1:SwEpochEdges(8),1),'-r','LineWidth',2)
plot(Swstart:newSwEpochEdges(8),S2allProbcorrect(1:SwEpochEdges(8),2),'--r','LineWidth',1)
plot(Swstart:newSwEpochEdges(8),S2allProbcorrect(1:SwEpochEdges(8),3),'--r','LineWidth',1)
for i = 1:8 %length(SwEpochEdges)
    if i == 3 || i==6 
line(newSwEpochEdges(i),[0:0.001:1],'Color','k')
    else
        line(newSwEpochEdges(i),[0:0.001:1],'Color',[0.5 0.5 0.5])
    end
end
xlabel('Trials')
ylabel('Probability of a Correct Response')

% plot Switch behavior alone with bumped trial numbers (to continue trial
% count from end of acquisition)
figure 
hold on
plot(Swstart:newSwEpochEdges(8),repmat(chance,length(Swstart:newSwEpochEdges(8)),1),'--','Color',[0.3 0.3 0.3],'LineWidth',1.5); %chance line
plot(Swstart:newSwEpochEdges(8),S1allProbcorrect(1:SwEpochEdges(8),1),'-b','LineWidth',2)
plot(Swstart:newSwEpochEdges(8),S1allProbcorrect(1:SwEpochEdges(8),2),'--b','LineWidth',1)
plot(Swstart:newSwEpochEdges(8),S1allProbcorrect(1:SwEpochEdges(8),3),'--b','LineWidth',1)
plot(Swstart:newSwEpochEdges(8),S2allProbcorrect(1:SwEpochEdges(8),1),'-r','LineWidth',2)
plot(Swstart:newSwEpochEdges(8),S2allProbcorrect(1:SwEpochEdges(8),2),'--r','LineWidth',1)
plot(Swstart:newSwEpochEdges(8),S2allProbcorrect(1:SwEpochEdges(8),3),'--r','LineWidth',1)
for i = 1:8 %length(SwEpochEdges)
    if i == 3 || i==6 
line(newSwEpochEdges(i),[0:0.001:1],'Color','k')
    else
        line(newSwEpochEdges(i),[0:0.001:1],'Color',[0.5 0.5 0.5])
    end
end

end

