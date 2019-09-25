% This script quantifies behavior from dio data. First, diodayprocess was
% run on all days (see preprocessing file).
% The way I quantify is as follows:
% For each reward given, I look at the time window starting from the 
% previous reward to just before the current reward was given (1000 samples 
% back from the current reward, so as not to include the poke invoking the new reward). 
% If during that interval there were only pokings corresponding to the previous 
% reward well, this is a correct trial (there are always multiple pokings to 
% get a single reward). If during that interval other wells were triggered, 
% this is an error trial. Errors and correct trials are stored per pair of
% wells in the following form: an error trial is stored in errorTypes(x,y),
% where x is the previous reward location, and y is the first well he
% erroneously poked in after that (notice, NOT where he was supposed to
% go). A correct trial is counted in correctTypes(x,y) where x is the previous 
% rewarded well and y is the current rewarded well.


% 11,13,16 are inputs (rat triggering wells)
% 30,31,32 are outputs (milk coming out)
% On day 10 no reward was given in sound well, so no bit 30

toSilentCorrectRatios=[];
toSoundCorrectRatios=[];
toHomeCorrectRatios=[];
allRewards=[];
sequentialTrials=[];

for day1=1:11
    if day1<10
        dayStr=['0' num2str(day1)];
    else
        dayStr=[num2str(day1)];
    end
        load(['/data15/gideon/Brg/brgDIO' dayStr '.mat'])

    d11=diopulses{11}.pulsetimes(:,1);
    d16=diopulses{16}.pulsetimes(:,1);
    d13=diopulses{13}.pulsetimes(:,1);
    try
        d30=diopulses{30}.pulsetimes(:,1);
    catch
        'No sound rewards, should be only day 10'
        d30=0;
    end
    d31=diopulses{31}.pulsetimes(:,1);
    d32=diopulses{32}.pulsetimes(:,1);
    % The correspondence changes by day
    
    % day 8:
    % 11-32 silent well
    % 13-31 sound well
    % 16-30 home well
  
        silentWell=d11;
        silentReward=d31;
        soundWell=d13;
        soundReward=d30;
        homeWell=d16;
        homeReward=d32;
        corresp1=[11 31;13 30;16 32];
        
    
    
    allWells={homeWell,silentWell,soundWell};
    
    % Deleting spurious rewards that sometimes appear as doubles
    soundReward=[soundReward(diff(soundReward)>5000); soundReward(end)];
    silentReward=[silentReward(diff(silentReward)>5000); silentReward(end)];
    homeReward=[homeReward(diff(homeReward)>5000); homeReward(end)];
    
    
    figure;
    h1=plot(homeWell,1,'bx');
    hold on
    plot(homeReward,1.1,'bo','markerfacecolor','b')
    
    h2=plot(silentWell,1.5,'rx');
    plot(silentReward,1.6,'ro','markerfacecolor','r')
    
    h3=plot(soundWell,2,'kx');
    plot(soundReward,2.1,'ko','markerfacecolor','k')
    
    errorTypes=zeros(3,3);
    correctTypes=zeros(3,3);
    
    [firstRew firstGivenIn]=min([homeReward(1),silentReward(1),soundReward(1)]);
    prevRew=firstRew;
    prevGivenIn=firstGivenIn;
    while true
        [nextRew nextGivenIn]=min([homeReward(find(homeReward>prevRew,1)),silentReward(find(silentReward>prevRew,1)),soundReward(find(soundReward>prevRew,1))]);
        notGivenIn=setdiff([1 2 3],prevGivenIn);
        
        % note that the error is written down as going from the prev reward to
        % where the rat actually poked, not where he was supposed to.
        
        % Only poking in the first unrewarded well during the interval
        if ~isempty(find(allWells{notGivenIn(1)}>prevRew & allWells{notGivenIn(1)}<nextRew-1000)) && ...
                isempty(find(allWells{notGivenIn(2)}>prevRew & allWells{notGivenIn(2)}<nextRew-1000))
                        errorTypes(prevGivenIn,notGivenIn(1))=errorTypes(prevGivenIn,notGivenIn(1))+1;
        sequentialTrials=[sequentialTrials 0];
            % Only poking in the second unrewarded well during the interval
        elseif ~isempty(find(allWells{notGivenIn(2)}>prevRew & allWells{notGivenIn(2)}<nextRew-1000)) && ...
                isempty(find(allWells{notGivenIn(1)}>prevRew & allWells{notGivenIn(1)}<nextRew-1000))
                      errorTypes(prevGivenIn,notGivenIn(2))=errorTypes(prevGivenIn,notGivenIn(2))+1;
                    sequentialTrials=[sequentialTrials 0];

            % if poking in two unrewarded wells, attribute the error to the first one
            % (for example, goes wrongly to the
            % silent arm, then tries the sound arm, and only then returning)
        elseif ~isempty(find(allWells{notGivenIn(2)}>prevRew & allWells{notGivenIn(2)}<nextRew-1000)) && ...
                ~isempty(find(allWells{notGivenIn(1)}>prevRew & allWells{notGivenIn(1)}<nextRew-1000))
            if min( allWells{notGivenIn(2)}(find(allWells{notGivenIn(2)}>prevRew & allWells{notGivenIn(2)}<nextRew-1000)))<...
                    min( allWells{notGivenIn(1)}(find(allWells{notGivenIn(1)}>prevRew & allWells{notGivenIn(1)}<nextRew-1000)))
                errorTypes(prevGivenIn,notGivenIn(2))=errorTypes(prevGivenIn,notGivenIn(2))+1;
            else
                errorTypes(prevGivenIn,notGivenIn(1))=errorTypes(prevGivenIn,notGivenIn(1))+1;
            end
                    sequentialTrials=[sequentialTrials 0];

                        % correct trial
        else
             correctTypes(prevGivenIn,nextGivenIn)=correctTypes(prevGivenIn,nextGivenIn)+1;
                sequentialTrials=[sequentialTrials 1];

        end
                
        prevRew=nextRew;
        prevGivenIn=nextGivenIn;
        if (isempty(find(homeReward>prevRew,1)) & isempty(find(silentReward>prevRew,1)) & isempty(find(soundReward>prevRew,1)))
            break
        end
        
    end
%    figure;imagesc(errorTypes)
%    figure;imagesc(correctTypes)
    
    % note, errorTypes(1,3) is the number of cases where he wronly went to 3,
    % in other words, when he should have gone to 2, so the correct ratio of
    % going to 2 is the amount he actually went there (correctTypes(1,2))
    % divided by the sum of that and errorTypes(1,3). Symmetrically for second
    % case.
    toSilentCorrectRatio=correctTypes(1,2)/(correctTypes(1,2)+errorTypes(1,3));
    toSoundCorrectRatio=correctTypes(1,3)/(correctTypes(1,3)+errorTypes(1,2));
    toHomeCorrectRatio=(correctTypes(2,1)+correctTypes(3,1))/(correctTypes(2,1)+correctTypes(3,1)+errorTypes(2,3)+errorTypes(3,2))
    
    toSilentCorrectRatios=[toSilentCorrectRatios toSilentCorrectRatio];
    toSoundCorrectRatios=[toSoundCorrectRatios toSoundCorrectRatio];
    toHomeCorrectRatios=[toHomeCorrectRatios toHomeCorrectRatio];
    allRewards=[allRewards length(homeReward)+length(silentReward)+length(soundReward)];
end
figure;
subplot(2,1,1)
imagesc(sequentialTrials)
title('correct/incorrect trials throughout experiment')
subplot(2,1,2)
plot(smooth(sequentialTrials,150),'k','linewidth',2)
xlabel('# trial')
ylabel('% correct, smoothed')
xlim([0 length(sequentialTrials)])


figure;
plot(allRewards,'r','linewidth',2)
xlabel('Day')
ylabel('# Trials')
axis([0 12 0 300])

figure;
plot(toSoundCorrectRatios,'linewidth',2)
hold on
plot(toSilentCorrectRatios,'r','linewidth',2)
plot((toSilentCorrectRatios+toSoundCorrectRatios)./2,'k','linewidth',2)
plot(toHomeCorrectRatios,'m','linewidth',2)
title('Behavior performance over days')
ylabel('Percent correct')
xlabel('Day')
legend('To Sound','To Silent','Sound/Silent average','To home')
xlim([0 12])
