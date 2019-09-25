%Create sequences to use for TI behavior tasks and save them as text files

%define the pairs you want to create.. 0 is always home and will always lead and interlace
%[ab bc de ef cd bd be ce af] OLD
%[ 1   2    3    4    5    6   7    8   9  ] OLD

%[ab bc de ef cd be af bd ce] CURRENT
%[ 1   2    3    4    5    6   7    8   9  ]

clear all; close all;
% pairs = [1 5]; %
% pairsweighted = [1 2 3 3 4 6 6 7];
pairsweighted = [1 2 3 3 4 5 6 6];
% pairsetname = 'ph3setsrnd'; %text file save name
% pairsetname = 'ph3_CDpair'; %text file save name
% pairsetname = 'ph4_allrnd'; %text file save name
% pairsetname = 'ph4_1-5_w5x4'; %text file save name %weighted 1 is weighting all new pairs x2
pairsetname = 'ph4_1-6_w3-6x2'; %text file save name %weighted 1 is weighting all new pairs x2
animal = 'T09'; %T000 for testing
lengthlist = 400; %1k is prob an order of mag more than enough
% starti = 2;
numtomake = 10;
plotbar = 1;

% runnum = 4;
dir =sprintf('/home/droumis/MATLAB/TI_behavior/%s/%s/',animal,pairsetname);
mkdir(dir);
cd(dir);
done = '';
%create a rand list of pairs with interp home
for i = 1:numtomake; %use to make multiple without a check
    %    while isempty(done);
    % i = 1;
    % if i = numtomake;
    %     break
    % end
    %
    rng('shuffle'); %reinitialize the random number generator
    
    % randnohome = randi(pairs,lengthlist, 1);
    % randhome = zeros(length(randnohome)*2,1); %create zero array
    % randhome(2:2:end) = randnohome; %interlace rand trials
    
    %create weighted distribution of randnums
    mult = floor(lengthlist/length(pairsweighted));
    ordnohome = [];
    for k = 1:length(pairsweighted);
        ordnohome = [ordnohome; ones(mult, 1)*pairsweighted(k)];
    end
    randnohome = ordnohome(randperm(length(ordnohome)));
    randhome = zeros(length(randnohome)*2,1); %create zero array
    randhome(2:2:end) = randnohome; %interlace rand trials
    
    %OLD
    % ab(i) = length(randhome(randhome(1:100) == 1));
    %     bc(i) = length(randhome(randhome(1:100) == 2));
    %     de(i) = length(randhome(randhome(1:100) == 3));
    %     ef(i) = length(randhome(randhome(1:100) == 4));
    %     CD(i) = length(randhome(randhome(1:100) == 5));
    %     bd(i) = length(randhome(randhome(1:100) == 6));
    %     be(i) = length(randhome(randhome(1:100) == 7));
    %     ce(i) = length(randhome(randhome(1:100) == 8));
    %     af(i) = length(randhome(randhome(1:100) == 9));
    
    ab(i) = length(randhome(randhome(1:160) == 1));
    bc(i) = length(randhome(randhome(1:160) == 2));
    de(i) = length(randhome(randhome(1:160) == 3));
    ef(i) = length(randhome(randhome(1:160) == 4));
    CD(i) = length(randhome(randhome(1:160) == 5));
    be(i) = length(randhome(randhome(1:160) == 6));
    af(i) = length(randhome(randhome(1:160) == 7));
    bd(i) = length(randhome(randhome(1:160) == 8));
    ce(i) = length(randhome(randhome(1:160) == 9));
    
    
    % %insert check to make sure that the pairs are reasonably balanced.. if not, don't save it!
    % if ab < be & ab < CD & bc <be & bc < CD & de < be & de < CD & ef <be & ef<CD & ab<af & bc<af & de < af & ef <af;
    %     done = 'yep';
    % end
    
    fid = fopen(sprintf('%d_%s_%s.txt',i,animal,pairsetname),'w'); %need to get rid of floating points
    fprintf(fid, '%d \n', randhome);
    fclose(fid);
    
end
% fid = fopen(sprintf('%d_%s_%s.txt',i,animal,pairsetname),'w'); %need to get rid of floating points
% fprintf(fid, '%d \n', randhome);
% fclose(fid);

if plotbar == 1;
    figure
    bar([1 2 3 4 5 6 7 8 9], [mean(ab) mean(bc) mean(de) mean(ef) mean(CD) mean(be) mean(af) mean(bd) mean(ce)]);  hold on;
    errorbar2([1 2 3 4 5 6 7 8 9], [mean(ab) mean(bc) mean(de) mean(ef) mean(CD) mean(be) mean(af) mean(bd) mean(ce) ],  [stderr(ab) stderr(bc) stderr(de) stderr(ef) stderr(CD) stderr(be) stderr(af) stderr(bd) stderr(ce)] , 0.3, 'k')
end

