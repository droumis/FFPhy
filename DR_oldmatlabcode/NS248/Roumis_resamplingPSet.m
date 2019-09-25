
%Demetris Roumis
%Pset Resampling
%11/3/2014

%% Problem 1a

N=10;
B = 1000;
[X,Y] = GetSample(N);
D = mean(X) - mean(Y);
d = zeros(B,1);
for b = 1:B;
    Xboot = X(ceil(N*rand(N,1)));
    Yboot = Y(ceil(N*rand(N,1)));
    d(b) = mean(Xboot) - mean(Yboot);
end
sQ = std(d);
CLQ = prctile(d,[2.5 97.5]);
sigt = D < CLQ(1) || D > CLQ(2);
%CLQ overlaps zero, probably not sig


%% problem 1b
dpair = zeros(B,1);
for b = 1:B;
    ind = ceil(N*rand(N,1));
    Xpairboot = X(ind);
    Ypairboot = Y(ind);
    dpair(b) = mean(Xpairboot) - mean(Ypairboot);
end
sQpair = std(dpair);
CLQpair = prctile(dpair,[2.5 97.5]);
sigtpair = D < CLQpair(1) || D > CLQpair(2);
%CLQpair does not overlap zero, probably sig

%% Problem 2a

Dperm = zeros(B,1);
Z = [X;Y];
for r = 1:B;
    Zperm = Z(randperm(size(Z)));
    Dperm(r) = mean(Zperm(1:N,:)) - mean(Zperm((N+1):(2*N),:));
end
pperm = mean(abs(Dperm)> abs(D)); %pperm = .32

%% Problem 2b
Zpair = [X Y];
Dpairperm = zeros(B,1);
for r = 1:B;
    rmat = round(rand(N,1));
    for j = 1:10;
        if rmat(j)
            Zpairperm(j,:) = fliplr(Zpair(j,:));
        else
            Zpairperm(j,:) = Zpair(j,:);
        end
    end
    Dpairperm(r) = mean(Zpairperm(:,1)) - mean(Zpairperm(:,2));
end
ppairperm = mean(abs(Dpairperm)> abs(D)); % = .006

%% Problem 3a
w = warning ('off','all');
Nsamps = [10 20];
Diffs = [0.12 0.25 0.5 1 2];
dSets = 200;
B = 500;
pperm = zeros(dSets,1);
power = zeros(i,j);

for i = 1:length(Nsamps);
    for j = 1:length(Diffs);
        for l = 1:dSets;
            [X, Y] = GetSample(Nsamps(i),Diffs(j));
            D = mean(X) - mean(Y);
            Z = [X;Y];
            Dperm = zeros(B,1);
            for k = 1:B;
                Zperm = Z(randperm(size(Z)));
                Dperm(k) = mean(Zperm(1:Nsamps(i),:)) - mean(Zperm((Nsamps(i)+1):(2*Nsamps(i)),:));
            end
            pperm(l) = (mean(abs(Dperm)> abs(D)))<0.05;
        end
        power(i,j) = mean(pperm);
    end
end
figure;
plot(Diffs,power(1,:),'-*b'); hold on;
plot(Diffs,power(2,:),'-or')
title('upaired')
xlabel('Effect Size')
ylabel('power')
legend('N=10', 'N=20')



%% Problem 3b
w = warning ('off','all');
Nsamps = [10 20];
Diffs = [0.12 0.25 0.5 1 2];
dSets = 200;
B = 500;
pperm = zeros(dSets,1);
power = zeros(i,j);

for i = 1:length(Nsamps);
    for j = 1:length(Diffs);
        for l = 1:dSets;
            [X, Y] = GetSample(Nsamps(i),Diffs(j));
            D = mean(X) - mean(Y);
            Zpair = [X Y];
            Dpairperm = zeros(B,1);
            for k = 1:B;
                rmat = round(rand(Nsamps(i),1));
                for t = 1:Nsamps(i);
                    if rmat(t)
                        Zpairperm(t,:) = fliplr(Zpair(t,:));
                    else
                        Zpairperm(t,:) = Zpair(t,:);
                    end
                end
                Dpairperm(k) = mean(Zpairperm(:,1)) - mean(Zpairperm(:,2));
            end
            pperm(l) = (mean(abs(Dpairperm)> abs(D)))<0.05;
        end
        power(i,j) = mean(pperm);
    end
end

figure;
plot(Diffs,power(1,:),'-*b'); hold on;
plot(Diffs,power(2,:),'-or')
title('paired')
xlabel('Effect Size')
ylabel('power')
legend('N=10', 'N=20')
































