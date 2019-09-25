% UPGMABT: Bootstraps a UPGMA to test the existence of a binomial 
%          distribution of sampled support values.
%

N = 4;                                        % Number of taxa (rooted)
tv = [1 2 4 5.5];                             % Taxon mean character values

P = [3:10,12,15:5:25,30:10:100];              % Numbers of characters
boot_iter = [50,100,200,500,1000];            % Numbers of bootstrap iterations
repl = 30;                                    % Number of replicates

nchar = length(P);
nboot = length(boot_iter);

maxchar = P(nchar);
XX = zeros(N,maxchar);
for i = 1:N                       % Construct data matrix
  XX(i,:) = tv(i)+randn(1,maxchar);
end;
XX

freq = zeros(nchar*nboot,repl);   % Allocate bootstrap-support frequencies
spec = zeros(nchar*nboot,2);
dist = zeros(N,N);                % Allocate distance matrix
row = 0;

for n1 = 1:(N-1)            % Pairwise Euclidean distances
  for n2 = 2:N
    d = XX(n1,:)-XX(n2,:);
    d = sum(d.*d);  
    dist(n1,n2) = d;
    dist(n2,n1) = d;
  end;
end;
dist
topology = upgma(dist,[],1)   % Cluster and find topology
fc = topology(1,1:2);         % First cluster

for c = 1:nchar
  ch = P(c);
  seq = randperm(maxchar);          % Random sample of chars from data matrix
  X = XX(:,seq(1:ch));

  for b = 1:nboot
    bi = boot_iter(b);

    row = row+1;
    spec(row,:) = [ch bi];          % Stash specifications for current run
    chars_iter = [ch bi]            % Screen display

    for r = 1:repl                  % Create a distribution of bootstrap
      s = 0;                        %   support frequencies

      for i = 1:bi
        Xb = bootsamp(X')';         % Bootstrap characters from the data matrix

        for n1 = 1:(N-1)            % Pairwise Euclidean distances
          for n2 = 2:N
            d = Xb(n1,:)-Xb(n2,:);
            d = sum(d.*d);  
            dist(n1,n2) = d;
            dist(n2,n1) = d;
          end;
        end;

        topology = upgma(dist,[],1);  % Cluster and find topology

        for n = 1:(N-2)               % Search tree topology
          if (topology(n,1:2)==fc)   % If contains cluster 1,2
            s = s+1;                    %   add to support
          end;
        end;
      end;

      freq(row,r) = s/bi;
    end;
  end;

  save waleed;
end;


