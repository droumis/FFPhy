% UPGMARND: Randomizes a UPGMA to test the existence of a binomial 
%           distribution of null support values.
%

N = 4;                                        % Number of taxa (rooted)
P = [3:10,12,15:5:25,30:10:100];              % Numbers of characters
boot_iter = [50,100,200,500,1000];            % Numbers of bootstrap iterations
repl = 30;                                    % Number of replicates

nchar = length(P);
nboot = length(boot_iter);

freq = zeros(nchar*nboot,repl);   % Allocate bootstrap-support frequencies
spec = zeros(nchar*nboot,2);
dist = zeros(N,N);                % Allocate distance matrix
row = 0;

for c = 1:nchar
  ch = P(c);

  for b = 1:nboot
    bi = boot_iter(b);

    row = row+1;
    spec(row,:) = [ch bi];          % Stash specifications for current run
    chars_iter = [ch bi]            % Screen display

    for r = 1:repl                  % Create a distribution of bootstrap
      s = 0;                        %   support frequencies
%      X = rand(N,ch);                 % Create random data matrix - continuous
      X = round(rand(N,ch));                 % Create random data matrix - discrete

      for i = 1:bi
        Xb = bootsamp(X')';         % Bootstrap characters from the data matrix

        t = randperm(N);              % Shuffle taxa for UPGMA
        t12 = find(t==1 | t==2);      % Keep track of positions of taxa 1,2
        Xb = Xb(t,:);

        for n1 = 1:(N-1)            % Pairwise distances
          for n2 = 2:N
%            d = Xb(n1,:)-Xb(n2,:);	% Euclidean
%            d = sum(d.*d);  
	    d = sum(abs(Xb(n1,:)-Xb(n2,:)));
            dist(n1,n2) = d;
            dist(n2,n1) = d;
          end;
        end;

        topology = upgma(dist,[],1);  % Cluster and find topology

        for n = 1:(N-2)               % Search tree topology
          if (topology(n,1:2)==t12)
%          if (topology(n,1:2)==[1,2])   % If contains cluster 1,2
            s = s+1;                    %   add to support
          end;
        end;
      end;

      freq(row,r) = s/bi;
    end;
  end;

  save waleed;
end;


