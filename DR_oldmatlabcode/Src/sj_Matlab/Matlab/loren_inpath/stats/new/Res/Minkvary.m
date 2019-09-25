% MINKVARY: Maps the values of a single continuous character onto a tree for 
%           varying values of k, the Minkowsky exponent, producing profiles of 
%           fitted node values, fitted branch lengths, and several tree 
%           properties: total tree length, variance of branch lengths, homoplasy 
%           (sum of absolute deviations between observed pairwise differences 
%           and patristic distances), and sum of absolute deviations between 
%           fitted and specified branch lengths.
%
%     Usage: [k,nodes,fbrlen,treelen,varbrlen,homoplasy,devbrlen] = ...
%                      minkvary(anc,X,{rooted},{brlen},{ktree},{kbnds},{nk})
%
%           anc =       ancestor function describing tree topology; the 
%                         root (arbitrary or real) must have ancestor 0.
%           X =         [n x p] data matrix for n taxa and p characters.
%           rooted =    optional boolean flag indicating whether tree root is (=1)
%                         or is not (=0) to be optimized; if not, the tree is
%                         treated as an arbitrarily rooted network [default = 1].
%           brlen =     optional extrinsic branch lengths for weighted 
%                         least-squares mapping, in same sequence as anc 
%                         [default = 1 for all branches].
%           ktree =     optional vector of values of Minkowsky k for which plots 
%                         of trees are to be produced [default = null].
%           kbnds =     optional 2-element vector of [kmin,kmax], the lower and 
%                         upper bounds of k to be mapped [default = [0.5,3]].
%           nk =        optional length of vector of k values [default = 150].
%           --------------------------------------------------------------------
%           k =         vector (length nk) of k values, from kmin to kmax.
%           nodes =     matrix of fitted node values as a function of k.
%           fbrlen =    matrix of fitted branch lengths as a function of k.
%           treelen =   total tree length as a function of k.
%           varbrlen =  variance in fitted branch lengths as a function of k.
%           homoplasy = total homoplasy as a function of k.
%           devbrlen =  total deviation of fitted branch lengths from extrinsic 
%                         branch lengths as a function of k.
%

% RE Strauss, 5/23/01

function [k,nodes,fbrlen,treelen,varbrlen,homoplasy,devbrlen] = ...
                                  minkvary(anc,X,rooted,brlen,ktree,kbnds,nk)

  if (nargin < 3) rooted = []; end;
  if (nargin < 4) brlen = []; end;
  if (nargin < 5) ktree = []; end;
  if (nargin < 6) kbnds = []; end;
  if (nargin < 7) nk = []; end;

  if (isempty(rooted))
    rooted = 1;
  end;
  if (isempty(nk))
    nk = 150;
  end;

  def_brlen = 0;
  if (isempty(brlen))
    brlen = ones(size(anc));
    def_brlen = 1;
  end;
  if (isempty(kbnds))
    kbnds = [0.5,3];
  elseif (length(kbnds)~=2)
    error('  MINKVARY: invalid kbnds vector.');
  end;

  if (length(brlen)~=length(anc))
    error('  MINKVARY: branch-length and ancestor-fn vectors incompatible.');
  end;
  if (any(brlen < 0))
    disp('  MINKVARY warning: negative branch lengths set to zero');
    i = find(brlen < 0);
    brlen(i) = zeros(length(i),1);
  end;

  if (~rooted)                          % If tree not rooted, remove root
    root = find(anc == 0);              % Find arbitrary root node
    desc = find(anc == root);           % Find its descendants
    d1 = desc(1);
    d2 = desc(2);
    anc(d1) = d2;                       % Form trichotomy
    if (~def_brlen)
      brlen(d1) = brlen(d1)+brlen(d2);  % Sum the branch lengths
    end;
    indx = find(anc);
    anc = anc(indx);                    % Remove root node
    brlen = brlen(indx);
  end;

  k = linspace(kbnds(1),kbnds(2),nk)';      % Minkowski k
  nnodes = length(X)-1;                     % Number of nodes
  nbranches = 2*length(X)-2;

  nodevals = zeros(nk,nnodes);            
  for i = 1:nk                            % Fit varying values of k
    nodevals(i,:) = contmap(anc,X,rooted,brlen,k(i),1)';
  end;

  XX = [ones(nk,1)*X nodevals];           % Extend data to match node sets
  fbrlen = zeros(nk,nbranches);           % Calculate fitted branch lengths
  for i = 1:nbranches
    fbrlen(:,i) = abs(XX(:,i)-XX(:,anc(i)));
  end;

  meantrue = mean(brlen);
  brlen = ones(nk,1)*brlen;                 % Duplicate true brlens for 'nk' rows
  meanfitted = mean(fbrlen')';
  sfbrlen = fbrlen;
  for i = 1:nk                              % Scale fitted to match true brlens
    sfbrlen(i,:) = fbrlen(i,:).*(meantrue/meanfitted(i));
  end;

  % Criteria

  homoplasy = zeros(nk,1);
  obsdists = eucl(X');                      % Observed Euclidean distances among taxa
  for i = 1:nk                              % Calculate patristic distances
    [step,p] = patrist(anc,fbrlen(i,:));      % Patristic distances
    devs = abs(p-obsdists);                   % Deviations of patristic from observed
    homoplasy(i) = sum(trilow(devs));         % Sum of patristic deviations
  end;

  treelen = rowsum(fbrlen);                 % Total tree length
  varbrlen = var(fbrlen')';                 % Variance among branch lengths

  devbrlen = rowsum(abs(sfbrlen-brlen));    % Deviations of fitted branch lengths from true


  % Plots

  ta = sprintf('%d,',anc);
  ta(end) = [];
  tb = sprintf('%d,',brlen);
  tb(end) = [];
  title = ['anc=[',ta,'], w=[',tb,']'];

  figure;                                   % Figure 1
  kk = k*ones(1,nnodes);                      % Bounds for node-value plots
  v = putbnds(kk(:),nodevals(:),0.05,1);
  for i = 1:nnodes                            % Plot fitted node values
    subplot(nnodes,1,i);
    plot(k,nodevals(:,i),'k');
    axis(v);
    puttick([],[],8);
    putylab(sprintf('Node %d',i),[],10);
    if (i==nnodes)
      putxlab('Minkowski k');
    end;
    if (i==1)
      puttitle(title,12);
    end;
  end;

  figure;                                   % Figure 2
  kk = k*ones(1,nbranches);                   % Bounds for fitted branch-length plots
  v = putbnds(kk(:),fbrlen(:),0.05,1);
  for i = 1:nbranches                         % Plot fitted branch lengths
    subplot(nbranches,1,i);
    plot(k,fbrlen(:,i),'k');
    axis(v);
    puttick([],[],8);
    putylab(sprintf('Br %d',i),[],10);
    if (i==nbranches)
      putxlab('Minkowski k');
    end;
    if (i==1)
      puttitle(title,12);
    end;
  end;

  figure;                                   % Figure 3
  subplot(4,1,1);
  plot(k,treelen,'k');
  putbnd(k,treelen);
  putylab('Tree len',[],12);
  puttick([],[],8);
  puttitle(title,12);

  subplot(4,1,2);
  plot(k,varbrlen,'k');
  putbnd(k,varbrlen);
  putylab('Brlen var',[],12);
  puttick([],[],8);

  subplot(4,1,3);
  plot(k,homoplasy,'k');
  if (range(homoplasy)<0.1)
    v(3:4) = [homoplasy(1)-0.5 homoplasy(1)+0.5];
    axis(v);
  else
    putbnd(k,homoplasy);
  end;
  putylab('Homoplasy',[],12);
  puttick([],[],8);

  subplot(4,1,4);
  plot(k,devbrlen,'k');
  if (range(devbrlen)<0.1)
    v(3:4) = [devbrlen(1)-0.5 devbrlen(1)+0.5];
    axis(v);
  else
    putbnd(k,devbrlen);
  end;
  putylab('Brlen dev',[],12);
  puttick([],[],8);

  putxlab('Minkowski k');

  return;



