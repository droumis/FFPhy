% ALLCLADE: Plots probabilities of all possible clade sizes for a rooted tree 
%           with T terminal taxa.
%
%     Usage: allclade(taxa)
%
%           taxa = vector of taxon sizes (T)
%

% RE Strauss, 8/19/98

function allclade(taxa)
  close all;

%  taxa = [5:15];
  lentaxa = length(taxa);
  cum = [];

  hold on;
  for t = 1:lentaxa
    T = taxa(t);
    p = zeros(T-2,1);
    for i = 2:(T-1)
      p(i-1) = cladeprb(i,T);
    end;
    plot([2:(T-1)]',p,'k');

    if (t==1)
      text(T-0.9,p(T-2),['T=',int2str(T)]);
    else
      text(T-0.9,p(T-2),int2str(T));
    end;

    cum = [cum; [2:(T-1)]',p];
  end;

  putbnd(cum(:,1),cum(:,2));
  putxbnd(1,T);
  putxlab('Clade Size (t)');
  putylab('Probability');
  hold off;

  figure;
  hold on;
  for t = 1:lentaxa
    T = taxa(t);
    p = zeros(T-2,1);
    for i = 2:(T-1)
      p(i-1) = cladeprb(i,T)*100;
    end;
    plot([2:(T-1)]',p,'k');

    if (t==1)
      text(T-0.9,p(T-2),['T=',int2str(T)]);
    else
      text(T-0.9,p(T-2),int2str(T));
    end;

    cum = [cum; [2:(T-1)]',p];
  end;

  putbnd(cum(:,1),cum(:,2));
  putxbnd(1,T);
  putxlab('Clade Size (t)');
  putylab('Expected BSF');
  hold off;

  return;
