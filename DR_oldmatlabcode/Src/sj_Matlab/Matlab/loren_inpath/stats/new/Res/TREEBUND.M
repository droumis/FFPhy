% TREEBUND: Plots a bundle of mapped trees on a 2-dimensional space.
%
%       Usage: treebund(X,grps,anc,ntrees)
%

% RE Strauss, 7/31/95
%   8/20/99 - changed plot colors for Matlab v5.

function treebund(X,grps,anc,ntrees)
  grp_name = uniquef(grps);              % Group identifiers
  ngrps = length(grp_name);             % Number of groups

  bm = zeros(ngrps,2);                  % Bootstrapped means

  clf;
  hold on;
  for g = 1:ngrps                         % Plot hulls and centroids
    j = (grps==g);
    c = mean(X(j,:));
    h = hull(X(j,:));
    plot(h(:,1),h(:,2),'k:');
    plot(c(1),c(2),'k*');
  end;

  for nt = 1:ntrees
    for g = 1:ngrps                     % Bootstrap the group means
      j = (grps==g);
      y = bootsamp(X(j,:));
      bm(g,:) = mean(y);
    end;
    node = contmap(anc,bm);

    plot([bm(1,1),node(3,1)],[bm(1,2),node(3,2)],'r-');
    plot([bm(2,1),node(3,1)],[bm(2,2),node(3,2)],'r-');
    plot([bm(3,1),node(2,1)],[bm(3,2),node(2,2)],'r-');
    plot([bm(4,1),node(1,1)],[bm(4,2),node(1,2)],'r-');
    plot([bm(5,1),node(1,1)],[bm(5,2),node(1,2)],'r-');
    plot([node(1,1),node(2,1)],[node(1,2),node(2,2)],'r-');
    plot([node(2,1),node(4,1)],[node(2,2),node(4,2)],'r-');
    plot([node(3,1),node(4,1)],[node(3,2),node(4,2)],'r-');

  end;
  
  return;
