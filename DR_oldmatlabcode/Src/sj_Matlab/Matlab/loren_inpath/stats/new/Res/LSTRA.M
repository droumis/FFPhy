% LSTRA:  Orthogonal least-squares theta-rho analysis for the Procrustes
%         mapping of one form onto another in any number of dimensions.
%
%     Usage: [map1,map2,mapa,sse] = lstra(form1,form2,{form2a},{noplot})
%
%         form1,form2 = [p x k] matrices of point configurations for p points
%                         and k dimensions; form2 is mapped onto form1.
%         form2a =      optional [q x k] set of auxiliary points to be mapped 
%                         along with form2.
%         noplot =      optional boolean flag indicating that no plots are to be 
%                         produced [default = 0].
%         ----------------------------------------------------------------------
%         map1,map2 =   [p x k] matrices of mapped point configurations,
%                         scaled and zero-centered.
%         mapa =        [q x k] matrix of mapped auxiliary points, to same scale 
%                         as map2.
%         sse =         sum of squared distances between corresponding points.
%

% Rohlf and Slice. 1990. Extensions of the Procrustes method for the optimal
%   superimposition of landmarks.  Syst. Zool. 39:40-59.

% RE Strauss, 9/1/95
%   11/24/98 - allow auxiliary set of points to be carried along with form2.
%    4/30/00 - changed default for plot.
%    6/11/03 - convert input vectors to column vectors;
%              plot only 2D points.

function [map1,map2,map2a,sse] = lstra(form1,form2,form2a,noplot)
  if (nargin < 3) form2a = []; end;
  if (nargin < 4) noplot = []; end;
  
  if (isvector(form1))
    form1 = form1(:);
  end;
  if (isvector(form2))
    form2 = form2(:);
  end;

  [p1,k] = size(form1);
  [p2,k2] = size(form2);

  if (isempty(form2a))
    auxil = 0;
    q = p2;
    k2a = k2;
    map2a = [];
  else
    [q,k2a] = size(form2a);
    auxil = 1;
  end;

  if (isempty(noplot))
    noplot = 0;
  end;

  if (any(abs([p1-p2,k-k2,k2-k2a])))    % Check conformability of forms
    error('  LSTRA: Forms must be identical in dimension and numbers of points');
  end;

  % Translate the centroids of the two forms to a common origin

  centr1 = mean(form1);
  centr2 = mean(form2);

  map1 = form1 - ones(p1,1)*centr1;
  map2 = form2 - ones(p2,1)*centr2;
  if (auxil)
    map2a = form2a - ones(q,1)*centr2;
  end;

  % Scale each form by the square root of the sum of the squared distances
  % from the centroid

  ip = eye(p1) - ones(p1)/p1;
  s1 = sqrt(trace(ip*map1*map1'*ip));
  s2 = sqrt(trace(ip*map2*map2'*ip));
  map1 = map1./s1;
  map2 = map2./s2;
  if (auxil)
    map2a = map2a./s2;
  end;

  % Rotate the second form to achieve the best approximation to the first.
  % If 2D, iteratively optimize the sse (otherwise too many parameters).

  [U,sigma,V] = svd(map1'*map2);
  S = sign(sigma);
  H = V*S*U';

  map2 = map2*H;
  if (auxil)
    map2a = map2a*H;
  end;

  % Compute sse

  sse = trace((map1-map2)*(map1-map2)');

  if (~noplot & k==2)                      % Plot 2D forms
    plot(map1(:,1),map1(:,2),'ro');          % Reference form
    hold on;
    plot(map2(:,1),map2(:,2),'bo');          % Mapped form
    for i = 1:size(map2,1)
      plot([map1(i,1) map2(i,1)],[map1(i,2) map2(i,2)],'k');
    end;

    xcrd = [map1(:,1) map2(:,1)];
    [maxxcrd,i] = max(xcrd');

    textcrd = map1;
    j = find(i==2);
    textcrd(j,:) = map2(j,:);
    textcrd(:,1) = textcrd(:,1) + 0.02*range([map1(:,1);map2(:,1)]);

    for i = 1:size(map1,1)
      text(textcrd(i,1),textcrd(i,2),tostr(i));
    end;

    if (auxil)
      plot(map2a(:,1),map2a(:,2),'ko');      % Auxiliary landmarks
    end;
    axis('square');
    axis(sqplot([map1;map2;map2a]));
    hold off;
  end;

  return;
