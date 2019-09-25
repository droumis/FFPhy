% RANDCIRC: Generate uniform-random points within a circle of given radius and
%           center.
%
%     Usage: crds = randcirc({npts},{r},{center},{doplot})
%
%           npts =   optional number of points to be generated [default=1].
%           r =      optional radius of circle [default=1].
%           center = optional 2-element vector specifying center coordinate of 
%                      circle [default = origin].
%           doplot = optional boolean flag indicating that plot of points and 
%                      circle is to be produced [default = 0].
%           ------------------------------------------------------------------
%           crds =   [npts x 2] matrix of cartesian point coordinates
%

% RE Strauss, 9/22/95
%   3/19/00 - miscellaneous improvements; added optional plot.
%   3/21/00 - changed call to circcrds().

function crds = randcirc(npts,r,center,doplot)
  if (nargin < 1) npts = []; end;
  if (nargin < 2) r = []; end;
  if (nargin < 3) center = []; end;
  if (nargin < 4) doplot = []; end;

  if (isempty(npts))
    npts = 1;
  end;
  if (isempty(r))
    r = 1;
  end;
  if (isempty(center))
    center = [0 0];
  end;
  if (isempty(doplot))
    doplot = 0;
  end;

  crds = zeros(npts,2);
  for i = 1:npts
    while (crds(i)==[0 0])
      p = 2*r*rand(1,2) - [r r];     % Generate points within square
      d = eucl(p,[0 0]);             %   calculate dists from center
      if (d <= r)                    %   and filter
        crds(i,:) = p;
      end;
    end;
  end;

  crds = crds + ones(npts,1)*center;

  if (doplot)
    ccrds = circcrds(r,center);

    plot(crds(:,1),crds(:,2),'ok');
    hold on;
      plot(ccrds(:,1),ccrds(:,2),'k');
    hold off;
    sqplot(ccrds);
  end;

  return;
