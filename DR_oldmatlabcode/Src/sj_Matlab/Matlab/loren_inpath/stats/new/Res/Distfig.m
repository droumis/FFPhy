% DISTFIG:  Plots a figure of landmarks and distances among landmarks, given the 
%             coordinates and distance specificiations.
%
%     Usage: distfig(specs,crds,{landmarks},{specimens},{individ})
%
%           specs =     [m x 2], [m x 3] or [m x 4] set of interpoint distance 
%                         specifications [see function distance()].
%           crds =      [n x p] set of point coordinates.
%           landmarks = optional vector specifying the indices of landmark 
%                         points to be plotted [default = all points in specs].
%           specimens = optional vector [length n] specifying k individuals, each 
%                         represented by the same number (np) of point coordinates 
%                         (so that np x k = n).  [Default = all coordinates for a  
%                         single individual].
%           individ =   optional vector of individuals (from the list of 
%                         specimens) to be ploted [default = first individual].
%

% RE Strauss, 4/20/00

function distfig(specs,crds,landmarks,specimens,individ)
  if (nargin < 3) landmarks = []; end;
  if (nargin < 4) specimens = []; end;
  if (nargin < 5) individ = []; end;

  [nspec,c] = size(specs);
  switch (c)
    case 2
    case 3,
      specs = specs(:,2:3);
    case 4,
      specs = specs(:,2:3);
    otherwise
      error('  DISTFIG: invalid specification matrix');
  end;
  if (isempty(landmarks))
    landmarks = uniquef(specs(:));
  else
    if (min(landmarks)<1 | max(landmarks)>max(specs(:)))
      error('  DISTFIG: invalid landmarks list');
    end;
  end;

  [ncrds,c] = size(crds);
  if (c~=2)
    error('  DISTFIG: invalid coordinate matrix');
  end;

  if (isempty(specimens))
    specimens = ones(ncrds,1);
  end;
  if (isempty(individ))
    individ = 1;
  end;
  nind = length(individ);

  ind_title = 0;
  if (nind>1)
    ind_title = 1;
  end;

  for ind = 1:nind
    [i,nf] = getobs(specimens,individ(ind));
    if (~nf)
      pts = crds(i,:);
      figure;
      plot(pts(landmarks,1),pts(landmarks,2),'ko');
      axis('equal');
      axis('off');
      hold on;
      for is = 1:nspec
        plot(pts(specs(is,:),1),pts(specs(is,:),2),'k');
      end;
      hold off;
      if (ind_title)
        puttitle(sprintf('Specimen %d',individ(ind)));
      end;
    else
      disp(sprintf('  DISTFIG warning: individual %d not found',individ(ind)));
    end;
  end;

  return;
