% Plotbino: Plots binomial distributions.  'Arrows' is an optional
%           vector of X values indicating histogram bars to be marked by 
%           downward-pointing arrows.
%
%     Syntax:  plotbino(N,p,arrows)
%
%           N,p =     parameters of binomial distribution.
%           arrows =  optional vector of X values indicating histogram bars to
%                       be marked by downward-pointing arrows.  Use negative 
%                       value for arrow marking mean of distribution.
%

% RE Strauss, 3/2/99
%   8/19/99 - change due to arrowdwn()
%   6/17/01 - change name from binoplot()

function plotbino(N,p,arrows)
  if (nargin<3) arrows = []; end;

  N = N(:);
  p = p(:);
  lenN = length(N);
  lenp = length(p);

  if (lenN ~= lenp)
    error('BINOPLOT: N and p vectors not compatible');
  end;

  for i = 1:lenN
    y = binopdf(0:N(i),N(i),p(i));
    x = [[0:N(i)]'; makegrps(0:N(i),round(y*1000))];
    if (isempty(arrows))
      histgram(x,[],[],N(i)+1,[],1,'rel');
    elseif (arrows<0)
      histgram(x,[],[],N(i)+1,[],0,'rel');
    else
      histgram(x,[],[],N(i)+1,[],1,'rel');
      ax = axis;
      maxyy = ax(4);

      hold on;
      alength = 0.12;                         % Proportional height of arrow shaft
      binheight = binopdf(arrows,N(i),p(i));
      for j = 1:length(arrows)
        position = binheight(j) + 0.25*(alength*maxyy);   % y-coord of tip
        maxyy = max([maxyy, position+(alength*maxyy)]);   % New maximum
        arrowdwn(arrows(j),position,alength*maxyy);       % Draw arrow
      end;
      hold off;
      
      ax(4) = maxyy + 0.05*maxyy;
      axis(ax);
    end;
  end;

  return;
