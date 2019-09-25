% CORRCIRC: Graphically displays a correlation value, plus/minus confidence
%           limits, in circular form.
%
%   Usage: corrcirc(low,corr,high)
%

% RE Strauss, 6/4/94
%   8/19/99 - changed plot colors and characteristics for Matlab v5.

function corrcirc(low,corr,high)
  r = .13;                              % Radius of circle in unit box
  w = .008;                             % Width of corr slice
  labels = 0;

  t = pi-(-1:.01:1)'*pi;                % Draw circle
  x = [0; r*cos(t); 0];
  y = [0; r*sin(t); 0];
  plot(x,y,'k');

  axis([-1,1,-1,1]);                    % Set plot characteristics
  axis('square');
  axis('equal');
  axis('off');
  if labels
    text(r+0.04,.07,'+1');
    text(r+0.04,-.07,'-1');
    text(-r-0.09,0,'0');
  end;
  hold on;

  plot([-r,0],[0,0],'k:');              % Draw dotted line at corr=0

  t = pi-(low:w:(corr-w))'*pi;          % Draw low confidence slice
  x = [0; r*cos(t); 0];
  y = [0; r*sin(t); 0];
  fill(x,y,'b');

  t = pi-((corr+w):w:high)'*pi;         % Draw high confidence slice
  x = [0; r*cos(t); 0];
  y = [0; r*sin(t); 0];
  fill(x,y,'b');

  t = pi-((corr-w):w/10:(corr+w))'*pi;  % Draw corr value
  x = [0; 1.1*r*cos(t); 0];
  y = [0; 1.1*r*sin(t); 0];
  plot(x,y,'k');

  hold off;
  return;

