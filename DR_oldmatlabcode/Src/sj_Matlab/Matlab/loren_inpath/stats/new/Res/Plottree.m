% PLOTTREE: Plots the topology of the tree (t1,(t2,(t3,t4))) for terminal 
%           taxa t1-t4, given a vector 'brlen' containing the six branch lengths: 
%             brlen(1) = (a1234)-t1
%             brlen(2) = (a1234)-(a234)
%             brlen(3) = (a234)-t2
%             brlen(4) = (a234)-(a34)
%             brlen(5) = (a34)-t3
%             brlen(6) = (a34)-t4
%
%     Syntax: plottree(brlen)
%

% RE Strauss, 3/5/97
%   5/2/03 - changed call to triangpt().

function plottree(brlen)
  t = [1:4];
  ax = [1:3];
  ay = [1:3];

  P = triangpt([t(3) 0; t(4) 0],brlen(5:6));
  ax(1) = P(2,1);
  ay(1) = P(2,2);
  
  P = triangpt([t(2) 0; ax(1) ay(1)],brlen(3:4));
  ax(2) = P(2,1);
  ay(2) = P(2,2);

  P = triangpt([t(1) 0; ax(2) ay(2)],brlen(1:2));
  ax(3) = P(2,1);
  ay(3) = P(2,2);


  plot([t(3) ax(1) t(4)],[0 ay(1) 0],'b');
  hold on;
  plot([t(2) ax(2) ax(1)],[0 ay(2) ay(1)],'b');
  plot([t(1) ax(3) ax(2)],[0 ay(3) ay(2)],'b');

  axis([1 t(4) ay(3) 0]);              % Axes
  axis square;
  axis off;

  hold off;
  return;
