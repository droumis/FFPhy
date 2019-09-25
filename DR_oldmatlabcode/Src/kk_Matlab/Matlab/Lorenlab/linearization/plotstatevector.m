function plotstatevector(statevector, lindist, time)
%plotstatevector(statevector, lindist, time)
%
%plots linear distance as a function of time. The plot is color coded for
%which trajectory the animal is on.
%

figure  
% ind = find(statevector == -1);
% plot(time(ind),lindist(ind),'k.');
% hold on;
% ind = find(statevector == 1);
% plot(time(ind),lindist(ind),'b.');
% ind = find(statevector == 2);
% plot(time(ind),lindist(ind),'c.');
% ind = find(statevector == 3);
% plot(time(ind),lindist(ind),'r.');
% ind = find(statevector == 4);
% plot(time(ind),lindist(ind),'m.');

for i = 1:(max(statevector))
   ind = find(statevector == i);
   %get all the linear locations when the animal was on the ith
   %trajectory
   switch mod(i,6)
    case 1
      col = 'b';
    case 2
      col = 'c';
    case 3
      col = 'r';
    case 4
      col = 'm';
    case 5
      col = 'g';
    case 0
      col = 'y';
   end
   dot = '.';
   if ((i > 6) & (i <= 12))
        dot = '+';
   elseif (i>=13)
        dot = '*';
   end
   plot((ind),lindist(ind),[col,dot]);
   hold on
   
end
ind = find(statevector == -1);
plot((ind),lindist(ind),'k.','MarkerSize',4);