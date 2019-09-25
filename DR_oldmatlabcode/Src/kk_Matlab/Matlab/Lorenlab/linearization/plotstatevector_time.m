function plotstatevector_time(statevector, lindist, time)
%plotstatevector(statevector, lindist, time)
%
%plots linear distance as a function of time. The plot is color coded for
%which trajectory the animal is on.
%

figure  
 ind = find(statevector == -1);
 plot(time(ind),lindist(ind),'k.');
 hold on;
 ind = find(statevector == 1);
 plot(time(ind),lindist(ind),'b.');
 ind = find(statevector == 2);
 plot(time(ind),lindist(ind),'c.');
 ind = find(statevector == 3);
 plot(time(ind),lindist(ind),'r.');
 ind = find(statevector == 4);
 plot(time(ind),lindist(ind),'m.');