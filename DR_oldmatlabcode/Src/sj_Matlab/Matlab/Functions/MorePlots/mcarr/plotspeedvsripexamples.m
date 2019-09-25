% Plot examples of CA1 eeg during ripples at slow, medium, and fast speeds

%Load raw eeg file
load '/data13/mcarr/Eig/EEG/Eigeeg01-2-25.mat'
e = eeg{1}{2}{25};
clear eeg
load '/data13/mcarr/Eig/EEG/Eigripple01-2-25.mat'
r = double(ripple{1}{2}{25}.data(:,1));
clear ripple

%Load ripple events and position info
load '/data13/mcarr/Eig/Eigripples01.mat'
rip = ripples{1}{2}{25};
clear ripples
load '/data13/mcarr/Eig/Eigpos01.mat'
p = pos{1}{2};
clear pos

% Look up speed for all ripple events
speed = p.data(rip.posind(~isnan(rip.posind)),5);

% Pick slow, medium, and fast ripples
slow = find(speed < 1/2);
slow = find(rip.maxthresh == max(rip.maxthresh(slow)));
medium = find (speed < 10 & speed > 5);
medium = find(rip.maxthresh == max(rip.maxthresh(medium)));
fast = find(speed > 15,3);
fast = find(rip.maxthresh == max(rip.maxthresh(fast)));

% Plot three events and spatial trajectories that go with them
figure(1)
hold on
figure(2)
hold on
offset = 0;
for a = [slow medium fast]
%     figure(1)
%     event = e.data(rip.midind(a)-100:rip.midind(a)+100);
%     revent = r(rip.midind(a)-100:rip.midind(a)+100);
%     time = rip.midind(a)-100:rip.midind(a)+100;
%     time = ((time - rip.midind(a))./rip.samprate)*1000;
%     plot(time,event + offset,'k',time + 150,revent+offset,'r')
    figure(2)
    plot(-offset/7 + p.data(:,2),p.data(:,3),'k')
    plot(-offset/7 + p.data(rip.posind(a)-150:rip.posind(a),2),p.data(rip.posind(a)-150:rip.posind(a),3),'r')
    plot(-offset/7 + p.data(rip.posind(a):rip.posind(a)+150,2),p.data(rip.posind(a):rip.posind(a)+150,3),'g')
    plot(-offset/7 + p.data(rip.posind(a),2),p.data(rip.posind(a),3),'c*')
    offset = offset - 750;
end

minima = min(event+offset+750);
scalex = [(time(100)+150-time(1)) (time(100)+150-time(76))];
scaley = [(minima-100) (minima - 300)];
box off
axis off
line(scalex,[scaley(2) scaley(2)],'Color','k')
line([scalex(1) scalex(1)], scaley,'Color','k')
set(gca,'xLim',[-75 225],'yLim',[-2100 500])
%legend('Speed < 0.5 cm/sec', '','Speed = 5 cm/sec','','Speed > 15 cm/sec','','50ms by 200uV')
% NOTE: scale bar is 50ms in x vs. 200uV in y

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/data13/mcarr/VelocityPaper/Ripples/%d_%d_%d_speedvsrip_trace_examples.ps', m, d, y);
print('-dpsc', savestring)
savestring = sprintf('/home/mcarr/Figures/%d_%d_%d_speedvsrip_trace_examples.jpg', m, d, y);
print('-djpeg',savestring)

% Save position figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/%d_%d_%d_speedvsrip_position_examples.pdf', m, d, y);
print('-dpdf',savestring)
