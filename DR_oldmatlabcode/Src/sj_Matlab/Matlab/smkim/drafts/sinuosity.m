
% count various stuffs for learning curves
labels = { 'M02' 'M03' 'M06' 'M12' 'M13' 'M14' 'M16' 'M17' 'M19' 'M20' 'M22' 'M24' 'M25' 'M26' };
islesion = [  1     1     0     1     1     1     1     1     1     1     1     0     0     0    ];

stepsizes = (1:1:100)';
vals = zeros(size(stepsizes));

%
figure(1);
set(gca,'YLim',[0 1]);
for s = 1:numel(labels)
    load([labels{s} '_Wtrack_smoothedpos.mat']);
    for i = 1:numel(stepsizes)
        vals(i) = pathsinuosity(smoothedpos{1}{1},stepsizes(i));
    end
    if islesion(s)
        line(stepsizes,vals,'Color','m','Marker','o');        
    else
        line(stepsizes,vals,'Color','g','Marker','o');        
    end
end
%

%{
figure(2);
hold on;
for s = 1:numel(labels)
    load([labels{s} '_Wtrack_smoothedpos.mat']);
    speed = hypot(smoothedpos{1}{1}.data(:,4),smoothedpos{1}{1}.data(:,5));
    [f,x] = ecdf(speed);
    if islesion(s)
        line(x,f,'Color','m','Marker','none','LineWidth',2);        
    else
        line(x,f,'Color','g','Marker','none','LineWidth',2);        
    end
end
%}

