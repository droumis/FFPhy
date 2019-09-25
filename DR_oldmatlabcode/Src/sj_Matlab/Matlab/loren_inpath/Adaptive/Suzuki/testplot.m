for i = 12:14
    trials = dataset{30}.condtrial{i};
    f = dataset{30}.cell{1}.fixrate(trials);
    s = dataset{30}.cell{1}.stimrate(trials);
    d = dataset{30}.cell{1}.delayrate(trials);
    r = dataset{30}.cell{1}.resprate(trials);
    plot(f, 'k');
    hold on
    plot(s, 'r');
    plot(d, 'g');
    plot(r, 'b');
    legend('fixation', 'stimulus', 'delay', 'response');
    pause
    clf
end

