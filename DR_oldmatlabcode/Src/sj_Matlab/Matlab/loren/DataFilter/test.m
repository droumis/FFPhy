% plot the ripple times



load /data/loren/mkarlsso/Fra/fraripples08


r = ripples{8}{1}{1}.data(:,1:2);

hold on
for i = 1:length(r)
    plot([r(i,1) r(i,2)], [0 0 ], 'b');
end
s = f(1).output{8}{1};
for t = 1:length(s)
    for c = 1:length(s{t})
	if ~isempty(s{t}{c})
	    plot(s{t}{c}.data(:,1), ...
	    	zeros(length(s{t}{c}.data(:,1)), 1), '.');
	end
    end
end

