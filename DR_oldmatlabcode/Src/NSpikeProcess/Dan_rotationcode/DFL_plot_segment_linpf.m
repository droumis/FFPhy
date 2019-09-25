function plot_segment_linpf(a,b,bin_coord,normalized_firing_rate)
color_map = hsv(20);
if b(1)-a(1) < 0
    a_temp = a;
    a = b;
    b = a_temp;
    normalized_firing_rate = normalized_firing_rate(end:-1:1);
end
hold off;
line([a(1) b(1)],[a(2) b(2)]);
hold on;
pt_prev = a;
for ii = 2:length(bin_coord)
    (b(2)-a(2))/(b(1)-a(1))
    pt_next = [bin_coord(ii-1)*cos(atan((b(2)-a(2))/(b(1)-a(1)))) + a(1), bin_coord(ii-1)*sin(atan((b(2)-a(2))/(b(1)-a(1)))) + a(2)];
    line([pt_prev(1) pt_next(1)],[pt_prev(2) pt_next(2)],'Color',ind2rgb(floor((normalized_firing_rate(ii-1)))+1,jet(20)),'LineWidth',3);
    pt_prev = pt_next;
end
