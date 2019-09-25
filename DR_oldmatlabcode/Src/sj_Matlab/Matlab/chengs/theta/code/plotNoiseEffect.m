% a=[5 10 20 30 60]';
a=[5 10 15 20 30]';
plotcol= hsv(10);
plotcol(1:6,:)=[[0 0 0]; [1 0 0]; [0 0 1]; [.7 0 1]; [0 1 0]; [0 1 1]];
x= [0:1:125];
clf
for i=1:length(a); 
    y= sqrt(x.^2+a(i)^2)-x;
%     h= plot(x,y/a(i));
    h= plot(x,y);
    set(h, 'Color', plotcol(i+2,:));
    hold on
end
axis tight
plot(get(gca,'XLim'), [2 2], 'k--')
xlabel('"real" std (ms)')
ylabel('\Delta std (ms)')
legend(num2str(a), 0)
myprint('large', 'jitterEffect')
