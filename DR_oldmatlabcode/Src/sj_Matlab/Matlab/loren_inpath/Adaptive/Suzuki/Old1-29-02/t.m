figure(1);
figure(2);
'move figure 2'
pause
lastfive = zeros(5,1);
for i=1:589,
  if (testID(i)),
    [xvals tvals] = instantspline(CP.x, CP.t, thetainit, thetahat, sum(ntimesteps(1:i)-fixationtimes(i)));
    [i thetahat(1,(i+1))]
    for j = 1:4
       lastfive(j) = lastfive(j+1);
    end
    if (cobj.response(i) ~= 0 ) 
      lastfive(5) = 1;
      figure(1);
      plot(CP.x(1:end-2), xvals(1:end-2));
      %axis([0 2 0 45]);
      set(gca, 'XLim', [0 2]);
      set(gca, 'XTick', [.3 .8 1.5], 'XTickLabel', char('300', '800', '1500'));
      title(sprintf('%2f%% correct', sum(lastfive)/5 * 100));
      figure(2);
      semilogx(CP.t(1:end-2), tvals(1:end-2));
    else
      lastfive(5) = 0;
      figure(1);
      plot(CP.x(1:end-2), xvals(1:end-2), 'r');
      set(gca, 'XLim', [0 2]);
      %axis([0 2 0 45]);
      set(gca, 'XTick', [.3 .8 1.5], 'XTickLabel', char('300', '800', '1500'));
      title(sprintf('%2f%% correct', sum(lastfive)/5 * 100));
      figure(2);
      semilogx(CP.t(1:end-2), tvals(1:end-2), 'r');
    end;
  end;
end; 
