figure
e1 = [288/ 612];
r1a = [256/ 580];
r1q = [373/ 1046];
e1ine2 = [182/ 442];
e2ine2 = [147/ 330];
subplot(1,2,1)
h = bar([e1 r1a r1q] * 100, 'k');
ylabel('Percent Significant');
set(gca, 'YLim', [0 50]);
xlim = get(gca, 'XLim');
set(gca, 'XTickLabel', char('E1', 'R a', 'R q'));

box off

subplot(1,2,2)
h = bar([e1ine2 e2ine2] * 100, 'k');
set(gca, 'YLim', [0 50], 'XLim', xlim);
set(gca, 'XTickLabel', char('E1', 'E2'));
box off
