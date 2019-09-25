load('/data14/jai/L1_/L1data01.mat');

rawpos = data{1,1}{1,2}.Pos.rawpos;
corrected = data{1,1}{1,2}.Pos.correcteddata;
figure;
plot(rawpos(:, 2), rawpos(:, 3),'k');

hold on
plot(corrected(:, 2).*2, corrected(:, 3).*2,'r')

