function showlinmeanvel(num)

global fmaux linmeanvel

load([fmaux.data2dir '/linmeanvel']);

if nargin < 1
    [d,e,t,c]= startCellList;
    while ~isempty(d)
	auxrun(d,e);
	[d,e]= getNextEpoch;
    end
else
    auxrun(num(1),num(2));
end

function auxrun(d, e)

global fmaux linmeanvel

figure;
subplot(4,1,1)
plot(linmeanvel{d}{e}.data{1}(:,1), linmeanvel{d}{e}.data{1}(:,2));

title(sprintf('%s, day %d, epoch %d', fmaux.prefix, d, e));

subplot(4,1,2)
plot(linmeanvel{d}{e}.data{2}(:,1), linmeanvel{d}{e}.data{2}(:,2));
subplot(4,1,3)
plot(linmeanvel{d}{e}.data{3}(:,1), linmeanvel{d}{e}.data{3}(:,2));
subplot(4,1,4)
plot(linmeanvel{d}{e}.data{4}(:,1), linmeanvel{d}{e}.data{4}(:,2));
