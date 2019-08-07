


%{
- label/color ntrodes, organize by area/layer
- filter ripples for only the included ones

%}

savefigs = 0;
pausefigs = 1;
%%
animal = 'D10';
animdef = animaldef(animal);
days = 4;
epochs = 2;
eeg = loadeegstruct(animdef{2}, animal, 'eeg', days, epochs, 1:30);
% animaldir, animalprefix, datatype, days, epochs, tetrodes
rips = loaddatastruct(animdef{2}, animal, 'ca1rippleskons');
pos = loaddatastruct(animdef{2}, animal, 'pos');
% plot lfp, rips, pos, speed
%
day = 4;
epoch = 2;
eegstack = cell2mat(cellfun(@(x) x.data', eeg{day}{epoch}, 'un', 0)');
srate = eeg{day}{epoch}{1}.samprate;
starttime = eeg{day}{epoch}{1}.starttime;
endtime = eeg{day}{epoch}{1}.endtime;
nsamps = length(eegstack(1,:));
time = repmat((starttime+(1:nsamps)/srate)',1,length(eegstack(:,1)));
gridmat = 600*repmat((1:length(eegstack(:,1)))', 1,length(eegstack(1,:)));
eeggrid = (eegstack + gridmat)';
%%
Pp=load_plotting_params({'defaults','dataExplore'});
% ---- init fig----
if savefigs && ~pausefigs
    close all
    ifig = figure('Visible','off','units','normalized','position', ...
        Pp.position);
else
    ifig = figure('units','normalized','position',Pp.position);
end
set(gcf,'color','white')

startT = 1;
endT = 600;
startIdx = knnsearch(time(:,1),startT+starttime);
endIdx = knnsearch(time(:,1),endT+starttime);

posd = pos{day}{epoch}.data(:,[6 7 8 9]); % x y h v
postime = pos{day}{epoch}.data(:,1);
postime = repmat(postime,1,length(posd(1,:)));
posstartIdx = knnsearch(postime(:,1),startT+starttime);
posendIdx = knnsearch(postime(:,1),endT+starttime);

% plot rips in win
winstarttime = postime(posstartIdx);
winendtime = postime(posendIdx);
ripindsinwin = find(rips{day}{epoch}{1}.starttime>winstarttime & ...
    rips{day}{epoch}{1}.endtime<winendtime);
ripsinwinTimes = [rips{day}{epoch}{1}.starttime(ripindsinwin) ...
    rips{day}{epoch}{1}.endtime(ripindsinwin)];
xs = ripsinwinTimes(:,1); % -lfpstack(ian).ripStartTime(irip);
xe = ripsinwinTimes(:,2); %-lfpstack(ian).ripStartTime(irip);
yl = ylim;

sfl = subaxis(6,1,[1:3], 'SpacingVert', Pp.SpVt, 'SpacingHoriz', Pp.SpHz, 'MarginLeft', ...
    Pp.MgLt, 'MarginRight', Pp.MgRt, 'MarginTop', Pp.MgTp, 'MarginBottom', Pp.MgBm);
plot(time(startIdx:endIdx,:), eeggrid(startIdx:endIdx,:), 'LineWidth', 1, 'Color', [.2 .2 .2])
hold on
h = zoom;
h.Motion = 'horizontal';
h.Enable = 'on';
axis tight
set(gca, 'xtick', []);
ylabel('lfp')
p = pan;
p.Motion = 'horizontal';
patch([xs'; xe'; xe'; xs'], repmat([yl(1) yl(1) yl(2) yl(2)]', 1, length(xe)),'y', ...
    'FaceAlpha',.15, 'edgecolor','none');
hold off

sfp = subaxis(6,1,4, 'SpacingVert', Pp.SpVt, 'SpacingHoriz', Pp.SpHz, 'MarginLeft', ...
    Pp.MgLt, 'MarginRight', Pp.MgRt, 'MarginTop', Pp.MgTp, 'MarginBottom', Pp.MgBm);
plot(postime(posstartIdx:posendIdx,[1 2]), posd(posstartIdx:posendIdx,[1 2]), ...
    'LineWidth', 2)
hold on
axis tight
set(gca, 'xtick', []);
ylabel('xypos')
yl = ylim;
patch([xs'; xe'; xe'; xs'], repmat([yl(1) yl(1) yl(2) yl(2)]', 1, length(xe)),'y', ...
    'FaceAlpha',.15, 'edgecolor','none');
hold off

sfp = subaxis(6,1,5,'SpacingVert', Pp.SpVt, 'SpacingHoriz', Pp.SpHz, 'MarginLeft', ...
    Pp.MgLt, 'MarginRight', Pp.MgRt, 'MarginTop', Pp.MgTp, 'MarginBottom', Pp.MgBm);
plot(postime(posstartIdx:posendIdx,1), posd(posstartIdx:posendIdx,[4]), 'LineWidth', 3, ...
    'Color', [.1 .4 .4])
hold on
axis tight
set(gca, 'xtick', []);
ylabel('velocity')
yl = ylim;
patch([xs'; xe'; xe'; xs'], repmat([yl(1) yl(1) yl(2) yl(2)]', 1, length(xe)),'y', ...
    'FaceAlpha',.15, 'edgecolor','none');
hold off

sfp = subaxis(6,1,6,'SpacingVert', Pp.SpVt, 'SpacingHoriz', Pp.SpHz, 'MarginLeft', ...
    Pp.MgLt, 'MarginRight', Pp.MgRt, 'MarginTop', Pp.MgTp, 'MarginBottom', Pp.MgBm);
plot(postime(posstartIdx:posendIdx,1), posd(posstartIdx:posendIdx,[3]), 'LineWidth', 2, ...
    'Color', [.5 .5 .5], 'LineStyle', ':', 'MarkerSize', .5)
hold on
axis tight
ylabel('head dir')
xlabel('time s')
yl = ylim;
patch([xs'; xe'; xe'; xs'], repmat([yl(1) yl(1) yl(2) yl(2)]', 1, length(xe)),'y', ...
    'FaceAlpha',.15, 'edgecolor','none');
hold off

allAxesInFigure = findall(ifig,'type','axes');
linkaxes(allAxesInFigure, 'x');

xlim([starttime, starttime+20])


