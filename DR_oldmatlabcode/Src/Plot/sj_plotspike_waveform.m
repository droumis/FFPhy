function h = sj_plotspike_waveform(daydirect,prefix,day,tet,cell,usech)

% Shantanu - May 2012
% Using matclust file and tt file, plot waveform for given cell

% eg. sj_plotspike_waveform('/data25/sjadhav/HPExpt/JW7/','JW7',5,8,3,[2:4]);
% sj_plotspike_waveform('/data25/sjadhav/HPExpt/HPa/','HPa',2,1,1,[1:4]);

if nargin<3,
    keyboard
    error('Please enter Folder, Expt Prefix and Day No!');
end
if nargin<4,
    tet=1;
end
if nargin<5,
    cell=1;
end
if nargin<6 || isempty(usech),
    usech=1:4;
end

set(0,'defaultaxesfontsize',24);
tfont = 28;
xfont = 20;
yfont = 20;
    
% Go to directory and get files
% ----------------------------
currdir = pwd;
cd(daydirect);
dayfolders = dir;
daystr = sprintf('%02d', day);
tetstr = sprintf('%02d', tet);
% Go to correct folder
for i = 3:length(dayfolders)
    % Go to Day Folder
    if dayfolders(i).isdir
        if strcmp(daystr,dayfolders(i).name(end-1:end))
            disp(upper(dayfolders(i).name))
            cd(dayfolders(i).name);
            % Go to Tet Folder
            tetfolders = dir;
            for j = 3:length(tetfolders)
                if tetfolders(j).isdir
                    if strcmp(tetstr,tetfolders(j).name(1:2))
                        disp(upper(tetfolders(j).name))
                        cd(tetfolders(j).name);
                    end
                end
            end
        end
    end
end

% Get Files and load raw and clust files
% --------------------------------------
% Raw File
% --------
ttfile = dir('*.tt');
if ~isempty(ttfile)
    disp(['      ',ttfile(1).name])
    dashfind = strfind(ttfile(1).name,'-');
    tetstring = ttfile(1).name(1:dashfind-1);
    [timestamps, waves] = readtt(ttfile(1).name);
end
% Clust File
% ----------
matclustfile = dir('matclust*');
if ~isempty(matclustfile)
    disp(['      ',matclustfile(1).name])
    load(matclustfile(1).name);
else
    disp('Error! There aint no clustered data');
    return
end

% Get tetrode info
% ------------------
currdir_tet = pwd;
cd ..
tetrodeinfo = getTetrodeConfig2;
thresh = tetrodeinfo(tet).thresh;
triggers = thresh;
cd(currdir_tet);

% Will index match the clustdata file - See if clustdata has some discarded spikes. You have to match indexes
% ----------------------------------------------------------------------------------------------------------
nspk_clust = size(clustdata.params,1);
% If not, get params again and see which ones were discarded by goodspikes
% THERE MUST BE A BETTER WAY. MATCLUST FILE SHOULD HAVE THAT INFORMATION
% ------------------------------------------------------------------------
if nspk_clust~=length(timestamps);
    MAXALLOWEDAMP = 2000; THRESH = 0;
    params = parmcalc(timestamps,waves,triggers)';
    goodspikes =  ((params(:,1)<MAXALLOWEDAMP)&(params(:,2)<MAXALLOWEDAMP)&(params(:,3)<MAXALLOWEDAMP)&(params(:,4)<MAXALLOWEDAMP)) & ...
        ((params(:,1)>THRESH)|(params(:,2)>THRESH)|(params(:,3)>THRESH)|(params(:,4)>THRESH)) ;
    waves=waves(:,:,goodspikes); % This will makes waves and clustdata indexes the same
end

% Get Spk Indexes from matclust file and get corresponding waveforms
% ------------------------------------------------------------------
spkidxs = double(clustattrib.clusters{cell}.index);
spkwaves = double(waves(:,:,spkidxs));

% % De-jitter waveforms using addpc code
% % ------------------------------------
[~, spkwaves_align] = sj_addpc(spkwaves, triggers, 0);

% Make a long vector out of the channels to be used for plotting
% --------------------------------------------------------------
spkwaves = spkwaves(:,usech,:);
spkwaves = reshape(spkwaves,size(spkwaves,1)*size(spkwaves,2),size(spkwaves,3));
spkwaves_align = spkwaves_align(:,usech,:);
spkwaves_align = reshape(spkwaves_align,size(spkwaves_align,1)*size(spkwaves_align,2),size(spkwaves_align,3));



% Plot subset
% -----

if length(spkwaves)>=50000
    randvec = randperm(50000);
    spkwaves_sub=spkwaves(rand,:,:);
    spkwaves_align_sub=spkwaves(rand,:,:);
end



figure; hold on;
%redimscreen_figforppt1;
plot(spkwaves,'b-');
plot(mean(spkwaves,2),'k-','LineWidth',3);
title([prefix,': Day',num2str(day),' Tet',num2str(tet),' Cell',num2str(cell)]);

%Plot - dejittered
%---------------------
figure; hold on;
%redimscreen_figforppt1;
if length(spkwaves)>=50000
    plot(spkwaves_align_sub,'-','Color',[0.5 0.5 0.5]);
else
    plot(spkwaves_align,'-','Color',[0.5 0.5 0.5]);
end
plot(mean(spkwaves_align,2),'r-','LineWidth',3);
title([prefix,': Day',num2str(day),' Tet',num2str(tet),' Cell',num2str(cell),' - Dejittered']);


% Density Plot
% figure; hold on;
% redimscreen_figforppt1;
% [n,x,y] = hist2d(spkwaves_align');
% imagesc(x,y,n); axis xy; colormap hot;
% title([prefix,': Day',num2str(day),' Tet',num2str(tet),' Cell',num2str(cell),' - Dejittered']);


% Density Plot
% figure; hold on;
% %redimscreen_figforppt1;
% [n,x,y] = hist2d(spkwaves');
% imagesc(x,y,n); axis xy; colormap hot;
% title([prefix,': Day',num2str(day),' Tet',num2str(tet),' Cell',num2str(cell) ],'FontSize',20,'Fontweight','normal');
% 
% 
% set(gca,'XLim',[0 160]); set(gca,'YLim',[-110 150]);
% 
% 
% % figdir = '/data25/sjadhav/HPExpt/HPa_direct/Figures/SpkWaves/'; prefix='HPa';
% % saveg=0;
% % if saveg==1
% %     figfile = [figdir,prefix,'PFCt17c6_SpkWave'];
% %     print('-dpdf', figfile); print('-djpeg', figfile); saveas(gcf,figfile,'fig');
% % end
% 









