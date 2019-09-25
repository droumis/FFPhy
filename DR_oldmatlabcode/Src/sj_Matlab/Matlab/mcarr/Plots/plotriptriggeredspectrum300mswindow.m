%% Plot examples with 300ms time windows
directoryname = '/data21/mcarr/Ten/';
fileprefix = 'ten';
day = 4;
epoch = 6;

params = {};
params.Fs = 1500;
params.fpass = [2 350];
params.trialave = 0;
win = [1 1];
cwin = [0.3 0.01];
cellfilter = 'isequal($area,''CA1'')';
tetfilter = 'isequal($area,''CA1'') | isequal($area,''CA3'')';

ripples = loaddatastruct(directoryname,fileprefix,'ripples',day);
cellinfo = loaddatastruct(directoryname,fileprefix,'tetinfo',day);
valid = evaluatefilter(cellinfo{day},tetfilter);
valid = unique(valid(:,2));

%Find valid riptimes
riptimes = getripples([day epoch], ripples, cellinfo, 'cellfilter', cellfilter,'minstd',3);
pos = loaddatastruct(directoryname,fileprefix,'pos',day);
speed = pos{day}{epoch}.data(lookup(riptimes(:,1),pos{day}{epoch}.data(:,1)),8);
minVelocity = 4;
riptimes = riptimes(speed < minVelocity,:); clear speed minVelocity

% Define triggering events as the start of each ripple
triggers = riptimes(:,1)-pos{day}{epoch}.data(1,1);

%Remove triggering events that are too close to the beginning or end
while triggers(1)<win(1)
    triggers(1) = [];
end
while triggers(end)> pos{day}{epoch}.data(end,1)-win(2)
    triggers(end) = [];
end
clear cellfilter tetfilter pos ripples

%Create list of tetrodes for this day
tmpflist = dir(sprintf('%s/EEGnonreference/*eeg%02d-%d-*.mat', directoryname, day,epoch));
flist = cell(size(tmpflist));
for i = 1:length(tmpflist)
    flist{i} = sprintf('%s/EEGnonreference/%s', directoryname, tmpflist(i).name);
end

ca1tetrodes = evaluatefilter(cellinfo{day}{epoch},'isequal($area,''CA1'') & $numcells > 1');
ca3tetrodes = evaluatefilter(cellinfo{day}{epoch},'isequal($area,''CA3'') & $numcells > 1');
spectrum = cell(length(flist),1);
%Go through each file in flist and compute the riptriggered spectrogram
for fnum = 1:length(flist)
    %get the terode number
    dash = find(flist{fnum} == '-');
    tet = str2double(flist{fnum}((dash(2)+1):dash(2)+3));

    %Determine if tetrode is valid
    if any(tet==valid)
        load(flist{fnum});
        e = eeg{day}{epoch}{tet}.data'; clear eeg

        % Calculate the event triggered spectrogram
        [S,t,f] = mtspecgramtrigc(e,triggers,[win(1) win(2)],[cwin(1) cwin(2)],params);

        % Compute a z-scored spectrogram using the mean and std for the entire session
        P = mtspecgramc(e,[cwin(1) cwin(1)],params); clear e
        meanP = mean(P);
        stdP = std(P);

        for i = 1:size(S,1)
            for j = 1:size(S,3)
                S(i,:,j) = (S(i,:,j) - meanP)./stdP;
            end
        end
        
        spectrum{tet}.spectrum = S; clear S meanP stdP P
        spectrum{tet}.time = t-win(1);
        spectrum{tet}.frequency = f;
    end        
end
clear cellinfo dash directoryname tet triggers tmpflist valid

tmp = [];
for i = ca1tetrodes'
    if isempty(tmp)
        tmp = spectrum{i}.spectrum;
        count = 1;
        time = spectrum{i}.time;
        freq = spectrum{i}.frequency;
    else
        tmp = tmp+spectrum{i}.spectrum;
        count = count+1;
    end
end
CA1spectrum = tmp./count;
tmp = [];
for i = ca3tetrodes'
    if isempty(tmp)
        tmp = spectrum{i}.spectrum;
        count = 1;
        time = spectrum{i}.time;
        freq = spectrum{i}.frequency;
    else
        tmp = tmp+spectrum{i}.spectrum;
        count = count+1;
    end
end
CA3spectrum = tmp./count;
clear ca1tetrodes ca3tetrodes count day endtime tmp

%CA1 spectrum
figure
surf(time,freq,mean(CA1spectrum,3)')
view(0,90); shading interp; colormap hot
set(gca,'clim',[-0.5 2],'xlim',[time(46) time(126)],'xtick',time(46:20:126),'xticklabel',-0.4:0.2:0.4,'ylim',freq([1 end]))
colorbar('ytick',-0.5:0.5:2)

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_ca1spectrum_300mswindows.png', m, d, y);
print('-dpng', savestring)


%CA3 spectrum
figure
surf(time,freq,mean(CA3spectrum,3)')
view(0,90); shading interp; colormap hot
set(gca,'clim',[-0.5 1.5],'xlim',[time(46) time(126)],'xtick',time(46:20:126),'xticklabel',-0.4:0.2:0.4,'ylim',freq([1 end]))
colorbar('ytick',-0.5:0.5:2)

%Save Figure
[y,m,d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_ca3spectrum_300mswindows.png', m, d, y);
print('-dpng', savestring)
